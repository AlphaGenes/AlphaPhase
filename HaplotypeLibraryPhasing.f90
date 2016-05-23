module HaplotypeLibraryPhasing
  implicit none
  
contains
  subroutine MakeHapLib(library, c)
    use HaplotypeLibraryDefinition
    use CoreDefinition
    use Random
    use Parameters
    implicit none

    type(HaplotypeLibrary), intent(in) :: library
    type(Core) :: c

    integer :: i
    integer :: id

    call library%resetHapFreq()
    call c%resetFullyPhased()
    call c%resetHapAnis()

    !THIS IS HORRIBLE!
    call library%initalise(c%getNCoreSnp(),500,500)

    do i = 1, c%getNAnisG()
      !Paternal Haps
      if (fullyPhased(c%getHaplotype(i,1))) then
	id = library%matchAddHap(c%getHaplotype(i,1))
	call c%setHapAnis(i,1,id)
	call c%setFullyPhased(i,1)
      endif

      !Maternal Haps
      if (fullyPhased(c%getHaplotype(i,2))) then
	!FUDGE FOR CONSISTENCY.  There should be no if statement here
	if (.not. consistent .or. (library%getSize() > 0)) then
	  call c%setFullyPhased(i,2)
	end if
	id = library%matchAddHap(c%getHaplotype(i,2))
	call c%setHapAnis(i,2,id)
      endif
    enddo

  end subroutine MakeHapLib

  subroutine ImputeFromLib(library, c, nGlobalHapsIter)
    ! Impute the phase for gametes that are not completely phased by LRP 
    ! by matching their phased loci to haplotypes in the Haplotype Library,
    ! following strategies listed in the section Step 2e of Hickey et al 2011.

    use Parameters, only : percgenohaplodisagree, consistent
    use HaplotypeLibraryDefinition
    use Constants
    use CoreDefinition
    use Clustering
    use Random
    implicit none

    type(HaplotypeLibrary), intent(in) :: library
    type(Core) :: c
    integer, intent(inout) :: nGlobalHapsIter

    integer(kind=1), dimension(:,:), pointer :: genos
    integer :: i, j, k, HapLibIter, nHapsOld, ErrorAllow, HapM, HapP, nCand, nCandPat, nCandMat
    integer :: WorkScaler, CountZero, CountOne
    integer(kind=1) :: first, val
    integer :: id
    integer, allocatable, dimension(:) :: CandGenos, WorkVec
    integer, pointer, dimension(:,:) :: CandPairs

    integer :: nAnisG, nSnp, nHaps

    integer, dimension(:), pointer :: compatHaps

    integer, allocatable, dimension (:,:) :: TempHapArray
    integer, allocatable, dimension (:) :: TempHapVector, ClusterMember
    integer :: nHapsCluster, rounds

    integer :: countA, countB, ErrorCountAB

    integer(kind=1), pointer, dimension(:) :: comp
    integer, pointer, dimension(:) :: CandHapsPat, CandHapsMat, CandHaps, matches

    logical :: singlePat, singleMat

    !!!! TEMP STUFF
    integer, pointer, dimension(:) :: CandHaps2

    nAnisG = c%getNAnisG()
    nSnp = c%getNCoreSnp()
    nHaps = library%getSize()


    allocate(CandGenos(nSnp))
    allocate(CandHaps(nAnisG * 2))
    allocate(WorkVec(nAnisG * 2))
    allocate(CandPairs(nAnisG * 2, 2))

    genos => c%getCoreGenos()

    ErrorAllow = int(PercGenoHaploDisagree * nSnp)
    ErrorCountAB = int(nSnp * 0.09)

    HapLibIter = 1
    if (nGlobalHapsIter == 1) then
      print*, "   ", "Iteration ", nGlobalHapsIter, "found ", nHaps, "haplotypes"
      nHapsOld = 0
    else
      nHapsOld = 0
    endif

    do while (nHapsOld /= nHaps)
      HapLibIter = HapLibIter + 1
      nGlobalHapsIter = nGlobalHapsIter + 1
      nHapsOld = nHaps

      do i = 1, nAnisG
	CandHaps = 0
	nCand = 0
	ErrorAllow = int(PercGenoHaploDisagree * count(Genos(i, :) /= MissingGenotypeCode))

	if ((.not. c%getFullyPhased(i,1)) .or. (.not. c%getFullyPhased(i,2)))  then
	  ! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): PATERNAL HAPLOTYPE
	  if (c%getFullyPhased(i, 1)) then
	    comp => complement(Genos(i,:), c%getHaplotype(i,1))
	    CandHaps2 => library%matchWithError(comp, ErrorAllow)

	    if (size(CandHaps2,1) > 1) then
	      do j = 1, nSNP
		if (agree(CandHaps2, library, j)) then
		  call c%setPhase(i, j, 2, library%getPhase(CandHaps2(1), j))
		end if
	      end do
	    endif

	    if (size(CandHaps2,1) == 1) then
	      call c%setHaplotype(i,2,library%getHap(CandHaps2(1)))
	      call c%setFullyPhased(i,2)
	      call library%incrementHapFreq(CandHaps2(1))
	      call c%setHapAnis(i, 2, CandHaps2(1))
	    end if

	    if (size(CandHaps2,1) == 0) then
	      do j = 1, nSnp
		  if ((comp(j) == 0).or.(comp(j) == 1)) then
		    call c%setPhase(i, j, 2, comp(j))
		  endif
	      enddo

	      if (fullyPhased(comp)) then
		id = library%matchAddHap(c%getHaplotype(i,2))
		call c%setHapAnis(i,2,id)
		call c%setFullyPhased(i,2)

		!!! TEMP HACK as nHaps is still used lower
		nHaps = library%getSize()
	      end if
	    endif
	  end if

	  ! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): MATERNAL HAPLOTYPE
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Haplotype 2 can get fully phased above and this will run despite both haplotypes now being phased
	  ! Affects results - not entirely sure why...
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if (c%getFullyPhased(i, 2)) then
	    comp => complement(Genos(i,:), c%getHaplotype(i,2))
	    CandHaps2 => library%matchWithError(comp, ErrorAllow)

	    if (size(CandHaps2,1) > 1) then
	      do j = 1, nSNP
		if (agree(CandHaps2, library, j)) then
		  call c%setPhase(i, j, 1, library%getPhase(CandHaps2(1), j))
		end if
	      end do
	    endif

	    if (size(CandHaps2,1) == 1) then
	      call c%setHaplotype(i,1,library%getHap(CandHaps2(1)))
	      call c%setFullyPhased(i,1)
	      call library%incrementHapFreq(CandHaps2(1))
	      call c%setHapAnis(i, 1, CandHaps2(1))
	    end if

	    if (size(CandHaps2,1) == 0) then
	      do j = 1, nSnp
		  if ((comp(j) == 0).or.(comp(j) == 1)) then
		    call c%setPhase(i, j, 1, comp(j))
		  endif
	      enddo

	      if (fullyPhased(comp)) then
		id = library%matchAddHap(c%getHaplotype(i,1))
		call c%setHapAnis(i,1,id)
		call c%setFullyPhased(i,1)

		!!! TEMP HACK as nHaps is still used lower
		nHaps = library%getSize()
	      end if
	    endif
	  end if

	  ! If neither of the gametes is completely phased (Section Step 2e.ii Hickey et al. 2011)
	  if ((.not. c%getFullyPhased(i,1)) .and. (.not. c%getFullyPhased(i,2))) then
	    HapP = 0
	    HapM = 0
	    CandHaps = 0
	    nCand = 0	  

	    if (consistent) then
	      allocate(compatHaps(nHaps))
	      do k = 1, nHaps
		compatHaps(k) = k
	      end do
	    else
	      compatHaps => library%getCompatHaps(genos(i,:))
	    end if

	    CandHapsPat => library%limitedMatchWithError(c%getHaplotype(i,1), ErrorAllow, compatHaps)

	    ! Update the number of candidates for paternal haplotype
	    nCandPat = size(CandHapsPat,1)

	    ! If only have one paternal candidate haplotype, then
	    ! the paternal haplotype is nCand
	    if (nCandPat == 1) HapP = CandHapsPat(1)

	    !Find all candiate maternal haps and then remove those that are already candidate paternal haps
	    CandHapsMat => uniqueHaps( &
	      library%limitedMatchWithError(c%getHaplotype(i,2), ErrorAllow, compatHaps), CandHapsPat)

	    ! If only have one maternal candidate haplotype, then
	    ! the maternal haplotype is nCand
	    nCandMat = size(CandHapsMat,1)
	    if (nCandMat == 1) HapM = CandHapsMat(1)

	    ! Make array of all haps
	    nCand = size(CandHapsPat) + size(CandHapsMat)
	    allocate(CandHaps(nCand))
	    CandHaps(1:size(CandHapsPat)) = CandHapsPat
	    CandHaps(size(CandHapsPat)+1:nCand) = CandHapsMat

	    ! If only one maternal candidate haplotype and many paternal candidate haplotypes
	    !! There's some odd logic here - if we have one paternal and maternal we end up keeping HapP as the
	    !! paternal even if it's not compitable with the mat.  May be checked later.
	    if ((nCandMat == 1).AND.(nCandPat > 0)) then
	      comp => complement(Genos(i,:),library%getHap(HapM))
	      matches => library%limitedMatchWithError(comp, ErrorAllow, CandHapsPat)
	      if (size(matches) == 1) then
		HapP = matches(1)
	      end if
	    end if

	    deallocate(compatHaps)

	    ! If only one paternal candidate haplotype and one / many maternal candidate haplotypes
	    if ((nCandPat == 1).and.(nCandMat > 0)) then
	      ! FUDGE TO DEAL WITH STRANGE LOGIC IN ORIGINAL CODE
	      HapM = 0

	      comp => complement(Genos(i,:),library%getHap(HapP))
	      matches => library%limitedMatchWithError(comp, ErrorAllow, CandHapsMat)
	      if (size(matches) == 1) then
		HapM = matches(1)
	      end if

	      ! FUDGE TO DEAL WITH STRANGE LOGIC IN ORIGINAL CODE
	      if (HapM /= CandHaps(nCand)) then
		HapM = 0
	      end if
	    end if

	  ! If only have one paternal candidate haplotype
	  if (HapP /= 0) then
	    call c%setHaplotype(i,1,library%getHap(HapP))
	    call c%setFullyPhased(i, 1)

	    ! Update the Library
	    call library%incrementHapFreq(HapP)
	    call c%setHapAnis(i, 1, HapP)

	    ! If no haplotype has been found for the maternal gamete, or 
	    ! there are more than one maternal candidate haplotype
	    if (HapM == 0) then
	      comp => complement(genos(i,:), c%getHaplotype(i,1))

	      do j = 1, nSnp
		  if ((comp(j) == 0).or.(comp(j) == 1)) then
		    call c%setPhase(i, j, 2, comp(j))
		  endif
	      enddo

	      if (fullyPhased(comp)) then
		id = library%matchAddHap(c%getHaplotype(i,2))
		call c%setHapAnis(i,2,id)
		call c%setFullyPhased(i,2)

		!!! TEMP HACK as nHaps is still used lower
		nHaps = library%getSize()
	      end if	    
	    end if
	  end if

	  ! If only have one maternal candidate haplotype 
	  if (HapM /= 0) then
	    ! Phase the maternal haplotype with the only paternal candidate haplotype
	    ! NOTE: This is always necessary because
	    !   - If there is only one paternal candidate haplotype, 
	    !     we have already check they are compatible (line 3307-3316) (Step 2e.ii.A)
	    !   - If there is more than one paternal candidate haplotype, 
	    !     the paternal gamete is phased as the complementary of the maternal gamete
	    call c%setHaplotype(i, 2, library%getHap(HapM))
	    call c%setFullyPhased(i, 2)

	    ! Update the Library
	    call library%incrementHapFreq(HapM)
	    call c%setHapAnis(i, 2, HapM)

	    ! If no haplotype has been found for the paternal gamete, or 
	    ! there are more than one paternal candidate haplotype
	    if (HapP == 0) then
	      comp => complement(genos(i,:), c%getHaplotype(i,2))

	      do j = 1, nSnp
		  if ((comp(j) == 0).or.(comp(j) == 1)) then
		    call c%setPhase(i, j, 1, comp(j))
		  endif
	      enddo

	      if (fullyPhased(comp)) then
		id = library%matchAddHap(c%getHaplotype(i,1))
		call c%setHapAnis(i,1,id)
		call c%setFullyPhased(i,1)

		!!! TEMP HACK as nHaps is still used lower
		nHaps = library%getSize()
	      end if	
	    end if
	  end if

	  ! If the paternal and maternal gamete cannot be identify without ambiguity
	  ! (more than one or none at all)           
	  if ((HapP == 0).and.(HapM == 0)) then
	    CandPairs => library%limitedCompatPairsWithError(Genos(i,:),ErrorAllow,CandHaps,nAnisG)

	    ! If only one pair agrees...
	    if (size(CandPairs,1) == 1) then
	      ! Phase the paternal haplotype and update the library with the new frequency 
	      HapP = CandPairs(1,1)
	      call c%setHaplotype(i, 1, library%getHap(HapP))
	      call c%setFullyPhased(i, 1)
	      call library%incrementHapFreq(HapP)
	      call c%setHapAnis(i, 1, HapP)

	      ! Phase the maternal haplotype and update the library with the new frequency 
	      HapM = CandPairs(1,2)
	      call c%setHaplotype(i, 2, library%getHap(HapM))
	      call c%setFullyPhased(i, 2)
	      call library%incrementHapFreq(HapM)
	      call c%setHapAnis(i, 2, HapM)
	    end if

	    ! If more than one pair agrees...                
	    if ((size(CandPairs,1) > 1).and.((size(CandPairs,1) * size(CandPairs,1)) < nAnisG)) then !Note the 200 number is a fudge

	      ! Check how many paternal candidates haplotypes
	      singlePat = .true.
	      id = CandPairs(1, 1)
	      do k = 2, Size(CandPairs,1)
		if (CandPairs(k, 1) /= id) then
		  SinglePat = .false.
		end if
	      end do

	      ! If there is only one paternal haplotype in all the candidate pairs
	      if (singlePat) then
		! Phase the paternal gamete with this haplotype
		call c%setHaplotype(i, 1, library%getHap(id))
		call c%setFullyPhased(i, 1)
		call library%incrementHapFreq(id)
		call c%setHapAnis(i, 1, id)

		! If only one haplotype is found for the paternal gamete 
		! and many for the maternal gamete, phase each loci only all pairs agree
		! (Step 2e.ii.B)
		do j = 1, nSNP
		  if (agree(CandPairs(:,2), library, j)) then
		    call c%setPhase(i, j, 2, library%getPhase(CandPairs(1, 2), j))
		  end if
		end do
	      end if

	      ! Check how many maternal candidates haplotypes
	      SingleMat = .true.
	      id = CandPairs(1, 2)
	      do k = 2, size(CandPairs,1)
		if (CandPairs(k, 2) /= id) then
		  SingleMat = .false.
		end if
	      end do

	      ! If there is only one maternal haplotype in all the candidate pairs
	      if (SingleMat) then
		call c%setHaplotype(i, 2, library%getHap(id))
		call c%setFullyPhased(i, 2)
		call library%incrementHapFreq(id)
		call c%setHapAnis(i, 2, id)

		! If only one haplotype is found for the paternal gamete 
		! and many for the maternal gamete, phase each loci only all pairs agree
		! (Step 2e.ii.C)
		do j = 1, nSNP
		  if (agree(CandPairs(:,1), library, j)) then
		    call c%setPhase(i, j, 1, library%getPhase(CandPairs(1,1), j))
		  endif
		enddo
	      endif

	      ! If proband is not completely phased and have more than one candidate 
	      ! for both paternal and maternal haplotype 
	      ! (Step 2e.iv)
	      if ( ((.not. c%getFullyPhased(i,1)) .or. (.not. (c%getFullyPhased(i,2)))) &
		.and. (.not. SinglePat) .and. (.not. SingleMat)) then

		! Initialize procedure of k-medoids
		WorkVec = 0
		do k = 1, size(CandPairs,1)
		  WorkVec(CandPairs(k, 1)) = 1
		  WorkVec(CandPairs(k, 2)) = 1
		enddo
		WorkScaler = sum(WorkVec(:))
		allocate(TempHapArray(WorkScaler, nSnp))
		allocate(TempHapVector(WorkScaler))
		nHapsCluster = 0

		! Clusterize with k-medoids
		! I think this is actually k-means!
		do k = 1, nAnisG * 2
		  if (WorkVec(k) == 1) then
		    nHapsCluster = nHapsCluster + 1
		    TempHapVector(nHapsCluster) = k
		    TempHapArray(nHapsCluster, 1:c%getNCoreSnp()) = &
		    library%getHap(k)
		  end if
		end do
		allocate(ClusterMember(nHapsCluster))
		do j = 1, nHapsCluster
		  if (mod(j, 2) == 0) then
		    ClusterMember(j) = 1
		  else
		    ClusterMember(j) = 2
		  endif
		end do
		rounds = cluster(TempHapArray, ClusterMember, 2, nMaxRounds, .false.)
		if (rounds <= nMaxRounds) then
		  if (count(ClusterMember(:) == 2) == 1) then
		    HapM = 0
		    do j = 1, nHapsCluster
		      if (ClusterMember(j) == 2) then
			HapM = TempHapVector(j)
			exit
		      endif
		    enddo
		    call c%setHaplotype(i,2, library%getHap(HapM))
		    call c%setFullyPhased(i, 2)
		    call library%incrementHapFreq(HapM)
		    call c%setHapAnis(i, 2, HapM)
		  end if
		  if (count(ClusterMember(:) == 1) == 1) then
		    HapP = 0
		    do j = 1, nHapsCluster
		      if (ClusterMember(j) == 1) then
			HapP = TempHapVector(j)
			exit
		      endif
		    enddo
		    call c%setHaplotype(i, 1, library%getHap(HapP))
		    call c%setFullyPhased(i, 1)
		    call library%incrementHapFreq(HapP)
		    call c%setHapAnis(i, 1, HapP)
		  end if
		  if ((count(ClusterMember(:) == 2) > 1).and.(count(ClusterMember(:) == 2) > 1)) then
		    call c%setHaplotypeToUnphased(i,1)
		    call c%setHaplotypeToUnphased(i,2)
		    do j = 1, nSnp
		      if (Genos(i, j) == 0) then
			call c%setPhase(i, j, 1, 0)
			call c%setPhase(i, j, 2, 0)
		      end if
		      if (Genos(i, j) == 2) then
			call c%setPhase(i, j, 1, 1)
			call c%setPhase(i, j, 2, 1)
		      end if
		      CountZero = 0
		      CountOne = 0
		      do k = 1, nHapsCluster
			if (ClusterMember(k) == 2) then
			  if (library%getPhase(TempHapVector(k), j) == 0)&
			  CountZero = CountZero + 1
			  if (library%getPhase(TempHapVector(k), j) == 1)&
			  CountOne = CountOne + 1
			endif
		      end do
		      if ((CountZero == 0).and.(CountOne > 0)) call c%setPhase(i, j, 2, 1)
		      if ((CountZero > 0).and.(CountOne == 0)) call c%setPhase(i, j, 2, 0)

		      CountZero = 0
		      CountOne = 0
		      do k = 1, nHapsCluster
			if (ClusterMember(k) == 1) then
			  if (library%getPhase(TempHapVector(k), j) == 0) CountZero = CountZero + 1                                                  
			  if (library%getPhase(TempHapVector(k), j) == 1) CountOne = CountOne + 1                                              
			endif
		      end do
		      if ((CountZero == 0).and.(CountOne > 0)) call c%setPhase(i, j, 1, 1)
		      if ((CountZero > 0).and.(CountOne == 0)) call c%setPhase(i, j, 1, 0)
		    end do
		  endif
		end if
		deallocate(ClusterMember)
		deallocate(TempHapArray)
		deallocate(TempHapVector)
	      endif
	    endif
	  endif
	endif
      end if   
    end do
    print*, "   ", "Iteration ", nGlobalHapsIter, "found ", nHaps, "haplotypes"
  enddo

  do i = 1, nAnisG
    CountA = 0
    CountB = 0
    do j = 1, nSnp
      if (Genos(i, j) /= MissingGenotypeCode) then
	if ((c%getPhase(i, j, 1) /= 9).and.(c%getPhase(i, j, 2) == 9)) then
	  val = Genos(i, j) - c%getPhase(i, j, 1)
	  if ((val == 0).or.(val == 1)) then !here 7th april 2011
	    call c%setPhase(i, j, 2, val)
	  else
	    CountA = CountA + 1
	  endif
	endif
	if ((c%getPhase(i, j, 2) /= 9).and.(c%getPhase(i, j, 1) == 9)) then
	  val = Genos(i, j) - c%getPhase(i, j, 2)
	  if ((val == 0).or.(val == 1)) then !here 7th april 2011
	    call c%setPhase(i, j, 1, val)
	  else
	    CountB = CountB + 1
	  endif
	endif
      end if
    end do
    if ((CountA > ErrorCountAB) .or. (CountB > ErrorCountAB)) then
      call c%setHaplotypeToUnphased(i,1)
      call c%setHaplotypeToUnphased(i,2)
      do j = 1, nSnp
	if (Genos(i, j) == 0) then
	  call c%setPhase(i, j, 1, 0)
	  call c%setPhase(i, j, 2, 0)
	end if
	if (Genos(i, j) == 2) then
	  call c%setPhase(i, j,1, 1)
	  call c%setPhase(i, j,1, 1)
	end if
      enddo
    endif
  end do

  do i = 1, nAnisG
    do j = 1, nSnp
      if (Genos(i, j) == 1) then
	if ((c%getPhase(i, j, 1) == 9).and.(c%getPhase(i, j, 2) /= 9)) call c%setPhase(i, j, 1, Genos(i, j) - c%getPhase(i, j, 2))
	if ((c%getPhase(i, j, 2) == 9).and.(c%getPhase(i, j, 1) /= 9)) call c%setPhase(i, j, 2, Genos(i, j) - c%getPhase(i, j, 1))

      endif
      if (Genos(i, j) == 0) then
	if ((c%getPhase(i, j, 1) == 9).and.(c%getPhase(i, j, 2) /= 9)) call c%setPhase(i, j, 1, 0)
	if ((c%getPhase(i, j, 2) == 9).and.(c%getPhase(i, j, 1) /= 9)) call c%setPhase(i, j, 2, 0)

      endif
      if (Genos(i, j) == 2) then
	if ((c%getPhase(i, j, 1) == 9).and.(c%getPhase(i, j, 2) /= 9)) call c%setPhase(i, j, 1, 1)
	if ((c%getPhase(i, j, 2) == 9).and.(c%getPhase(i, j, 1) /= 9)) call c%setPhase(i, j, 2, 1)
      endif
    enddo
  enddo

  call library%resetHapFreq()
  call c%resetFullyPhased()
  call c%resetHapAnis()

  deallocate(CandGenos)
  deallocate(CandHaps)
  deallocate(WorkVec)
  deallocate(CandPairs)

  end subroutine ImputeFromLib

  function agree(candidates, library, position) result (a)
    use HaplotypeLibraryDefinition
    
    integer, dimension(:), intent(in) :: candidates
    class(HaplotypeLibrary), intent(in) :: library
    integer, intent(in) :: position

    logical :: a
    integer(kind = 1), dimension(size(candidates,1)) :: phases
    integer :: i

    do i = 1, size(candidates,1)
      phases(i) = library%getPhase(candidates(i),position)
    end do

    a = (all(phases == 0) .and. (all(phases == 1)))
  end function agree
end module HaplotypeLibraryPhasing