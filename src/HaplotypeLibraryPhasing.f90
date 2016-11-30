module HaplotypeLibraryPhasing
  use Constants
  implicit none

  integer, parameter, private :: nMaxRounds = 100
  
contains
  subroutine MakeHapLib(c, library, consistent)
    use HaplotypeLibraryDefinition
    use CoreDefinition
    use Random
    
    type(Core), intent(in) :: c
    logical, intent(in) :: consistent
    type(HaplotypeLibrary), intent(in) :: library

    integer :: i, id

    do i = 1, c % getNAnisG()
      !Paternal Haps
      if (fullyPhased(c % getHaplotype(i, 1))) then
	call newHaplotype(c, i, 1, library)
      endif

      !Maternal Haps
      if (fullyPhased(c % getHaplotype(i, 2))) then
	!FUDGE FOR CONSISTENCY.  There should be no if statement here and then we could replace this with a call to newHaplotype
	if (.not.consistent .or. (library % getSize() > 0)) then
	  call c % setFullyPhased(i, 2)
	end if
	id = library % matchAddHap(c % getHaplotype(i, 2))
	call c % setHapAnis(i, 2, id)
      endif
    enddo

  end subroutine MakeHapLib

  subroutine ImputeFromLib(library, c, nGlobalHapsIter, PercGenoHaploDisagree, minHapFreq, consistent)
    ! Impute the phase for gametes that are not completely phased by LRP 
    ! by matching their phased loci to haplotypes in the Haplotype Library,
    ! following strategies listed in the section Step 2e of Hickey et al 2011.

    use HaplotypeLibraryDefinition
    use Constants
    use CoreDefinition
    use Clustering
    use Random
    
    type(HaplotypeLibrary), intent(in) :: library
    class(Core) :: c
    integer, intent(inout) :: nGlobalHapsIter
    double precision, intent(in) :: PercGenoHaploDisagree
    integer, intent(in) :: minHapFreq
    logical, intent(in) :: consistent

    integer :: i, j, k, nHapsOld, ErrorAllow, HapM, HapP
    integer, pointer, dimension(:,:) :: CandPairs

    integer, dimension(:), pointer :: compatHaps

    integer(kind = 1), pointer, dimension(:) :: comp
    integer, pointer, dimension(:) :: CandHapsPat, CandHapsMat, AllCandHapsMat, CandHaps, matches

    logical :: singlePat, singleMat

    ErrorAllow = int(PercGenoHaploDisagree * c % getNCoreSnp())

!    if ((nGlobalHapsIter == 1) .and. (ItterateType .eq. "Off")) then
    if (nGlobalHapsIter == 1) then
      print*, "   ", "Iteration ", nGlobalHapsIter, "found ", library % getSize(), "haplotypes"
    endif
    nHapsOld = 0

    do while (nHapsOld /= library % getSize())
      nGlobalHapsIter = nGlobalHapsIter + 1
      nHapsOld = library % getSize()

      do i = 1, c % getNAnisG()
	ErrorAllow = int(PercGenoHaploDisagree * c % numNotMissing(i))

	! If only one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): PATERNAL HAPLOTYPE
	if (c % getFullyPhased(i, 1) .and. (.not.c % getFullyPhased(i, 2))) then
	  call processComplement(c, i, library, 1, PercGenoHaploDisagree)
	end if

	! If only one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): MATERNAL HAPLOTYPE
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Haplotype 2 can get fully phased above and this will run despite both haplotypes now being phased
	! Affects results, likely due to error being allowed when phasing paternal from maternal in this step
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (c % getFullyPhased(i, 2) .and. (.not.c % getFullyPhased(i, 1))) then
	  call processComplement(c, i, library, 2, PercGenoHaploDisagree)
	end if

	! If neither of the gametes is completely phased (Section Step 2e.ii Hickey et al. 2011)
	if ((.not.c % getFullyPhased(i, 1)) .and. (.not.c % getFullyPhased(i, 2))) then
	  if (consistent) then
	    allocate(compatHaps(library % getSize()))
	    do k = 1, library % getSize()
	      compatHaps(k) = k
	    end do
	  else
	    compatHaps => library % getCompatHapsFreq(c % getSingleCoreGenos(i),minHapFreq, PercGenoHaploDisagree)
	  end if

	  CandHapsPat => library % limitedMatchWithError(c % getHaplotype(i, 1), ErrorAllow, compatHaps)
	  !Find all candiate maternal haps and then remove those that are already candidate paternal haps
	  AllCandHapsMat => library % limitedMatchWithError(c % getHaplotype(i, 2), ErrorAllow, compatHaps)
	  CandHapsMat => uniqueHaps(AllCandHapsMat, CandHapsPat)
	  deallocate(AllCandHapsMat)

	  deallocate(compatHaps)

	  HapP = 0
	  HapM = 0

	  ! If only have one paternal candidate haplotype, then the paternal haplotype is that one
	  if (size(CandHapsPat, 1) == 1) HapP = CandHapsPat(1)
	  ! If only have one maternal candidate haplotype, then the maternal haplotype is that one
	  if (size(CandHapsMat, 1) == 1) HapM = CandHapsMat(1)

	  ! If only one maternal candidate haplotype and many paternal candidate haplotypes
	  !! There's some odd logic here - if we have one paternal and maternal we end up keeping HapP as the
	  !! paternal even if it's not compitable with the mat.  May be checked later.
	  if ((size(CandHapsMat, 1) == 1).AND.(size(CandHapsPat, 1) > 0)) then
	    comp => complement(c % getSingleCoreGenos(i), library % getHap(HapM))
	    matches => library % limitedMatchWithError(comp, ErrorAllow, CandHapsPat)
	    if (size(matches) == 1) then
	      HapP = matches(1)
	    end if
	    deallocate(comp)
	    deallocate(matches)
	  end if

	  ! If only one paternal candidate haplotype and one / many maternal candidate haplotypes
	  if ((size(CandHapsPat, 1) == 1).and.(size(CandHapsMat, 1) > 0)) then
	    ! FUDGE TO DEAL WITH STRANGE LOGIC IN ORIGINAL CODE
	    HapM = 0

	    comp => complement(c % getSingleCoreGenos(i), library % getHap(HapP))
	    matches => library % limitedMatchWithError(comp, ErrorAllow, CandHapsMat)
	    if (size(matches) == 1) then
	      HapM = matches(1)
	    end if

	    ! FUDGE TO DEAL WITH STRANGE LOGIC IN ORIGINAL CODE
	    if (HapM /= CandHapsMat(size(CandHapsMat, 1))) then
	      HapM = 0
	    end if
	    
	    deallocate(comp)
	    deallocate(matches)
	  end if

	  ! If only have one paternal candidate haplotype
	  if (HapP /= 0) then
	    call matchedHaplotype(c, i, 1, library, HapP)

	    ! If no haplotype has been found for the maternal gamete, or 
	    ! there are more than one maternal candidate haplotype
	    if (HapM == 0) then
	      comp => complement(c % getSingleCoreGenos(i), c % getHaplotype(i, 1))

	      do j = 1, c % getNCoreSnp()
		if ((comp(j) == 0).or.(comp(j) == 1)) then
		  call c % setPhase(i, j, 2, comp(j))
		endif
	      enddo

	      if (fullyPhased(comp)) then
		call newHaplotype(c, i, 2, library)
	      end if
	      
	      deallocate(comp)
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
	    call matchedHaplotype(c, i, 2, library, HapM)

	    ! If no haplotype has been found for the paternal gamete, or 
	    ! there are more than one paternal candidate haplotype
	    if (HapP == 0) then
	      comp => complement(c % getSingleCoreGenos(i), c % getHaplotype(i, 2))

	      do j = 1, c % getNCoreSnp()
		if ((comp(j) == 0).or.(comp(j) == 1)) then
		  call c % setPhase(i, j, 1, comp(j))
		endif
	      enddo

	      if (fullyPhased(comp)) then
		call newHaplotype(c, i, 1, library)
	      end if
	      
	      deallocate(comp)
	    end if
	  end if

	  ! If the paternal and maternal gamete cannot be identify without ambiguity
	  ! (more than one or none at all)           
	  if ((HapP == 0).and.(HapM == 0)) then
	    ! Make array of all haps
	    allocate(CandHaps(size(CandHapsPat) + size(CandHapsMat)))
	    CandHaps(1:size(CandHapsPat)) = CandHapsPat
	    CandHaps(size(CandHapsPat) + 1:size(CandHaps)) = CandHapsMat

	    ! Get pairs of haplotypes that agree with the genotype
	    CandPairs => library % limitedCompatPairsWithError(c % getSingleCoreGenos(i), ErrorAllow, CandHaps, c % getNAnisG())

	    ! If we've found matching pairs                
	    if ((size(CandPairs, 1) > 0).and.((size(CandPairs, 1) * size(CandPairs, 1)) < c % getNAnisG())) then !Note the 200 number is a fudge

	      ! Check how many paternal candidates haplotypes
	      SinglePat = all(CandPairs(:, 1) == CandPairs(1, 1))

	      ! If there is only one paternal haplotype in all the candidate pairs
	      if (SinglePat) then
		! Phase the paternal gamete with this haplotype - if they're all the same it doesn't matter what pair we use
		call matchedHaplotype(c, i, 1, library, CandPairs(1, 1))

		! If only one haplotype is found for the paternal gamete 
		! and many for the maternal gamete, phase each loci only all pairs agree
		! (Step 2e.ii.B)
		do j = 1, c % getNCoreSnp()
		  if (agree(CandPairs(:, 2), library, j)) then
		    call c % setPhase(i, j, 2, library % getPhase(CandPairs(1, 2), j))
		  end if
		end do
		call c%setSwappable(i, 1)
	      end if

	      ! Check how many maternal candidates haplotypes
	      SingleMat = all(CandPairs(:, 2) == CandPairs(1, 2))

	      ! If there is only one maternal haplotype in all the candidate pairs
	      if (SingleMat) then
		call matchedHaplotype(c, i, 2, library, CandPairs(1, 2))

		! If only one haplotype is found for the paternal gamete 
		! and many for the maternal gamete, phase each loci only all pairs agree
		! (Step 2e.ii.C)
		do j = 1, c % getNCoreSnp()
		  if (agree(CandPairs(:, 1), library, j)) then
		    call c % setPhase(i, j, 1, library % getPhase(CandPairs(1, 1), j))
		  endif
		enddo
		call c%setSwappable(i, 2)
	      endif

	      ! If proband is not completely phased and have more than one candidate 
	      ! for both paternal and maternal haplotype 
	      ! (Step 2e.iv)
	      if (((.not.c % getFullyPhased(i, 1)) .or. (.not.c % getFullyPhased(i, 2))) &
		.and. (.not.SinglePat) .and. (.not.SingleMat)) then
		call clusterAndPhase(c, library, CandPairs, i)
		call c%setSwappable(i, 3)
	      endif
	    endif
	    deallocate(CandHaps)
	    deallocate(CandPairs)
	  endif
	  deallocate(CandHapsMat)
	  deallocate(CandHapsPat)
	endif
      end do
      print*, "   ", "Iteration ", nGlobalHapsIter, "found ", library % getSize(), "haplotypes"
    enddo

!    call abErrors(c)
    call complementPhaseSnps(c)

    call library % resetHapFreq()
    call c % resetFullyPhased()
    call c % resetHapAnis()

  end subroutine ImputeFromLib

  function agree(candidates, library, position) result (a)
    use HaplotypeLibraryDefinition

    integer, dimension(:), intent(in) :: candidates
    class(HaplotypeLibrary), intent(in) :: library
    integer, intent(in) :: position

    logical :: a
    integer(kind = 1), dimension(size(candidates, 1)) :: phases
    integer :: i

    do i = 1, size(candidates, 1)
      phases(i) = library % getPhase(candidates(i), position)
    end do

    a = (all(phases == 0) .or. (all(phases == 1)))
  end function agree

  subroutine processComplement(c, animal, library, fully, percgenohaplodisagree)
    use Constants
    use CoreDefinition
    use HaplotypeLibraryDefinition

    type(Core), intent(in) :: c
    integer, intent(in) :: animal, fully
    double precision, intent(in) :: percgenohaplodisagree
    type(HaplotypeLibrary) :: library

    integer(kind = 1), dimension(:), pointer :: comp
    integer, dimension(:), pointer :: CandHaps
    integer :: i, ErrorAllow, id, notfully

    ! Maps 1 -> 2 and 2 -> 1
    notfully = 3 - fully

    ErrorAllow = int(PercGenoHaploDisagree * c % numNotMissing(animal))

    comp => complement(c % getSingleCoreGenos(animal), c % getHaplotype(animal, fully))
    CandHaps => library % matchWithError(comp, ErrorAllow)

    if (size(CandHaps, 1) > 1) then
      do i = 1, c % getNCoreSnp()
	if (agree(CandHaps, library, i)) then
	  call c % setPhase(animal, i, notfully, library % getPhase(CandHaps(1), i))
	end if
      end do
    endif

    if (size(CandHaps, 1) == 1) then
      call c % setHaplotype(animal, notfully, library % getHap(CandHaps(1)))
      call c % setFullyPhased(animal, notfully)
      call library % incrementHapFreq(CandHaps(1))
      call c % setHapAnis(animal, notfully, CandHaps(1))
    end if

    if (size(CandHaps, 1) == 0) then
      do i = 1, c % getNCoreSnp()
	if ((comp(i) == 0).or.(comp(i) == 1)) then
	  call c % setPhase(animal, i, notfully, comp(i))
	endif
      enddo

      if (fullyPhased(comp)) then
	id = library % matchAddHap(c % getHaplotype(animal, notfully))
	call c % setHapAnis(animal, notfully, id)
	call c % setFullyPhased(animal, notfully)

      end if
    endif
    
    deallocate(comp)
    deallocate(CandHaps)
  end subroutine processComplement

  function complement(genos, haplotype) result (comp)
    use Constants
    integer(kind = 1), dimension(:), intent(in) :: genos, haplotype
    integer(kind = 1), dimension(:), pointer :: comp

    integer :: i

    allocate(comp(size(genos, 1)))

    do i = 1, size(genos, 1)
      if (genos(i) == MissingGenotypeCode) then
	comp(i) = 9
      else
	comp(i) = genos(i) - haplotype(i)
      end if
    end do
  end function complement

  function HapsToCluster(CandPairs) result (ToCluster)
    integer, dimension(:,:), intent(in) :: CandPairs
    integer, dimension(:), pointer :: ToCluster

    integer, dimension(:), allocatable :: tempToCluster
    integer :: numToCluster, i, j

    integer(8) :: length, s

    allocate(tempToCluster(2 * size(CandPairs)))
    tempToCluster = 0
    numToCluster = 0

    do i = 1, size(CandPairs, 1)
      do j = 1, 2
	if (all(tempToCluster /= CandPairs(i, j))) then
	  numToCluster = numToCluster + 1
	  tempToCluster(numToCluster) = CandPairs(i, j)
	end if
      end do
    end do

    allocate(toCluster(numToCluster))
    toCluster = tempToCluster(1:numToCluster)



    length = size(toCluster, 1)
    s = 4
    call qsort(toCluster, length, s, cmp_function)
  end function HapsToCluster

  function cmp_function(a1, a2) result (cmp)
    integer(2), intent(in) :: a1, a2
    integer(2) :: cmp

    cmp = a1 - a2
  end function

  function findloc(array, val) result(loc)
    integer, dimension(:), intent(in) :: array
    integer, intent(in) :: val
    integer :: loc

    integer :: i

    loc = 0
    do i = 1, size(array)
      if (array(i) == val) then
	loc = i
	exit
      end if
    end do
  end function findloc

  subroutine matchedHaplotype(c, animal, phase, library, id)
    use CoreDefinition
    use HaplotypeLibraryDefinition

    class(Core) :: c
    integer, intent(in) :: animal, phase, id
    class(HaplotypeLibrary) :: library

    call c % setHaplotype(animal, phase, library % getHap(id))
    call c % setFullyPhased(animal, phase)
    call library % incrementHapFreq(id)
    call c % setHapAnis(animal, phase, id)
  end subroutine matchedHaplotype

  subroutine newHaplotype(c, animal, phase, library)
    use CoreDefinition
    use HaplotypeLibraryDefinition

    class(Core) :: c
    integer, intent(in) :: animal, phase
    class(HaplotypeLibrary) :: library

    integer :: id

    id = library % matchAddHap(c % getHaplotype(animal, phase))
    call c % setHapAnis(animal, phase, id)
    call c % setFullyPhased(animal, phase)
  end subroutine newHaplotype

  subroutine clusterAndPhase(c, library, CandPairs, animal)
    use CoreDefinition
    use HaplotypeLibraryDefinition
    use Clustering
    use Constants

    class(Core) :: c
    class(HaplotypeLibrary) :: library
    integer, dimension(:,:), intent(in) :: CandPairs
    integer, intent(in) :: animal

    integer :: nSnp

    integer(kind = 1), dimension (:,:), pointer :: TempHapArray
    integer, dimension(:), pointer :: ClusterMember
    integer, dimension(:), pointer :: TempHapVector
    integer :: rounds

    integer :: countZero, countOne

    integer :: j, k

    nSNp = c % getNCoreSnp()

    ! Initialize procedure of k-medoids
    TempHapVector => HapsToCluster(CandPairs)
    TempHapArray => library % getHaps(TempHapVector)
    ClusterMember => initialAssignmentAlternate(size(TempHapVector, 1))
    rounds = cluster(TempHapArray, ClusterMember, 2, nMaxRounds, .false.)

    if (rounds <= nMaxRounds) then
      if (count(ClusterMember(:) == 2) == 1) then
	call matchedHaplotype(c, animal, 2, library, TempHapVector(findloc(ClusterMember, 2)))
      end if
      if (count(ClusterMember(:) == 1) == 1) then
	call matchedHaplotype(c, animal, 1, library, TempHapVector(findloc(ClusterMember, 1)))
      end if
      if ((count(ClusterMember(:) == 2) > 1).and.(count(ClusterMember(:) == 2) > 1)) then
	call c % setHaplotypeToUnphased(animal, 1)
	call c % setHaplotypeToUnphased(animal, 2)
	do j = 1, nSnp
	  if (c % getCoreGeno(animal, j) == 0) then
	    call c % setPhase(animal, j, 1, 0)
	    call c % setPhase(animal, j, 2, 0)
	  end if
	  if (c % getCoreGeno(animal, j) == 2) then
	    call c % setPhase(animal, j, 1, 1)
	    call c % setPhase(animal, j, 2, 1)
	  end if
	  CountZero = 0
	  CountOne = 0
	  do k = 1, size(ClusterMember, 1)
	    if (ClusterMember(k) == 2) then
	      if (library % getPhase(TempHapVector(k), j) == 0)&
	      CountZero = CountZero + 1
	      if (library % getPhase(TempHapVector(k), j) == 1)&
	      CountOne = CountOne + 1
	    endif
	  end do
	  if ((CountZero == 0).and.(CountOne > 0)) call c % setPhase(animal, j, 2, 1)
	  if ((CountZero > 0).and.(CountOne == 0)) call c % setPhase(animal, j, 2, 0)

	  CountZero = 0
	  CountOne = 0
	  do k = 1, size(ClusterMember, 1)
	    if (ClusterMember(k) == 1) then
	      if (library % getPhase(TempHapVector(k), j) == 0) CountZero = CountZero + 1
	      if (library % getPhase(TempHapVector(k), j) == 1) CountOne = CountOne + 1
	    endif
	  end do
	  if ((CountZero == 0).and.(CountOne > 0)) call c % setPhase(animal, j, 1, 1)
	  if ((CountZero > 0).and.(CountOne == 0)) call c % setPhase(animal, j, 1, 0)
	end do
      endif
    end if
    deallocate(ClusterMember)
    deallocate(TempHapArray)
    deallocate(TempHapVector)
  end subroutine clusterAndPhase

  subroutine abErrors(c)
    use CoreDefinition
    use Constants

    class(Core) :: c

    integer :: countA, countB, ErrorCountAB
    integer :: i, j
    integer(kind = 1) :: val

    ErrorCountAB = int(c % getNCoreSnp() * 0.09)

    do i = 1, c % getNAnisG()
      CountA = 0
      CountB = 0
      do j = 1, c % getNCoreSnp()
	if (c % getCoreGeno(i, j) /= MissingGenotypeCode) then
	  if ((c % getPhase(i, j, 1) /= MissingPhaseCode).and.(c % getPhase(i, j, 2) == MissingPhaseCode)) then
	    val = c % getCoreGeno(i, j) - c % getPhase(i, j, 1)
	    if ((val == 0).or.(val == 1)) then !here 7th april 2011
	      call c % setPhase(i, j, 2, val)
	    else
	      CountA = CountA + 1
	    endif
	  endif
	  if ((c % getPhase(i, j, 2) /= MissingPhaseCode).and.(c % getPhase(i, j, 1) == MissingPhaseCode)) then
	    val = c % getCoreGeno(i, j) - c % getPhase(i, j, 2)
	    if ((val == 0).or.(val == 1)) then !here 7th april 2011
	      call c % setPhase(i, j, 1, val)
	    else
	      CountB = CountB + 1
	    endif
	  endif
	end if
      end do
      if ((CountA > ErrorCountAB) .or. (CountB > ErrorCountAB)) then
	call c % setHaplotypeToUnphased(i, 1)
	call c % setHaplotypeToUnphased(i, 2)
	do j = 1, c % getNCoreSnp()
	  if (c % getCoreGeno(i, j) == 0) then
	    call c % setPhase(i, j, 1, 0)
	    call c % setPhase(i, j, 2, 0)
	  end if
	  if (c % getCoreGeno(i, j) == 2) then
	    call c % setPhase(i, j, 1, 1)
	    call c % setPhase(i, j, 1, 1)
	  end if
	enddo
      endif
    end do
  end subroutine abErrors

  subroutine complementPhaseSnps(c)
    use CoreDefinition

    class(Core) :: c

    integer :: i, j


    do i = 1, c % getNAnisG()
      do j = 1, c % getNCoreSnp()
	if (c % getCoreGeno(i, j) == 1) then
	  if ((c % getPhase(i, j, 1) == MissingPhaseCode).and.(c % getPhase(i, j, 2) /= MissingPhaseCode)) then
	    call c % setPhase(i, j, 1, c % getCoreGeno(i, j) - c % getPhase(i, j, 2))
	  end if
	  if ((c % getPhase(i, j, 2) == MissingPhaseCode).and.(c % getPhase(i, j, 1) /= MissingPhaseCode)) then
	    call c % setPhase(i, j, 2, c % getCoreGeno(i, j) - c % getPhase(i, j, 1))
	  end if
	endif
	if (c % getCoreGeno(i, j) == 0) then
	  if ((c % getPhase(i, j, 1) == MissingPhaseCode).and.(c % getPhase(i, j, 2) /= MissingPhaseCode)) then
	    call c % setPhase(i, j, 1, 0)
	  end if
	  if ((c % getPhase(i, j, 2) == MissingPhaseCode).and.(c % getPhase(i, j, 1) /= MissingPhaseCode)) then
	    call c % setPhase(i, j, 2, 0)
	  end if
	endif
	if (c % getCoreGeno(i, j) == 2) then
	  if ((c % getPhase(i, j, 1) == MissingPhaseCode).and.(c % getPhase(i, j, 2) /= MissingPhaseCode)) then
	    call c % setPhase(i, j, 1, 1)
	  end if
	  if ((c % getPhase(i, j, 2) == MissingPhaseCode).and.(c % getPhase(i, j, 1) /= MissingPhaseCode)) then
	    call c % setPhase(i, j, 2, 1)
	  end if
	endif
      enddo
    enddo
  end subroutine complementPhaseSnps
  
  function fullyPhased(haplotype) result (fully)
    integer(kind=1), dimension(:), intent(in) :: haplotype
    logical :: fully

    fully = all((haplotype == 0) .or. (haplotype == 1))
  end function fullyPhased
  
  function uniquehaps(haps1, haps2) result (uniq)
    integer, dimension(:), intent(in) :: haps1, haps2
    integer, dimension(:), pointer :: uniq

    integer, dimension(:), allocatable :: tempU

    integer :: i, p

    allocate(tempU(size(haps1)))
    p = 0
    do i = 1, size(haps1)
      if (.not. any (haps1(i) == haps2)) then
	p = p + 1
	tempU(p) = haps1(i)
      end if
    end do

    allocate (uniq(p))
    uniq = tempU(1:p)
    deallocate(tempU)
  end function uniquehaps

end module HaplotypeLibraryPhasing