module HaplotypeLibraryPhasing
  use ConstantModule
  implicit none

  integer, parameter, private :: nMaxRounds = 100
  
contains
  subroutine UpdateHapLib(c, library)
    use HaplotypeLibraryDefinition
    use HaplotypeModule
    use CoreDefinition
    use Random
    
    type(Core), intent(in) :: c
    type(HaplotypeLibrary), intent(in) :: library
    
    type(Haplotype) :: hap

    integer :: i
    
    do i = 1, c % getNAnisG()
      !Paternal Haps
      hap = c % phase(i, 1)
      if (hap % fullyPhased() ) then
	call newHaplotype(c, i, 1, library)
      endif

      !Maternal Haps
      hap = c % phase(i, 2)
      if (hap%fullyPhased()) then
	call newHaplotype(c, i, 2, library)
      endif
    enddo
    
  end subroutine UpdateHapLib

  subroutine ImputeFromLib(library, c, nGlobalHapsIter, PercGenoHaploDisagree, minHapFreq, quiet)
    ! Impute the phase for gametes that are not completely phased by LRP 
    ! by matching their phased loci to haplotypes in the Haplotype Library,
    ! following strategies listed in the section Step 2e of Hickey et al 2011.

    use HaplotypeLibraryDefinition
    use CoreDefinition
    use ClusteringModule
    use Random
    use GenotypeModule
    use HaplotypeModule
    
    type(HaplotypeLibrary), intent(in) :: library
    class(Core) :: c
    integer, intent(inout) :: nGlobalHapsIter
    double precision, intent(in) :: PercGenoHaploDisagree
    integer, intent(in) :: minHapFreq
    logical, intent(in) :: quiet

    integer :: i, nHapsOld, ErrorAllow, HapM, HapP
    integer, pointer, dimension(:,:) :: CandPairs

    integer, dimension(:), pointer :: compatHaps

    type(Haplotype) :: comp
    integer, pointer, dimension(:) :: CandHapsPat, CandHapsMat, AllCandHapsMat, CandHaps, matches
    
    type(Genotype), pointer :: geno
    type(Haplotype), pointer :: hap1, hap2

    logical :: singlePat, singleMat, oneNotTwo, twoNotOne

    ErrorAllow = int(PercGenoHaploDisagree * c % getNCoreSnp())

    if ((nGlobalHapsIter == 1) .and. (.not. quiet)) then
      print*, "   ", "Iteration ", nGlobalHapsIter, "found ", library % getSize(), "haplotypes"
    end if
    nHapsOld = 0

    do while (nHapsOld /= library % getSize())
      nGlobalHapsIter = nGlobalHapsIter + 1
      nHapsOld = library % getSize()

      do i = 1, c % getNAnisG()
	geno => c%coreGenos(i)
	hap1 => c%phase(i,1)
	hap2 => c%phase(i,2)
	ErrorAllow = int(PercGenoHaploDisagree * geno % numNotMissing())
	
	oneNotTwo = c % getFullyPhased(i, 1) .and. (.not.c % getFullyPhased(i, 2))
	twoNotOne = c % getFullyPhased(i, 2) .and. (.not.c % getFullyPhased(i, 1))

	! If only one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): PATERNAL HAPLOTYPE
	if (oneNotTwo) then
	  call processComplement(c, i, library, 1, PercGenoHaploDisagree)
	end if

	! If only one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): MATERNAL HAPLOTYPE
	if (twoNotOne) then
	  call processComplement(c, i, library, 2, PercGenoHaploDisagree)
	end if

	! If neither of the gametes is completely phased (Section Step 2e.ii Hickey et al. 2011)
	if ((.not.c % getFullyPhased(i, 1)) .and. (.not.c % getFullyPhased(i, 2))) then
          compatHaps => library % getCompatHapsFreq(geno,minHapFreq, PercGenoHaploDisagree)

	  CandHapsPat => library % limitedMatchWithError(c % phase(i, 1), ErrorAllow, compatHaps)
	  !Find all candiate maternal haps and then remove those that are already candidate paternal haps
	  AllCandHapsMat => library % limitedMatchWithError(c % phase(i, 2), ErrorAllow, compatHaps)
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
	  if ((size(CandHapsMat, 1) == 1).AND.(size(CandHapsPat, 1) > 1)) then
	    comp = geno%complement(library%getHap(HapM))
	    matches => library % limitedMatchWithError(comp, ErrorAllow, CandHapsPat)
	    if (size(matches) == 1) then
	      HapP = matches(1)
	    end if
	    deallocate(matches)
	  end if

	  ! If only one paternal candidate haplotype and many maternal candidate haplotypes
	  if ((size(CandHapsPat, 1) == 1).and.(size(CandHapsMat, 1) > 1)) then
	    comp = geno%complement(library%getHap(HapP))
	    matches => library % limitedMatchWithError(comp, ErrorAllow, CandHapsMat)
	    if (size(matches) == 1) then
	      HapM = matches(1)
	    end if
	    deallocate(matches)
	  end if
	  
	  ! If one paternal candidate and one maternal candidate and they are incompatible use neither
	  if ((size(CandHapsPat, 1) == 1).and.(size(CandHapsMat, 1) == 1)) then
	    if (.not. geno%compatibleHaplotypes(library%getHap(HapP), library%getHap(HapM), ErrorAllow)) then
	      HapP = 0
	      HapM = 0
	    end if
	  end if

	  ! If only have one paternal candidate haplotype
	  if (HapP /= 0) then
	    call matchedHaplotype(c, i, 1, library, HapP)

	    ! If no haplotype has been found for the maternal gamete, or 
	    ! there are more than one maternal candidate haplotype
	    if (HapM == 0) then
	      comp = geno%complement(c%phase(i,1))
	      call hap2%setFromOther(comp)

	      if (comp%fullyPhased()) then
		call newHaplotype(c, i, 2, library)
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
	    call matchedHaplotype(c, i, 2, library, HapM)

	    ! If no haplotype has been found for the paternal gamete, or 
	    ! there are more than one paternal candidate haplotype
	    if (HapP == 0) then
	      comp = geno%complement(c%phase(i,2))
	      call hap1%setFromOther(comp)

	      if (comp%fullyPhased()) then
		call newHaplotype(c, i, 1, library)
	      end if
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
	    CandPairs => library % limitedCompatPairsWithError(c % coreGenos(i), ErrorAllow, CandHaps, c % getNAnisG())

	    ! If we've found matching pairs                
	    if ((size(CandPairs, 1) > 0).and.((size(CandPairs, 1) * size(CandPairs, 1)) < c % getNAnisG())) then

	      ! Check how many paternal candidates haplotypes
	      SinglePat = all(CandPairs(:, 1) == CandPairs(1, 1))

	      ! If there is only one paternal haplotype in all the candidate pairs
	      if (SinglePat) then
		! Phase the paternal gamete with this haplotype - if they're all the same it doesn't matter what pair we use
		call matchedHaplotype(c, i, 1, library, CandPairs(1, 1))

		! If only one haplotype is found for the paternal gamete 
		! and many for the maternal gamete, phase each loci only all pairs agree
		! (Step 2e.ii.B)
		call hap2%setOneBits(library%allOne(CandPairs(:,2)))
		call hap2%setZeroBits(library%allZero(CandPairs(:,2)))
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
		call hap1%setOneBits(library%allOne(CandPairs(:,1)))
		call hap1%setZeroBits(library%allZero(CandPairs(:,1)))
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
      if (.not. quiet) then
	print*, "   ", "Iteration ", nGlobalHapsIter, "found ", library % getSize(), "haplotypes"
      end if
    enddo

    call complementPhaseSnps(c)

    call library % resetHapFreq()
    call c % resetFullyPhased()
    call c % resetHapAnis()

  end subroutine ImputeFromLib

  subroutine processComplement(c, animal, library, fully, percgenohaplodisagree)
    use CoreDefinition
    use HaplotypeLibraryDefinition
    use GenotypeModule
    use HaplotypeModule

    type(Core), intent(in) :: c
    integer, intent(in) :: animal, fully
    double precision, intent(in) :: percgenohaplodisagree
    type(HaplotypeLibrary) :: library

    type(Haplotype) :: comp
    type(Haplotype), pointer :: hap
    type(Genotype), pointer :: geno
    integer, dimension(:), pointer :: CandHaps
    integer :: ErrorAllow, notfully

    ! Maps 1 -> 2 and 2 -> 1
    notfully = 3 - fully

    geno => c %coreGenos(animal)
    ErrorAllow = int(PercGenoHaploDisagree * geno % numNotMissing())
    comp = geno%complement(c % phase(animal, fully))
    CandHaps => library % matchWithError(comp, ErrorAllow)
    
    hap => c%phase(animal, notfully)

    if (size(CandHaps, 1) > 1) then
      call hap%setOneBits(library%allOne(CandHaps))
      call hap%setZeroBits(library%allZero(CandHaps))
    endif

    if (size(CandHaps, 1) == 1) then
      call c % setHaplotype(animal, notfully, library % getHap(CandHaps(1)))
      call c % setFullyPhased(animal, notfully)
      call library % incrementHapFreq(CandHaps(1))
      call c % setHapAnis(animal, notfully, CandHaps(1))
    end if

    if (size(CandHaps, 1) == 0) then
      call hap%setFromOther(comp)

      if (comp%fullyPhased()) then
	call newHaplotype(c, animal, notfully, library)

      end if
    endif
    
    deallocate(CandHaps)
  end subroutine processComplement

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
    use HaplotypeModule
    use CoreDefinition
    use HaplotypeLibraryDefinition

    class(Core) :: c
    integer, intent(in) :: animal, phase
    class(HaplotypeLibrary) :: library

    integer :: id

    id = library % matchAddHap(c % phase(animal, phase))
    call c % setHapAnis(animal, phase, id)
    call c % setFullyPhased(animal, phase)
  end subroutine newHaplotype

  subroutine clusterAndPhase(c, library, CandPairs, animal)
    use CoreDefinition
    use HaplotypeLibraryDefinition
    use HaplotypeModule
    use GenotypeModule
    use ClusteringModule

    class(Core) :: c
    class(HaplotypeLibrary) :: library
    integer, dimension(:,:), intent(in) :: CandPairs
    integer, intent(in) :: animal

    integer :: nSnp

    integer(kind = 1), dimension (:,:), allocatable :: TempHapArray
    integer, dimension(:), pointer :: ClusterMember
    integer, dimension(:), pointer :: TempHapVector
    integer :: rounds


    integer :: i, i1, i2
    type(Haplotype) :: h
    type(Haplotype), pointer :: hap1, hap2
    type(Genotype), pointer :: geno
    integer, dimension(:), allocatable :: ids1, ids2

    nSNp = c % getNCoreSnp()
    
    geno => c%coreGenos(animal)
    hap1 => c%phase(animal,1)
    hap2 => c%phase(animal,2)

    ! Initialize procedure of k-medoids
    TempHapVector => HapsToCluster(CandPairs)
    allocate(TempHapArray(size(TempHapVector),nSnp))
    do i = 1, size(TempHapVector)
      h = library%getHap(TempHapVector(i))
      TempHapArray(i,:) = h%toIntegerArray()
    end do
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
	call geno%setHaplotype(hap1)
	call geno%setHaplotype(hap2)
	
	i1 = 0
	i2 = 0
	allocate(ids1(count(ClusterMember(:) == 1)))
	allocate(ids2(count(ClusterMember(:) == 2)))
	do i = 1, size(ClusterMember,1)
	  if (ClusterMember(i) == 1) then
	    i1 = i1 + 1
	    ids1(i1) = TempHapVector(i)
	  else
	    i2 = i2 + 1
	    ids2(i2) = TempHapVector(i)
	  end if
	end do
		
	call hap1%setZeroBits(library%oneZeroNoOnes(ids1))
	call hap1%setOneBits(library%oneOneNoZeros(ids1))
	
	call hap2%setZeroBits(library%oneZeroNoOnes(ids2))
	call hap2%setOneBits(library%oneOneNoZeros(ids2))
      endif
    end if
    deallocate(ClusterMember)
    deallocate(TempHapArray)
    deallocate(TempHapVector)
  end subroutine clusterAndPhase

  subroutine complementPhaseSnps(c)
    use CoreDefinition
    use GenotypeModule
    use HaplotypeModule

    class(Core) :: c

    integer :: i
    
    type(Genotype), pointer :: geno
    type(Haplotype), pointer :: hap1, hap2
    type(Haplotype) :: hap1complement, hap2complement
    
    do i = 1, c % getNAnisG()
      geno => c%getCoreGenos(i)
      hap1 => c%phase(i,1)
      hap2 => c%phase(i,2)
      
      hap1complement = geno%complement(hap1)
      hap2complement = geno%complement(hap2)      
      
      call Hap1%setFromOtherIfMissing(hap2complement)
      call Hap2%setFromOtherIfMissing(hap1complement)
    enddo
  end subroutine complementPhaseSnps
  
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