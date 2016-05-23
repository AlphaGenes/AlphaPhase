module HaplotypeLibrary
  implicit none
  public

  type, public :: HapLib
    private
    integer(kind = 1), dimension (:,:), allocatable :: store
    integer, dimension(:), allocatable :: hapFreq
    integer :: size
    integer :: nSnps
    integer :: storeSize, stepSize
    integer, dimension(:), allocatable :: randomOrder
  contains
    private
    procedure, public :: initalise
    procedure, public :: hasHap
    procedure :: addHap
    procedure, public :: getHap
    procedure, public :: matchWithError
    procedure, public :: limitedMatchWithError
    procedure, public :: limitedCompatPairsWithError
    procedure, public :: getPhase
    procedure, public :: getSize
    procedure, public :: getHapRel
    procedure, public :: getNumSnps
    procedure, public :: resetHapFreq
    procedure, public :: incrementHapFreq
    procedure, public :: getHapFreq
    procedure, public :: getCompatHaps    
    procedure, public :: matchAddHap
  end type HapLib

contains
  subroutine initalise(library, nSnps, storeSize, stepSize)
    use Random
    class(HapLib) :: library
    integer, intent(in) :: nSnps
    integer, intent(in) :: storeSize
    integer, intent(in) :: stepSize

    integer :: nCount, secs

    library % nSnps = nSnps
    library % size = 0
    library % storeSize = storeSize
    library % stepSize = stepSize
    if (allocated(library%store)) then
      deallocate(library%store)
      deallocate(library%hapFreq)
    end if
    allocate(library % store(storeSize, nSnps))
    allocate(library % hapFreq(storeSize))
    library % store = 0
    library % hapFreq = 0
    call system_clock(nCount)
    secs = mod(nCount, int(1e6))
    if (allocated(library%randomOrder)) then
      deallocate(library%randomOrder)
    end if
    allocate(library % randomOrder(library % nSnps))
    call RandomOrder(library % randomOrder, library % nSnps, 1, -abs(secs))
  end subroutine initalise

  function hasHap(library, haplotype) result(id)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer :: id

    integer :: i, j
    logical :: match

    id = 0

    do i = 1, library % size
      match = .true.
      do j = 1, library % nSnps
	if (library % store(i, library % RandomOrder(j)) /= haplotype(library % RandomOrder(j))) then
	  match = .false.
	  exit
	end if
      end do
      if (match) then
	id = i
	exit
      end if
    end do
  end function hasHap

  function addHap(library, haplotype) result(id)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer :: id

    integer :: newStoreSize
    integer(kind = 1), dimension(:,:), allocatable :: tempStore
    integer, dimension(:), allocatable :: tempHapFreq

    if (library % Size == library % storeSize) then
      newStoreSize = library % storeSize + library % stepSize
      
      allocate(tempStore(library % storeSize, library % nSnps))
      tempStore = library%store
      deallocate(library%store)
      allocate(library%store(newStoreSize, library % nSnps))
      library % store = 0
      library % store(1:library % Size,:) = tempStore
      deallocate(tempStore)
      
      allocate(tempHapFreq(library % storeSize))
      tempHapFreq = library%hapFreq
      deallocate(library%hapFreq)
      allocate(library%hapFreq(newStoreSize))
      library % hapFreq = 0
      library % hapFreq(1:library % Size) = tempHapFreq
      deallocate(tempHapFreq)
      
      library % StoreSize = newStoreSize
    end if

    library % Size = library % Size + 1
    library % Store(library % Size,:) = haplotype
    
    library%hapfreq(library%size) = 1
    id = library%size
  end function addHap
  
  function matchAddHap(library, haplotype) result (id)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer :: id
    
    id = library%hasHap(haplotype)
    if (id == 0) then
      id = library%addHap(haplotype)
    else
      library%hapfreq(id) = library%hapfreq(id) + 1
    end if    
  end function matchAddHap

  function matchWithError(library, haplotype, allowedError) result(matches)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer, intent(in) :: allowedError
    integer, dimension(:), pointer :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, j, e, num
    logical :: match

    allocate(tempMatches(library % size))

    num = 0

    do i = 1, library % size
      e = 0
      match = .true.
      do j = 1, library % nSnps
	if (haplotype(library % RandomOrder(j)) /= 9) then
	  if (library % store(i, library % RandomOrder(j)) /= haplotype(library % RandomOrder(j))) then
	    e = e + 1
	    if (e > allowedError) then
	      match = .false.
	      exit
	    endif
	  end if
	end if
      end do
      if (match) then
	num = num + 1
	tempMatches(num) = i
      end if
    end do
    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function matchWithError
  
  function limitedMatchWithError(library, haplotype, allowedError, limit) result(matches)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer, intent(in) :: allowedError
    integer, dimension(:), intent(in) :: limit
    integer, dimension(:), pointer :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, j, k, e, num
    logical :: match

    allocate(tempMatches(library % size))

    num = 0

    do k = 1, size(limit)
      i = limit(k)
      e = 0
      match = .true.
      do j = 1, library % nSnps
	if (haplotype(library % RandomOrder(j)) /= 9) then
	  if (library % store(i, library % RandomOrder(j)) /= haplotype(library % RandomOrder(j))) then
	    e = e + 1
	    if (e > allowedError) then
	      match = .false.
	      exit
	    endif
	  end if
	end if
      end do
      if (match) then
	num = num + 1
	tempMatches(num) = i
      end if
    end do
    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function limitedMatchWithError
  
  function limitedCompatPairsWithError(library, genos, ErrorAllow, limit, nAnisG) result(pairs)
    use Constants
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: genos
    integer, intent(in) :: ErrorAllow
    integer, dimension(:), intent(in) :: limit
    integer, intent(in) :: nAnisG
    integer, dimension(:,:), pointer :: pairs
    
    integer, dimension(:,:), pointer :: tempPairs
    integer :: i, j, k, p, e
    logical :: match
    
    allocate(tempPairs(nAnisG*2,2))
    
    p = 0
    i = 1
    do while ((i <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
      j = i + 1
      do while ((j <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
	match = .true.
	do k = 1, library%nSnps
	  if (Genos(library % randomOrder(k)) /= MissingGenotypeCode)then
	    if ((library%getPhase(limit(i), library % randomOrder(k)) &
	    + library%getPhase(limit(j), library % randomOrder(k))) &
	    /= Genos(library%randomOrder(k))) then
	      e = e + 1
	      if (e > ErrorAllow) then
		match = .false.
		exit
	      end if
	    end if
	  endif
	end do
	
	if (match) then
	  p = p + 1
	  tempPairs(p,1) = limit(i)
	  tempPairs(p,2) = limit(j)
	end if
	j = j + 1
      end do
      i = i + 1
    end do
    
    allocate(pairs(p,2))
    pairs = tempPairs(1:p,:)
    deallocate(tempPairs)
  end function limitedCompatPairsWithError
	

  function getHap(library, id) result(hap)
    class(HapLib) :: library
    integer, intent(in) :: id
    integer(kind = 1), dimension(:), allocatable :: hap

    allocate(hap(library % nSnps))
    hap = library % store(id,:)
  end function getHap

  function getPhase(library, id, snp) result(phase)
    class(HapLib) :: library
    integer, intent(in) :: id, snp
    integer(kind = 1) :: phase
  
    phase = library % store(id,snp)
  end function getPhase

  function getSize(library) result(size)
    class(HapLib) :: library
    integer :: size

    size = library % size
  end function getSize

  function getHapRel(library) result (rel)
    class(HapLib) :: library
    integer, allocatable, dimension (:,:) :: rel

    integer :: i, j, k

    integer :: counter

    allocate(rel(library % size, library % size))

    rel = 0
    do i = 1, library % size
      do j = i + 1, library % size
	counter = 0

	do k = 1, library % nSnps
	  if (library % store(i, k) == library % store(j, k)) counter = counter + 1
	end do
	rel(i, j) = counter
	rel(j, i) = counter
      end do
    end do
  end function getHapRel


  function getNumSnps(library) result (num)
    class(HapLib) :: library
    integer :: num

    num = library % nSnps
  end function getNumSnps
  
  subroutine resetHapFreq(library)
    class(HapLib) :: library
    
    library%HapFreq = 0
  end subroutine resetHapFreq
  
  subroutine incrementHapFreq(library, id)
    class(HapLib) :: library
    integer, intent(in) :: id
    
    library%HapFreq(id) = library%HapFreq(id) + 1
  end subroutine incrementHapFreq
  
  function getHapFreq(library, id) result (freq)
    class(HapLib) :: library
    integer, intent(in) :: id
    integer :: freq
    
    freq = library%hapFreq(id)
  end function getHapFreq
  
  function getCompatHaps(library, genos) result (compatHaps)
    use parameters, only : percgenohaplodisagree
    
    class(HapLib) :: library
    integer(kind=1), dimension(:), intent(in) :: genos
    integer, dimension(:), pointer :: compatHaps
    
    integer, dimension(:), allocatable :: tempCompatHaps
    integer :: k, j, numCompatHaps, disagree, ErrorAllow
    
    ErrorAllow = int(PercGenoHaploDisagree * library%nSnps)
    allocate(tempCompatHaps(library%size))
    numCompatHaps = 0
    do k = 1, library%size
      Disagree = 0
      do j = 1, library%nSnps
	if ((library%getPhase(k, library%randomOrder(j)) == 0) .and. (genos(library%randomOrder(j)) == 2)) then
	  Disagree = Disagree + 1
	  if (Disagree > ErrorAllow) then
	    exit
	  end if
	end if
	if ((library%getPhase(k, library%randomOrder(j)) == 1) .and. (genos(library%randomOrder(j)) == 0)) then
	  Disagree = Disagree + 1
	  if (Disagree > ErrorAllow) then
	    exit
	  end if
	end if
      end do
      if (Disagree <= ErrorAllow) then
	numCompatHaps = numCompatHaps + 1
	tempCompatHaps(numCompatHaps) = k
      end if
    end do
    allocate(compatHaps(numCompatHaps))
    compatHaps = tempCompatHaps(1:numCompatHaps)
    deallocate(tempCompatHaps)
  end function getCompatHaps
  
  
  
  
  
  
  !MESSY! SHOULDN'T BE HERE!!!!
  !subroutine CheckCompatHapGeno(genos, phase)
  subroutine CheckCompatHapGeno(c)
  use Parameters, only: PercGenoHaploDisagree
  use Constants
  use CoreSubsetDefinition
  implicit none

  class(CoreSubset) :: c

  integer(kind=1), dimension(:,:), pointer :: genos
  
  integer :: i, j, CountError, ErrorAllow, counter, counterMissing, nAnisG, nCoreSnp
  double precision :: value, Yield

  nAnisG  = c%getNAnisG()
  nCoreSnp = c%getNCoreSnp()
  
  genos => c%getCoreGenos()
  
  ErrorAllow = int(PercGenoHaploDisagree * nCoreSnp)

  do i = 1, nAnisG
    CountError = 0
    counterMissing = 0
    do j = 1, nCoreSnp
      if ((c%getPhase(i, j, 1) /= 9).and.(c%getPhase(i, j, 2) /= 9)) then
	counterMissing = counterMissing + 1
	if ((Genos(i, j) /= MissingGenotypeCode).and.(c%getPhaseGeno(i,j)  /= Genos(i, j))) CountError = CountError + 1
      end if
    end do
    ErrorAllow = int(PercGenoHaploDisagree * counterMissing)
    if (CountError >= ErrorAllow) then
      do j = 1, nCoreSnp
	if (Genos(i, j) /= MissingGenotypeCode) then
	  if ((c%getPhase(i, j, 1) /= 9).and.(c%getPhase(i, j, 2) /= 9).and.(c%getPhaseGeno(i, j) /= Genos(i, j))) then
	    if (Genos(i, j) == 1) call c%setPhase(i, j, 2, 9)
	    if (Genos(i, j) == MissingGenotypeCode) call c%setPhase(i, j, 2, 9)
	    if (Genos(i, j) == 0) then
	      call c%setPhase(i, j, 1, 0)
	      call c%setPhase(i, j, 2, 0)
	    end if
	    if (Genos(i, j) == 2) then
	      call c%setPhase(i, j, 1, 1)
	      call c%setPhase(i, j, 2, 1)
	    end if
	  endif
	endif
      enddo
    endif
  end do

  print*, " "
  write (*, '(a3,f6.2,a45)') "  ", c%getYield(1), "% was the Paternal allele yield for this core"
  write (*, '(a3,f6.2,a45)') "  ", c%getYield(2), "% was the Maternal allele yield for this core"

end subroutine CheckCompatHapGeno

subroutine MakeHapLib(library, c)
  use CoreDefinition
  use Random
  use Parameters
  implicit none

  type(HapLib), intent(in) :: library
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
  use Constants
  use CoreDefinition
  use Clustering
  use Random
  implicit none
  
  type(HapLib), intent(in) :: library
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

subroutine WriteHapLib(library, currentcore, c)
  use Parameters, only: fullfileoutput
  use Constants
  use CoreDefinition
  implicit none
  
  type(HapLib), intent(in) :: library
  type(Core), intent(in) :: c
  integer, intent(in) :: currentcore

  integer :: i, j, k, counter, SizeCore, nHaps !, nAnisG
  character(len = 300) :: filout
  
  SizeCore = library%getNumSnps()
  
  nHaps = library%getSize()

  if (FullFileOutput == 1) then
    if (WindowsLinux == 1) then
      write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".txt")') currentcore
      open (unit = 24, FILE = filout, status = 'unknown')
    else
      write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".txt")') currentcore
      open (unit = 24, FILE = filout, status = 'unknown')
    endif
  endif
  if (WindowsLinux == 1) then
    write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".bin")') currentcore
    open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
  else
    write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') currentcore
    open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
  endif


  write (34) nHaps, SizeCore
  do i = 1, nHaps
    if (FullFileOutput == 1)&
    write (24, '(2i6,a2,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)') &
    i, library%getHapFreq(i), " ", library%getHap(i)
    write (34) library%getHap(i)
  end do

  if (FullFileOutput == 1) then
    close(24)
  end if
  close(34)
  
  print*, "   ", "Final iteration found ", nHaps, "haplotypes"

  print*, ""
  write (*, '(a4,a30,f5.2,a1)') "   ", "Final yield for this core was ", c%getTotalYield(), "%"

  if (WindowsLinux == 1) then
    open (unit = 29, file = ".\PhasingResults\PhasingYield.txt", status = "unknown", position = "append")
  else
    open (unit = 29, file = "./PhasingResults/PhasingYield.txt", status = "unknown", position = "append")
  endif
  
  write (29, '(i10,f7.2)') CurrentCore, c%getTotalYield()
  
  close(29)

end subroutine WriteHapLib

function fullyPhased(haplotype) result (fully)
  integer(kind=1), dimension(:), intent(in) :: haplotype
  logical :: fully
  
  fully = all((haplotype == 0) .or. (haplotype == 1))
end function fullyPhased

function complement(genos, haplotype) result (comp)
  use Constants
  integer(kind=1), dimension(:), intent(in) :: genos, haplotype
  integer(kind=1), dimension(:), pointer :: comp
  
  integer :: i
  
  allocate(comp(size(genos,1)))
  
  do i = 1, size(genos,1)
    if (genos(i) == MissingGenotypeCode) then
      comp(i) = 9
    else
      comp(i) = genos(i) - haplotype(i)
    end if
  end do
end function complement

function agree(candidates, library, position) result (a)
  integer, dimension(:), intent(in) :: candidates
  class(HapLib), intent(in) :: library
  integer, intent(in) :: position
  
  logical :: a
  integer(kind = 1), dimension(size(candidates,1)) :: phases
  integer :: i
  
  do i = 1, size(candidates,1)
    phases(i) = library%getPhase(candidates(i),position)
  end do
  
  a = (all(phases == 0) .and. (all(phases == 1)))
end function agree

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
end function uniquehaps

end module HaplotypeLibrary
