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
    procedure, public :: addHap
    procedure, public :: getHap
    procedure, public :: matchWithError
    procedure, public :: getPhase
    procedure, public :: getSize
    procedure, public :: getHapRel
    procedure, public :: getNumSnps
    procedure, public :: resetHapFreq
    procedure, public :: incrementHapFreq
    procedure, public :: getHapFreq
    
    procedure, public :: matchAddHap
    
    procedure :: addHap2
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

  subroutine addHap(library, haplotype)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype

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
  end subroutine addHap
  
  !!! Temporary Hack !!!
  function addHap2(library, haplotype) result(id)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    
    integer :: id
    
    call library%addHap(haplotype)
    library%hapfreq(library%size) = 1
    id = library%size
  end function addHap2
  
  function matchAddHap(library, haplotype) result (id)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer :: id
    
    logical :: match
    integer :: i, j
    
    id = 0
    do i = 1, library%size
      match = .true.
      do j = 1, library%nSnps
	if (library%store(i, library%randomOrder(j)) /= haplotype(library%randomOrder(j))) then
	  match = .false.
	  exit
	end if
      end do
      if (match) then
	id = i
	library%hapfreq(i) = library%hapfreq(i) + 1
	exit
      end if
    end do
    
    if (id == 0) then
      id = library%addHap2(haplotype)
    end if    
  end function matchAddHap

  function matchWithError(library, haplotype, allowedError) result(matches)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    integer, intent(in) :: allowedError
    integer, dimension(:), allocatable :: matches

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
	exit
      end if
    end do
    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function matchWithError

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
  
  !MESSY! SHOULDN'T BE HERE!!!!
  !subroutine CheckCompatHapGeno(genos, phase)
  subroutine CheckCompatHapGeno(c)
  use Parameters, only: PercGenoHaploDisagree
  use Constants
  use CoreSubsetDefinition
  implicit none

  class(CoreSubset) :: c

  integer(kind=1), dimension(:,:), pointer :: genos
  
  integer :: i, j, CountError, SizeCore, ErrorAllow, Disagree, counter, counterMissing, nAnisG, nCoreSnp
  double precision :: value, Yield

  nAnisG  = c%getNAnisG()
  nCoreSnp = c%getNCoreSnp()
  
  genos => c%getCoreGenos()
  
  !Refactor out!
  SizeCore = nCoreSnp
  
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
  
  integer :: i, j, truth
  integer :: nSNPcore, nCount, nAnisG
  integer, allocatable :: Shuffle(:)
  integer :: secs
  
  integer :: id

  nAnisG = c%getNAnisG()

  call library%resetHapFreq()
  call c%resetFullyPhased()
  call c%resetHapAnis()

  call system_clock(nCount)
  secs = mod(nCount, int(1e6))

  ! Create random indexes
  nSnpCore = c%getNCoreSnp() ! Total number of markers in the core
  
  !THIS IS HORRIBLE!
  call library%initalise(nSNPcore,500,500)
  
  allocate(Shuffle(nSNPcore))
  call RandomOrder(Shuffle, nSNPcore, 1, -abs(secs))

  do i = 1, nAnisG
    !Paternal Haps
    truth = 0
    do j = 1, nSNPcore
      if (c%getPhase(i, Shuffle(j), 1) == 9) then
	truth = 1
	exit
      endif
    enddo
    if (truth == 0) then
      id = library%matchAddHap(c%getHaplotype(i,1))
      call c%setHapAnis(i,1,id)
      call c%setFullyPhased(i,1)
    endif
    
    !Maternal Haps
    truth = 0
    do j = 1, nSNPcore
      if (c%getPhase(i, Shuffle(j), 2) == 9) then
	truth = 1
	exit
      endif
    enddo
    if (truth == 0) then
      !FUDGE FOR CONSISTENCY.  There should be no if statement here
      if (.not. consistent .or. (library%getSize() > 0)) then
	call c%setFullyPhased(i,2)
      end if
      id = library%matchAddHap(c%getHaplotype(i,2))
      call c%setHapAnis(i,2,id)
    endif
  enddo

  deallocate(Shuffle)
  
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
  integer :: i, j, k, l, m, truth, truth1, HapLibIter, nHapsOld, Disagree, SizeCore, ErrorAllow, HapM, HapP, nCand, nCandPat, Miss, nHapsTmp
  integer :: CompatPairs, WorkScaler, CountZero, CountOne, Switch
  integer(kind=1) :: value
  integer :: id
  integer, allocatable, dimension(:) :: CandGenos, CandHaps, WorkVec
  integer, allocatable, dimension(:,:) :: CandPairs
  character(len = 300) :: filout

  integer :: nSNPcore, nCount, nAnisG, nSnp, secs, nHaps
  integer, allocatable :: Shuffle(:)
  
  integer, allocatable, dimension(:) :: compatHaps
  integer :: numCompatHaps
  
  integer, allocatable, dimension (:,:) :: TempHapArray
  integer, allocatable, dimension (:) :: TempHapVector, ClusterMember
  integer :: nHapsCluster, rounds
  
  integer :: countA, countB, ErrorCountAB
  
  nAnisG = c%getNAnisG()
  nSnp = c%getNCoreSnp()
  nHaps = library%getSize()


  allocate(CandGenos(nSnp))
  allocate(CandHaps(nAnisG * 2))
  allocate(WorkVec(nAnisG * 2))
  allocate(CandPairs(nAnisG * 2, 2))
  
  genos => c%getCoreGenos()

  ! Create a seed for RNG
  call system_clock(nCount)
  secs = mod(nCount, int(1e6))

  ! Create random indexes
  nSnpCore = nSnp ! Total number of markers in the core
  allocate(Shuffle(nSNPcore))
  call RandomOrder(Shuffle, nSNPcore, 1, -abs(secs))

  SizeCore = nSnp
  ErrorAllow = int(PercGenoHaploDisagree * SizeCore)
  ErrorCountAB = int(SizeCore * 0.09)

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
	  CandHaps = 0
	  nCand = 0
	  truth1 = 0
	  HapM = 0
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    do j = 1, nSNPcore
	      if (Genos(i, Shuffle(j)) /= MissingGenotypeCode) then
		if (library%getPhase(k, Shuffle(j)) + c%getPhase(i, Shuffle(j), 1) /= Genos(i, Shuffle(j))) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    truth = 1
		    exit
		  endif
		end if
	      end if
	    end do
	    if (truth == 0) then
	      HapM = k
	      nCand = nCand + 1
	      CandHaps(nCand) = k
	    end if
	  end do

	  if (nCand > 1) then
	    do j = 1, nSNPcore
	      Disagree = 0
	      do k = 1, nCand
		do l = k, nCand
		  if (library%getPhase(CandHaps(k), Shuffle(j)) /= library%getPhase(CandHaps(l), Shuffle(j))) then
		    Disagree = 1
		  end if
		end do
	      end do
	      if (Disagree == 0) then
		call c%setPhase(i, Shuffle(j), 2, library%getPhase(CandHaps(1), Shuffle(j)))
	      end if
	    end do
	  endif

	  if (nCand == 1) then
	    call c%setHaplotype(i,2,library%getHap(HapM))
	    call c%setFullyPhased(i,2)
	    call library%incrementHapFreq(HapM)
	    call c%setHapAnis(i, 2, HapM)
	  end if

	  if (nCand == 0) then
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - c%getPhase(i, j, 1)
		if ((value == 0).or.(value == 1)) then
		  call c%setPhase(i, j, 2, value)
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      call c%setFullyPhased(i,2)
	      truth = 0
	      do k = 1, nHaps
		Disagree = 0
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= c%getPhase(i, Shuffle(j), 2)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  call library%incrementHapFreq(k)
		  call c%setHapAnis(i, 2, k)
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		call library%addHap(c%getHaplotype(i, 2))
		call library%incrementHapFreq(nHaps)
		call c%setHapAnis(i, 2, nHaps)
	      end if
	    end if
	  endif
	end if

	! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): MATERNAL HAPLOTYPE
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Haplotype 2 can get fully phased above and this will run despite both haplotypes now being phased
	! Affects results - not entirely sure why...
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (c%getFullyPhased(i, 2)) then
	  truth1 = 0
	  HapP = 0
	  CandHaps = 0
	  nCand = 0
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    do j = 1, nSNPcore
	      if (Genos(i, Shuffle(j)) /= MissingGenotypeCode) then
		if (library%getPhase(k, Shuffle(j)) + c%getPhase(i, Shuffle(j), 2) /= Genos(i, Shuffle(j))) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    truth = 1
		    exit
		  endif
		end if
	      end if
	    enddo
	    if (truth == 0) then
	      HapP = k
	      nCand = nCand + 1
	      CandHaps(nCand) = k
	    endif
	  enddo

	  if (nCand > 1) then
	    do j = 1, nSNPcore
	      Disagree = 0
	      do k = 1, nCand
		do l = k, nCand
		  if (library%getPhase(CandHaps(k), Shuffle(j)) /= library%getPhase(CandHaps(l), Shuffle(j))) then
		    Disagree = 1
		  end if
		end do
	      end do
	      if (Disagree == 0) then
		call c%setPhase(i, Shuffle(j), 1, library%getPhase(CandHaps(1), Shuffle(j)))
	      end if
	    enddo
	  endif

	  if (nCand == 1) then
	    call c%setHaplotype(i, 1, library%getHap(HapP))
	    call c%setFullyPhased(i, 1)
	    call library%incrementHapFreq(HapP)
	    call c%setHapAnis(i, 1, HapP)
	  endif

	  if (nCand == 0) then
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - c%getPhase(i, j, 2)
		if ((value == 0).or.(value == 1)) then
		  call c%setPhase(i, j, 1, value)
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      call c%setFullyPhased(i, 1)
	      truth = 0
	      do k = 1, nHaps
		Disagree = 0
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= c%getPhase(i, Shuffle(j), 1)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  call library%incrementHapFreq(k)
		  call c%setHapAnis(i, 1, k)
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		call library%addHap(c%getHaplotype(i, 1))
		call library%incrementHapFreq(nHaps)
		call c%setHapAnis(i, 1, nHaps)
	      end if
	    endif
	  endif
	end if

	! If neither of the gametes is completely phased (Section Step 2e.ii Hickey et al. 2011)
	if ((.not. c%getFullyPhased(i,1)) .and. (.not. c%getFullyPhased(i,2))) then
	  HapP = 0
	  HapM = 0
	  CandHaps = 0
	  nCand = 0
	  
	  allocate(compatHaps(nHaps))
	  numCompatHaps = 0
	  do k = 1, nHaps
	    if (consistent) then
	      numCompatHaps = numCompatHaps + 1
	      compatHaps(numCompatHaps) = k
	    else	      
	      Disagree = 0
	      do j = 1, nSNPcore
		if ((library%getPhase(k, Shuffle(j)) == 0) .and. (genos(i, Shuffle(j)) == 2)) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    exit
		  end if
		end if
		if ((library%getPhase(k, Shuffle(j)) == 1) .and. (genos(i, Shuffle(j)) == 0)) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    exit
		  end if
		end if
	      end do
	      if (Disagree <= ErrorAllow) then
		numCompatHaps = numCompatHaps + 1
		compatHaps(numCompatHaps) = k
	      end if
	    end if
	  end do

	  do k = 1, numCompatHaps
	    truth = 0
	    Disagree = 0
	    do j = 1, nSNPcore
	      if ((c%getPhase(i, Shuffle(j), 1) /= 9).and.(c%getPhase(i, Shuffle(j), 1) /= library%getPhase(compatHaps(k), Shuffle(j)))) Disagree = Disagree + 1
	      if (Disagree > ErrorAllow) then
		truth = 1
		exit
	      end if
	    end do

	    ! If the there is no disagreement, we've found a new candidate haplotype
	    if (truth == 0) then
	      nCand = nCand + 1
	      CandHaps(nCand) = compatHaps(k)
	    end if
	  enddo
	  
	  ! Update the number of candidates for paternal haplotype
	  nCandPat = nCand

	  ! If only have one paternal candidate haplotype, then
	  ! the paternal haplotype is nCand
	  if (nCand == 1) HapP = CandHaps(nCand)
	  
	  ! Find candidates for maternal haplotype
	  do k = 1, numCompatHaps
	    truth = 0
	    Disagree = 0
	    do j = 1, nSNPcore
	      if ((c%getPhase(i, Shuffle(j), 2) /= 9).and.(c%getPhase(i, Shuffle(j), 2) /= library%getPhase(compatHaps(k), Shuffle(j)))) Disagree = Disagree + 1
	      if (Disagree > ErrorAllow) then
		truth = 1
		exit
	      end if
	    end do

	    ! If the there is no disagreement, we've found a new candidate haplotype
	    if (truth == 0) then
	      truth1 = 0
	      do j = 1, nCand
		if (k == CandHaps(j)) then
		  truth1 = 1
		  exit
		end if
	      enddo
	      if (truth1 == 0) then
		nCand = nCand + 1
		CandHaps(nCand) = compatHaps(k)
	      endif
	    endif
	  enddo

	  deallocate(compatHaps)
	  
	  ! If only have one maternal candidate haplotype, then
	  ! the maternal haplotype is nCand
	  if ((nCand - nCandPat) == 1) HapM = CandHaps(nCand)
	  
	  ! If only one maternal candidate haplotype and many paternal candidate haplotypes
	  if ((HapM > 0).AND.(HapP == 0)) then
	    truth1 = 0
	    do k = 1, nCandPat
	      Disagree = 0
	      truth = 1
	      do j = 1, nSNPcore
		if ((Genos(i, Shuffle(j)) /= MissingGenotypeCode).and.&
		  (Genos(i, Shuffle(j)) /= (library%getPhase(HapM, Shuffle(j)) + library%getPhase(CandHaps(k), Shuffle(j))))) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    truth = 0
		    exit
		  end if
		endif
	      enddo
	      if (truth == 1) then
		truth1 = truth1 + 1
		HapP = CandHaps(k)
		if (truth1 > 1) then
		  HapP = 0
		end if
	      end if
	    end do
	  end if

	  ! If only one paternal candidate haplotype and one / many maternal candidate haplotypes 
	  if ((nCandPat == 1).and.(nCand - nCandPat > 0)) then
	    truth1 = 0
	    do k = nCandPat + 1, nCand
	      Disagree = 0
	      truth = 1
	      do j = 1, nSNPcore
		if ((Genos(i, Shuffle(j)) /= MissingGenotypeCode).and.(Genos(i, Shuffle(j)) /= (library%getPhase(HapP, Shuffle(j)) &
		  + library%getPhase(CandHaps(k),Shuffle(j))))) then                                                             
		Disagree = Disagree + 1
		if (Disagree > ErrorAllow) then
		  truth = 0
		  exit
		end if
	      endif
	    enddo
	    ! Really not sure about this next if.  Seems to me to be saying that if we've found one compatible hap but then find
	    ! one that isn't we set HapM to none - despite so far only having one match...  Also why the difference from above?
	    ! Think this is here so that if we only have one candidate for maternal it's compatible but think it has undesirable
	    ! side effects
	    if (truth == 0) HapM = 0
	    if (truth == 1) then
	      truth1 = truth1 + 1
	      HapM = CandHaps(k)
	      if (truth1 > 1) then
		HapM = 0
	      end if
	    end if
	  end do
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
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - c%getPhase(i, j, 1)
		if ((value == 0).or.(value == 1)) then
		  call c%setPhase(i, j, 2, value)
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      call c%setFullyPhased(i, 2)
	      truth = 0

	      ! Update (if necessary) Haplotype Library with the new maternal gamete found
	      do k = 1, nHaps
		Disagree = 0
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= c%getPhase(i, Shuffle(j), 2)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  call library%incrementHapFreq(k)
		  call c%setHapAnis(i, 2, k)
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		call library%addHap(c%getHaplotype(i, 2))
		call library%incrementHapFreq(nHaps)
		call c%setHapAnis(i, 2, nHaps)
	      end if
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
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - c%getPhase(i, j, 2)
		if ((value == 0).or.(value == 1)) then
		  call c%setPhase(i, j, 1, value)
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      call c%setFullyPhased(i, 1)
	      truth = 0

	      ! Update (if necessary) Haplotype Library with the new maternal gamete found
	      do k = 1, nHaps
		Disagree = 0
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= c%getPhase(i, Shuffle(j), 1)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  call library%incrementHapFreq(k)
		  call c%setHapAnis(i, 1, k)
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		call library%addHap(c%getHaplotype(i, 1))
		call library%incrementHapFreq(nHaps)
		call c%setHapAnis(i, 1, nHaps)
	      end if
	    end if
	  end if
	end if

	! If the paternal and maternal gamete cannot be identify without ambiguity
	! (more than one or none at all)           
	if ((HapP == 0).and.(HapM == 0)) then
	  CandPairs = 0
	  CompatPairs = 0
	  do k = 1, nCand
	    do l = (k + 1), nCand
	      Disagree = 0
	      truth = 1

	      ! Check agreement between pairs
	      do j = 1, nSNPcore
		if (Genos(i, Shuffle(j)) /= MissingGenotypeCode)then
		  if ((library%getPhase(CandHaps(k), Shuffle(j)) + library%getPhase(CandHaps(l), Shuffle(j))) /= Genos(i, Shuffle(j))) then
		    Disagree = Disagree + 1
		    if (Disagree > ErrorAllow) then
		      truth = 0
		      exit
		    end if
		  end if
		endif
	      end do

	      ! If there is a pair that agrees, 
	      ! the haplotypes are consider to be the paternal and maternal gametes (arbitrarily)
	      if (truth == 1) then
		CompatPairs = CompatPairs + 1
		CandPairs(CompatPairs, 1) = CandHaps(k)
		CandPairs(CompatPairs, 2) = CandHaps(l)
		HapP = CandHaps(k)
		HapM = CandHaps(l)
		if ((CompatPairs * CompatPairs) > (nAnisG - 1)) exit
	      end if
	    end do
	    if ((CompatPairs * CompatPairs) > (nAnisG - 1)) exit
	  end do

	  ! If only one pair agrees...
	  if (CompatPairs == 1) then
	    ! Phase the paternal haplotype and update the library with the new frequency 
	    call c%setHaplotype(i, 1, library%getHap(HapP))
	    call c%setFullyPhased(i, 1)
	    call library%incrementHapFreq(HapP)
	    call c%setHapAnis(i, 1, HapP)
	    ! Phase the maternal haplotype and update the library with the new frequency 
	    call c%setHaplotype(i, 2, library%getHap(HapM))
	    call c%setFullyPhased(i, 2)
	    call library%incrementHapFreq(HapM)
	    call c%setHapAnis(i, 2, HapM)
	  end if

	  ! If more than one pair agrees...                
	  if ((CompatPairs > 1).and.((CompatPairs * CompatPairs) < nAnisG)) then !Note the 200 number is a fudge
	    truth = 1

	    ! Check how many paternal candidates haplotypes
	    id = CandPairs(1, 1)
	    do k = 2, CompatPairs
	      if (CandPairs(k, 1) /= id) truth = 0
	    end do
	    Switch = 0

	    ! If there is only one paternal haplotype in all the candidate pairs
	    if (truth == 1) then
	      ! Phase the paternal gamete with this haplotype
	      call c%setHaplotype(i, 1, library%getHap(id))
	      call c%setFullyPhased(i, 1)
	      call library%incrementHapFreq(id)
	      call c%setHapAnis(i, 1, id)

	      ! If only one haplotype is found for the paternal gamete 
	      ! and many for the maternal gamete, phase each loci only all pairs agree
	      ! (Step 2e.ii.B)
	      do j = 1, nSNPcore
		value = library%getPhase(CandPairs(1, 2), Shuffle(j))
		truth1 = 1
		do k = 2, CompatPairs
		  if (library%getPhase(CandPairs(k, 2), Shuffle(j)) /= value) then
		    truth1 = 0
		    exit
		  end if
		end do
		if (truth1 == 1) then
		  call c%setPhase(i, Shuffle(j), 2, value)
		end if
	      end do
	      Switch = 1
	    end if

	    ! Check how many maternal candidates haplotypes
	    truth = 1
	    id = CandPairs(1, 2)
	    do k = 2, CompatPairs
	      if (CandPairs(k, 2) /= id) truth = 0
	    end do

	    ! If there is only one maternal haplotype in all the candidate pairs
	    if (truth == 1) then
	      call c%setHaplotype(i, 2, library%getHap(id))
	      call c%setFullyPhased(i, 2)
	      call library%incrementHapFreq(id)
	      call c%setHapAnis(i, 2, id)

	      ! If only one haplotype is found for the paternal gamete 
	      ! and many for the maternal gamete, phase each loci only all pairs agree
	      ! (Step 2e.ii.C)
	      do j = 1, nSNPcore
		value = library%getPhase(CandPairs(1, 1), Shuffle(j))
		truth1 = 1
		do k = 2, CompatPairs
		  if (library%getPhase(CandPairs(k, 1), Shuffle(j)) /= value) then
		    truth1 = 0
		    exit
		  end if
		end do
		if (truth1 == 1) then
		  call c%setPhase(i, Shuffle(j), 1, value)
		endif
	      enddo
	      Switch = 1
	    endif

	    ! If proband is not completely phased and have more than one candidate 
	    ! for both paternal and maternal haplotype 
	    ! (Step 2e.iv)
	    if ( ((.not. c%getFullyPhased(i,1)) .or. (.not. (c%getFullyPhased(i,2)))) .and. (Switch == 0)) then

	      ! Initialize procedure of k-medoids
	      WorkVec = 0
	      do k = 1, CompatPairs
		WorkVec(CandPairs(k, 1)) = 1
		WorkVec(CandPairs(k, 2)) = 1
	      enddo
	      WorkScaler = sum(WorkVec(:))
	      allocate(TempHapArray(WorkScaler, SizeCore))
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
		  do j = 1, nSnpCore
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
		  end do
		  do j = 1, nSnpCore
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
  do j = 1, nSnpCore
    if (Genos(i, j) /= MissingGenotypeCode) then
      !if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) == 9)) then
      if ((c%getPhase(i, j, 1) /= 9).and.(c%getPhase(i, j, 2) == 9)) then
	!value = Genos(i, j) - Phase(i, j, 1)
	value = Genos(i, j) - c%getPhase(i, j, 1)
	if ((value == 0).or.(value == 1)) then !here 7th april 2011
	  !Phase(i, j, 2) = value
	  call c%setPhase(i, j, 2, value)
	else
	  CountA = CountA + 1
	endif
      endif
      !if ((Phase(i, j, 2) /= 9).and.(Phase(i, j, 1) == 9)) then
      if ((c%getPhase(i, j, 2) /= 9).and.(c%getPhase(i, j, 1) == 9)) then
	!value = Genos(i, j) - Phase(i, j, 2)
	value = Genos(i, j) - c%getPhase(i, j, 2)
	if ((value == 0).or.(value == 1)) then !here 7th april 2011
	  !Phase(i, j, 1) = value
	  call c%setPhase(i, j, 1, value)
	else
	  CountB = CountB + 1
	endif
      endif
    end if
  end do
  if ((CountA > ErrorCountAB) .or. (CountB > ErrorCountAB)) then
    !Phase(i, :, :) = 9
    call c%setHaplotypeToUnphased(i,1)
    call c%setHaplotypeToUnphased(i,2)
    do j = 1, nSnpCore
      !if (Genos(i, j) == 0) Phase(i, j,:) = 0
      !if (Genos(i, j) == 2) Phase(i, j,:) = 1
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
!!  if (CountB > ErrorCountAB) then
!!    Phase(i, : ,:) = 9
!!    do j = 1, nSnpCore
!!      if (Genos(i, j) == 0) Phase(i, j,:) = 0
!!      if (Genos(i, j) == 2) Phase(i, j,:) = 1
!!    enddo
!!  endif
end do

do i = 1, nAnisG
  do j = 1, nSnpCore
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

  ! This should be refactored out at some point - same as currentcore
  integer :: OutputPoint
  
  integer :: i, j, k, counter, SizeCore, nHaps !, nAnisG
  character(len = 300) :: filout
  
  SizeCore = library%getNumSnps()
  OutputPoint = currentcore
  
  nHaps = library%getSize()

  if (FullFileOutput == 1) then
    if (WindowsLinux == 1) then
      write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".txt")') OutputPoint
      open (unit = 24, FILE = filout, status = 'unknown')
    else
      write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".txt")') OutputPoint
      open (unit = 24, FILE = filout, status = 'unknown')
    endif
  endif
  if (WindowsLinux == 1) then
    write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".bin")') OutputPoint
    open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
  else
    write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') OutputPoint
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

end module HaplotypeLibrary
