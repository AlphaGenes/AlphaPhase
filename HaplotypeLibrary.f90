module HaplotypeLibrary
  implicit none
  public

  type, public :: HapLib
    private
    !integer(kind = 1), dimension (:,:), pointer :: store => null()
    integer(kind = 1), dimension (:,:), allocatable :: store
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
    end if
    allocate(library % store(storeSize, nSnps))
    library % store = 0
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

  !function addHap(library, haplotype) result(id)
  subroutine addHap(library, haplotype)
    class(HapLib) :: library
    integer(kind = 1), dimension(:), intent(in) :: haplotype
    !integer :: id

    integer :: newStoreSize
!    integer(kind = 1), dimension(:,:), allocatable, target :: newStore
    integer(kind = 1), dimension(:,:), allocatable :: tempStore

    if (library % Size == library % storeSize) then
!      newStoreSize = library % storeSize + library % stepSize
!      allocate(newStore(newStoreSize, library % nSnps))
!      newStore = 0
!      newStore(1:library % Size,:) = library % Store
!      !deallocate(library%Store)
!      library % Store => newStore
!      library % StoreSize = newStoreSize
      newStoreSize = library % storeSize + library % stepSize
      allocate(tempStore(library % storeSize, library % nSnps))
      tempStore = library%store
      deallocate(library%store)
      allocate(library%store(newStoreSize, library % nSnps))
      library % store = 0
      library % store(1:library % Size,:) = tempStore
      deallocate(tempStore)
      library % StoreSize = newStoreSize
    end if

    library % Size = library % Size + 1
    library % Store(library % Size,:) = haplotype
    !id = library%Size

    !end function addHap
  end subroutine addHap

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

    !integer(kind = 1), dimension(:), allocatable :: hap

    !hap = library % getHap(id)
    !phase = hap(snp)
    !if ((id == 13) .and. (snp == 1)) then
    !  print *, "ERR:", library % store(id,snp)
    !end if
    
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
	!do k = StartCoreSnp, EndCoreSnp
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
  
  !MESSY!
  subroutine CheckCompatHapGeno(genos, phase)
  use Global, only: MissingGenotypeCode, PercGenoHaploDisagree
  implicit none

  integer(kind=1), dimension(:,:), intent(in) :: genos
  integer(kind=1), dimension(:,:,:), intent(inout) :: phase
  
  integer :: i, j, CountError, SizeCore, ErrorAllow, Disagree, counter, counterMissing, nAnisG, nCoreSnp
  double precision :: value, Yield
  
  nAnisG = size(genos,1)
  nCoreSnp = size(genos,2)
  
  !Refactor out!
  SizeCore = nCoreSnp
  
   ErrorAllow = int(PercGenoHaploDisagree * nCoreSnp)

  do i = 1, nAnisG
    CountError = 0
    counterMissing = 0
    do j = 1, nCoreSnp
      if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) /= 9)) then
	counterMissing = counterMissing + 1
	if ((Genos(i, j) /= MissingGenotypeCode).and.(sum(Phase(i, j,:)) /= Genos(i, j))) CountError = CountError + 1
      end if
    end do
    ErrorAllow = int(PercGenoHaploDisagree * counterMissing)
    if (CountError >= ErrorAllow) then
      do j = 1, nCoreSnp
	if (Genos(i, j) /= MissingGenotypeCode) then
	  if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) /= 9).and.(sum(Phase(i, j,:)) /= Genos(i, j))) then
	    if (Genos(i, j) == 1) Phase(i, j, 2) = 9
	    if (Genos(i, j) == MissingGenotypeCode) Phase(i, j, 2) = 9
	    if (Genos(i, j) == 0) Phase(i, j,:) = 0
	    if (Genos(i, j) == 2) Phase(i, j,:) = 1
	  endif
	endif
      enddo
    endif
  end do

  counter = count(Phase(:, :, 1) == 0)
  counter = count(Phase(:, :, 1) == 1) + counter
  Yield = (float(counter)/(nAnisG * SizeCore)) * 100
  print*, " "
  write (*, '(a3,f6.2,a45)') "  ", Yield, "% was the Paternal allele yield for this core"
  counter = count(Phase(:, :, 2) == 0)
  counter = count(Phase(:, :, 2) == 1) + counter
  Yield = (float(counter)/(nAnisG * SizeCore)) * 100
  write (*, '(a3,f6.2,a45)') "  ", Yield, "% was the Maternal allele yield for this core"

end subroutine CheckCompatHapGeno

subroutine MakeHapLib(library, phase, fullyPhased, hapfreq, hapanis)
  use Global, only: nGlobalHaps
  use GlobalClusteringHaps
  implicit none

  type(HapLib), intent(in) :: library
  integer(kind=1), dimension(:,:,:), intent(in) :: phase
  logical, dimension(:,:), intent(inout) :: fullyPhased
  integer, dimension(:), intent(inout) :: hapfreq
  integer, dimension(:,:), intent(inout) :: hapanis
  
  integer :: i, j, k, l, m, truth, truth1
  integer :: nSNPcore, nCount, nAnisG
  integer, allocatable :: Shuffle(:)
  integer :: nHaps
  integer :: secs


  INTERFACE
    subroutine RandomOrder(order, n, start, idum)
      !     Generate a random ordering of the integers 1 ... n.

      integer, INTENT(IN) :: n, start
      integer, allocatable, INTENT(OUT) :: order(:)
    end subroutine RandomOrder
  END INTERFACE
  
  nAnisG = size(phase,1)

  nHaps = 0
  HapFreq = 0
  FullyPhased = .false.
  HapAnis = -99

  ! Create a seed for RNG
  call system_clock(nCount)
  secs = mod(nCount, int(1e6))

  ! Create random indexes
  nSnpCore = size(phase,2) ! Total number of markers in the core

  !THIS IS HORRIBLE!
  call library%initalise(nSNPcore,500,500)
  
  allocate(Shuffle(nSNPcore))
  call RandomOrder(Shuffle, nSNPcore, 1, -abs(secs))

  do i = 1, nAnisG
    !Paternal Haps
    truth = 0
    ! do j=StartCoreSnp,EndCoreSnp
    ! if (Phase(i,j,1)==9) then
    do j = 1, nSNPcore
      if (Phase(i, Shuffle(j), 1) == 9) then
	truth = 1
	exit
      endif
    enddo
    if (truth == 0) then
      if (nHaps == 0) then
	nHaps = 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	HapAnis(i, 1) = nHaps
	call library%addHap(Phase(i, :, 1))
      else
	do k = 1, nHaps
	  Truth1 = 0
	  ! do j=StartCoreSnp,EndCoreSnp
	  !         if (HapLib(k,j)/=Phase(i,j,1)) then
	  do j = 1, nSNPcore
	    if (library%getPhase(k, Shuffle(j)) /= Phase(i, Shuffle(j), 1)) then
	      Truth1 = 1
	      exit
	    end if
	  end do
	  if (Truth1 == 0) then
	    HapFreq(k) = HapFreq(k) + 1
	    HapAnis(i, 1) = k
	    exit
	  end if
	end do
	if (Truth1 == 1) then
	nHaps = nHaps + 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	call library%addHap(Phase(i, :, 1))
	HapAnis(i, 1) = nHaps
	end if
      end if
      FullyPhased(i, 1) = .true.
    endif
    !Maternal Haps
    truth = 0
    ! do j=StartCoreSnp,EndCoreSnp
    !         if (Phase(i,j,2)==9) then
    do j = 1, nSNPcore
      if (Phase(i, Shuffle(j), 2) == 9) then
	truth = 1
	exit
      endif
    enddo
    if (truth == 0) then
      if (nHaps == 0) then
	nHaps = 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	HapAnis(i, 2) = nHaps
	call library%addHap(Phase(i, :, 2))
      else
	do k = 1, nHaps
	  Truth1 = 0
	  ! do j=StartCoreSnp,EndCoreSnp
	  !         if (HapLib(k,j)/=Phase(i,j,2)) then
	  do j = 1, nSNPcore
	    if (library%getPhase(k, Shuffle(j)) /= Phase(i, Shuffle(j), 2)) then
	      Truth1 = 1
	      exit
	    end if
	  end do
	  if (Truth1 == 0) then
	    HapFreq(k) = HapFreq(k) + 1
	    HapAnis(i, 2) = k
	    exit
	  end if
	end do
	if (Truth1 == 1) then
	nHaps = nHaps + 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	call library%addHap(Phase(i, :, 2))
	HapAnis(i, 2) = nHaps
	end if
	FullyPhased(i, 2) = .true.
      endif
    endif
  enddo
  nGlobalHaps = nHaps

  deallocate(Shuffle)
  
end subroutine MakeHapLib

subroutine ImputeFromLib(library, genos, phase, fullyphased, hapfreq, hapanis)
  ! Impute the phase for gametes that are not completely phased by LRP 
  ! by matching their phased loci to haplotypes in the Haplotype Library,
  ! following strategies listed in the section Step 2e of Hickey et al 2011.

  use Global, only : percgenohaplodisagree, missinggenotypecode, nglobalhapsiter, nglobalhaps
  use GlobalClusteringHaps
  implicit none
  
  type(HapLib), intent(in) :: library
  integer(kind=1), dimension(:,:), intent(in) :: genos
  integer(kind=1), dimension(:,:,:), intent(inout) :: phase
  logical, dimension(:,:) :: fullyphased
  integer, dimension(:), intent(inout) :: hapfreq
  integer, dimension(:,:), intent(inout) ::hapanis

  integer :: i, j, k, l, m, truth, truth1, HapLibIter, nHapsOld, Disagree, SizeCore, ErrorAllow, HapM, HapP, nCand, nCandPat, Miss, nHapsTmp
  integer :: CompatPairs, value, WorkScaler, CountZero, CountOne, Switch
  integer :: CountA, CountB, ErrorCountAB
  integer, allocatable, dimension(:) :: CandGenos, CandHaps, WorkVec!,ErrorAllow
  integer, allocatable, dimension(:,:) :: CandPairs
  character(len = 300) :: filout

  integer :: nSNPcore, nCount, nAnisG, nSnp, secs, nHaps
  integer, allocatable :: Shuffle(:)
  
  
  INTERFACE
    subroutine RandomOrder(order, n, start, idum)
      !     Generate a random ordering of the integers 1 ... n.

      integer, INTENT(IN) :: n, start
      integer, allocatable, INTENT(OUT) :: order(:)
    end subroutine RandomOrder
  END INTERFACE

  nAnisG = size(genos,1)
  nSnp = size(genos,2)
  nHaps = library%getSize()


  allocate(CandGenos(nSnp))
  allocate(CandHaps(nAnisG * 2))
  allocate(WorkVec(nAnisG * 2))
  allocate(CandPairs(nAnisG * 2, 2))

  ! Create a seed for RNG
  call system_clock(nCount)
  secs = mod(nCount, int(1e6))

  ! Create random indexes
  nSnpCore = nSnp ! Total number of markers in the core
  allocate(Shuffle(nSNPcore))
  call RandomOrder(Shuffle, nSNPcore, 1, -abs(secs))

  SizeCore = nSnp
  ErrorAllow = int(PercGenoHaploDisagree * SizeCore)
  SnpInCore = SizeCore
  ErrorCountAB = int(SizeCore * 0.09)

  HapLibIter = 1
  nGlobalHaps = nHaps
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
      if ((.not. FullyPhased(i,1)) .or. (.not. FullyPhased(i,2)))  then

	! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): PATERNAL HAPLOTYPE
	if (FullyPhased(i, 1)) then
	  CandHaps = 0
	  nCand = 0
	  truth1 = 0
	  HapM = 0
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if (Genos(i,j)/=MissingGenotypeCode) then
	      ! if (HapLib(k,j)+Phase(i,j,1)/=Genos(i,j)) then
	      if (Genos(i, Shuffle(j)) /= MissingGenotypeCode) then
		if (library%getPhase(k, Shuffle(j)) + Phase(i, Shuffle(j), 1) /= Genos(i, Shuffle(j))) then
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
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      Disagree = 0
	      do k = 1, nCand
		do l = k, nCand
		  ! if (HapLib(CandHaps(k),j)/=HapLib(CandHaps(l),j)) then
		  if (library%getPhase(CandHaps(k), Shuffle(j)) /= library%getPhase(CandHaps(l), Shuffle(j))) then
		    Disagree = 1
		  end if
		end do
	      end do
	      if (Disagree == 0) then
		! Phase(i,j,2)=HapLib(CandHaps(1),j)
		Phase(i, Shuffle(j), 2) = library%getPhase(CandHaps(1), Shuffle(j))
	      end if
	    end do
	  endif

	  if (nCand == 1) then
	    Phase(i, :, 2) = library%getHap(HapM)
	    FullyPhased(i, 2) = .true.
	    HapFreq(HapM) = HapFreq(HapM) + 1
	    HapAnis(i, 2) = HapM
	  end if

	  if (nCand == 0) then
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 1)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 2) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 2) = .true.
	      truth = 0
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		! if (HapLib(k,j)/=Phase(i,j,2)) then
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= Phase(i, Shuffle(j), 2)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 2) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, :, 2))
		HapAnis(i, 2) = nHaps
	      end if
	    end if
	  endif
	end if

	! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): MATERNAL HAPLOTYPE
	if (FullyPhased(i, 2)) then
	  truth1 = 0
	  HapP = 0
	  CandHaps = 0
	  nCand = 0
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if (Genos(i,j)/=MissingGenotypeCode) then
	      !     if (HapLib(k,j)+Phase(i,j,2)/=Genos(i,j)) then
	      if (Genos(i, Shuffle(j)) /= MissingGenotypeCode) then
		if (library%getPhase(k, Shuffle(j)) + Phase(i, Shuffle(j), 2) /= Genos(i, Shuffle(j))) then
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
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      Disagree = 0
	      do k = 1, nCand
		do l = k, nCand
		  ! if (HapLib(CandHaps(k),j)/=HapLib(CandHaps(l),j)) then
		  if (library%getPhase(CandHaps(k), Shuffle(j)) /= library%getPhase(CandHaps(l), Shuffle(j))) then
		    Disagree = 1
		  end if
		end do
	      end do
	      if (Disagree == 0) then
		! Phase(i,j,1)=HapLib(CandHaps(1),j)
		Phase(i, Shuffle(j), 1) = library%getPhase(CandHaps(1), Shuffle(j))
	      end if
	    enddo
	  endif

	  if (nCand == 1) then
	    Phase(i, :, 1) = library%getHap(HapP)
	    FullyPhased(i, 1) = .true.
	    HapFreq(HapP) = HapFreq(HapP) + 1
	    HapAnis(i, 1) = HapP
	  endif

	  if (nCand == 0) then
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 2)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 1) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 1) = .true.
	      truth = 0
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		do j = 1, nSNPcore
		  ! if (HapLib(k,j)/=Phase(i,j,1)) then
		  if (library%getPhase(k, Shuffle(j)) /= Phase(i, Shuffle(j), 1)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 1) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, :, 1))
		HapAnis(i, 1) = nHaps
	      end if
	    endif
	  endif
	end if

	! If neither of the gametes is completely phased (Section Step 2e.ii Hickey et al. 2011)
	if ((.not. FullyPhased(i,1)) .and. (.not. FullyPhased(i,2))) then
	  HapP = 0
	  HapM = 0
	  CandHaps = 0
	  nCand = 0

	  ! Find candidates for paternal haplotype
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if ((Phase(i,j,1)/=9).and.(Phase(i,j,1)/=HapLib(k,j))) Disagree=Disagree+1
	      if ((Phase(i, Shuffle(j), 1) /= 9).and.(Phase(i, Shuffle(j), 1) /= library%getPhase(k, Shuffle(j)))) Disagree = Disagree + 1
	      if (Disagree > ErrorAllow) then
		truth = 1
		exit
	      end if
	    end do

	    ! If the there is no disagreement, we've found a new candidate haplotype
	    if (truth == 0) then
	      nCand = nCand + 1
	      CandHaps(nCand) = k
	    end if
	  enddo

	  ! Update the number of candidates for paternal haplotype
	  nCandPat = nCand

	  ! If only have one paternal candidate haplotype, then
	  ! the paternal haplotype is nCand
	  if (nCand == 1) HapP = CandHaps(nCand)
	  
	  ! Find candidates for maternal haplotype
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if ((Phase(i,j,2)/=9).and.(Phase(i,j,2)/=HapLib(k,j))) Disagree=Disagree+1
	      if ((Phase(i, Shuffle(j), 2) /= 9).and.(Phase(i, Shuffle(j), 2) /= library%getPhase(k, Shuffle(j)))) Disagree = Disagree + 1
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
		CandHaps(nCand) = k
	      endif
	    endif
	  enddo

	  ! If only have one maternal candidate haplotype, then
	  ! the maternal haplotype is nCand
	  if ((nCand - nCandPat) == 1) HapM = CandHaps(nCand)
	  
	  ! If only one maternal candidate haplotype and many paternal candidate haplotypes
	  if ((HapM > 0).AND.(HapP == 0)) then
	    truth1 = 0
	    do k = 1, nCandPat
	      Disagree = 0
	      truth = 1
	      ! do j=StartCoreSnp,EndCoreSnp
	      ! if ((Genos(i,j)/=MissingGenotypeCode).and.&
	      !         (Genos(i,j)/=(HapLib(HapM,j)+HapLib(CandHaps(k),j)))) then
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

	  ! If only one paternal candidate haplotype and many maternal candidate haplotypes 
	  if ((nCandPat == 1).and.(nCand - nCandPat > 0)) then
	    truth1 = 0
	    do k = nCandPat + 1, nCand
	      Disagree = 0
	      truth = 1
	      ! do j=StartCoreSnp,EndCoreSnp
	      ! if ((Genos(i,j)/=MissingGenotypeCode).and.(Genos(i,j)/=(HapLib(HapP,j)+HapLib(CandHaps(k),j)))) then                                                             
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
	  Phase(i, :, 1) = library%getHap(HapP)
	  FullyPhased(i, 1) = .true.

	  ! Update the Library
	  HapFreq(HapP) = HapFreq(HapP) + 1
	  HapAnis(i, 1) = HapP

	  ! If no haplotype has been found for the maternal gamete, or 
	  ! there are more than one maternal candidate haplotype
	  if (HapM == 0) then
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 1)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 2) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 2) = .true.
	      truth = 0

	      ! Update (if necessary) Haplotype Library with the new maternal gamete found
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		! if (HapLib(k,j)/=Phase(i,j,2)) then
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= Phase(i, Shuffle(j), 2)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 2) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, :, 2))
		HapAnis(i, 2) = nHaps
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
	  Phase(i, :, 2) = library%getHap(HapM)
	  FullyPhased(i, 2) = .true.

	  ! Update the Library
	  HapFreq(HapM) = HapFreq(HapM) + 1
	  HapAnis(i, 2) = HapM

	  ! If no haplotype has been found for the paternal gamete, or 
	  ! there are more than one paternal candidate haplotype
	  if (HapP == 0) then
	    Miss = 0
	    do j = 1, nSnpCore
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 2)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 1) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 1) = .true.
	      truth = 0

	      ! Update (if necessary) Haplotype Library with the new maternal gamete found
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		! if (HapLib(k,j)/=Phase(i,j,1)) then
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j)) /= Phase(i, Shuffle(j), 1)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 1) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, :, 1))
		HapAnis(i, 1) = nHaps
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
	      ! do j=StartCoreSnp,EndCoreSnp
	      ! if (Genos(i,j)/=MissingGenotypeCode)then
	      !     if ((HapLib(CandHaps(k),j)+HapLib(CandHaps(l),j))/=Genos(i,j)) then
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
	    Phase(i, :, 1) = library%getHap(HapP)
	    FullyPhased(i, 1) = .true.
	    HapFreq(HapP) = HapFreq(HapP) + 1
	    HapAnis(i, 1) = HapP
	    ! Phase the maternal haplotype and update the library with the new frequency 
	    Phase(i, :, 2) = library%getHap(HapM)
	    FullyPhased(i, 2) = .true.
	    HapFreq(HapM) = HapFreq(HapM) + 1
	    HapAnis(i, 2) = HapM
	  end if

	  ! If more than one pair agrees...                
	  if ((CompatPairs > 1).and.((CompatPairs * CompatPairs) < nAnisG)) then !Note the 200 number is a fudge
	    truth = 1

	    ! Check how many paternal candidates haplotypes
	    value = CandPairs(1, 1)
	    do k = 2, CompatPairs
	      if (CandPairs(k, 1) /= value) truth = 0
	    end do
	    Switch = 0

	    ! If there is only one paternal haplotype in all the candidate pairs
	    if (truth == 1) then
	      ! Phase the paternal gamete with this haplotype
	      Phase(i, :, 1) = library%getHap(value)
	      FullyPhased(i, 1) = .true.
	      HapFreq(value) = HapFreq(value) + 1
	      HapAnis(i, 1) = value

	      ! If only one haplotype is found for the paternal gamete 
	      ! and many for the maternal gamete, phase each loci only all pairs agree
	      ! (Step 2e.ii.B)
	      ! do j=StartCoreSnp,EndCoreSnp
	      do j = 1, nSNPcore
		! value=HapLib(CandPairs(1,2),j)
		value = library%getPhase(CandPairs(1, 2), Shuffle(j))
		truth1 = 1
		do k = 2, CompatPairs
		  ! if (HapLib(CandPairs(k,2),j)/=value) then
		  if (library%getPhase(CandPairs(k, 2), Shuffle(j)) /= value) then
		    truth1 = 0
		    exit
		  end if
		end do
		if (truth1 == 1) then
		  ! Phase(i,j,2)=value
		  Phase(i, Shuffle(j), 2) = value
		end if
	      end do
	      Switch = 1
	    end if

	    ! Check how many maternal candidates haplotypes
	    truth = 1
	    value = CandPairs(1, 2)
	    do k = 2, CompatPairs
	      if (CandPairs(k, 2) /= value) truth = 0
	    end do

	    ! If there is only one maternal haplotype in all the candidate pairs
	    if (truth == 1) then
	      Phase(i, :, 2) = library%getHap(value)
	      FullyPhased(i, 2) = .true.
	      HapFreq(value) = HapFreq(value) + 1
	      HapAnis(i, 2) = value

	      ! If only one haplotype is found for the paternal gamete 
	      ! and many for the maternal gamete, phase each loci only all pairs agree
	      ! (Step 2e.ii.C)
	      ! do j=StartCoreSnp,EndCoreSnp
	      do j = 1, nSNPcore
		! value=HapLib(CandPairs(1,1),j)
		value = library%getPhase(CandPairs(1, 1), Shuffle(j))
		truth1 = 1
		do k = 2, CompatPairs
		  ! if (HapLib(CandPairs(k,1),j)/=value) then
		  if (library%getPhase(CandPairs(k, 1), Shuffle(j)) /= value) then
		    truth1 = 0
		    exit
		  end if
		end do
		if (truth1 == 1) then
		  ! Phase(i,j,1)=value
		  Phase(i, Shuffle(j), 1) = value
		endif
	      enddo
	      Switch = 1
	    endif

	    ! If proband is not completely phased and have more than one candidate 
	    ! for both paternal and maternal haplotype 
	    ! (Step 2e.iv)
	    if ( ((.not. FullyPhased(i,1)) .or. (.not. (FullyPhased(i,2)))) .and. (Switch == 0)) then

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
	      do k = 1, nAnisG * 2
		if (WorkVec(k) == 1) then
		  nHapsCluster = nHapsCluster + 1
		  TempHapVector(nHapsCluster) = k
		  TempHapArray(nHapsCluster, 1:SnpInCore) = &
		  library%getHap(k)
		end if
	      end do
	      allocate(Medoids(nClusters, SnpInCore))
	      allocate(ClusterMember(nHapsCluster))
	      allocate(MinClust(nHapsCluster))
	      do j = 1, nHapsCluster
		if (mod(j, 2) == 0) then
		  ClusterMember(j) = 1
		else
		  ClusterMember(j) = 2
		endif
	      end do
	      call EvaluateMedoidsHaps
	      Change = 0
	      MinClust = 1
	      rounds = 1
	      call RePartitionHaps
	      do j = 1, nMaxRounds
		call EvaluateMedoidsHaps
		Change = 0
		call RePartitionHaps
		if (Change == 0) exit
	      enddo
	      if (rounds <= nMaxRounds) then
		if (count(ClusterMember(:) == 2) == 1) then
		  HapM = 0
		  do j = 1, nHapsCluster
		    if (ClusterMember(j) == 2) then
		      HapM = TempHapVector(j)
		      exit
		    endif
		  enddo
		  Phase(i, :, 2) = library%getHap(HapM)
		  FullyPhased(i, 2) = .true.
		  HapFreq(HapM) = HapFreq(HapM) + 1
		  HapAnis(i, 2) = HapM
		end if
		if (count(ClusterMember(:) == 1) == 1) then
		  HapP = 0
		  do j = 1, nHapsCluster
		    if (ClusterMember(j) == 1) then
		      HapP = TempHapVector(j)
		      exit
		    endif
		  enddo
		  Phase(i, :, 1) = library%getHap(HapP)
		  FullyPhased(i, 1) = .true.
		  HapFreq(HapP) = HapFreq(HapP) + 1
		  HapAnis(i, 1) = HapP
		end if
		if ((count(ClusterMember(:) == 2) > 1).and.(count(ClusterMember(:) == 2) > 1)) then
		  Phase(i, :,:) = 9
		  do j = 1, nSnpCore
		    if (Genos(i, j) == 0) Phase(i, j,:) = 0
		    if (Genos(i, j) == 2) Phase(i, j,:) = 1
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
		    if ((CountZero == 0).and.(CountOne > 0)) Phase(i, j, 2) = 1
		    if ((CountZero > 0).and.(CountOne == 0)) Phase(i, j, 2) = 0
		  end do
		  do j = 1, nSnpCore
		    if (Genos(i, j) == 0) Phase(i, j,:) = 0
		    if (Genos(i, j) == 2) Phase(i, j,:) = 1
		    CountZero = 0
		    CountOne = 0
		    do k = 1, nHapsCluster
		      if (ClusterMember(k) == 1) then
			if (library%getPhase(TempHapVector(k), j) == 0) CountZero = CountZero + 1                                                  
			if (library%getPhase(TempHapVector(k), j) == 1) CountOne = CountOne + 1                                              
		      endif
		    end do
		    if ((CountZero == 0).and.(CountOne > 0)) Phase(i, j, 1) = 1
		    if ((CountZero > 0).and.(CountOne == 0)) Phase(i, j, 1) = 0
		  end do
		endif
	      end if
	      deallocate(ClusterMember)
	      deallocate(MinClust)
	      deallocate(Medoids)
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
      if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) == 9)) then
	value = Genos(i, j) - Phase(i, j, 1)
	if ((value == 0).or.(value == 1)) then !here 7th april 2011
	  Phase(i, j, 2) = value
	else
	  CountA = CountA + 1
	endif
      endif
      if ((Phase(i, j, 2) /= 9).and.(Phase(i, j, 1) == 9)) then
	value = Genos(i, j) - Phase(i, j, 2)
	if ((value == 0).or.(value == 1)) then !here 7th april 2011
	  Phase(i, j, 1) = value
	else
	  CountB = CountB + 1
	endif
      endif
    end if
  end do
  if (CountA > ErrorCountAB) then
    Phase(i, :, :) = 9
    do j = 1, nSnpCore
      if (Genos(i, j) == 0) Phase(i, j,:) = 0
      if (Genos(i, j) == 2) Phase(i, j,:) = 1
    enddo
  endif
  if (CountB > ErrorCountAB) then
    Phase(i, : ,:) = 9
    do j = 1, nSnpCore
      if (Genos(i, j) == 0) Phase(i, j,:) = 0
      if (Genos(i, j) == 2) Phase(i, j,:) = 1
    enddo
  endif
end do

do i = 1, nAnisG
  do j = 1, nSnpCore
    if (Genos(i, j) == 1) then
      if ((Phase(i, j, 1) == 9).and.(Phase(i, j, 2) /= 9)) Phase(i, j, 1) = Genos(i, j) - Phase(i, j, 2)
      if ((Phase(i, j, 2) == 9).and.(Phase(i, j, 1) /= 9)) Phase(i, j, 2) = Genos(i, j) - Phase(i, j, 1)
    endif
    if (Genos(i, j) == 0) then
      if ((Phase(i, j, 1) == 9).and.(Phase(i, j, 2) /= 9)) Phase(i, j, 1) = 0
      if ((Phase(i, j, 2) == 9).and.(Phase(i, j, 1) /= 9)) Phase(i, j, 2) = 0
    endif
    if (Genos(i, j) == 2) then
      if ((Phase(i, j, 1) == 9).and.(Phase(i, j, 2) /= 9)) Phase(i, j, 1) = 1
      if ((Phase(i, j, 2) == 9).and.(Phase(i, j, 1) /= 9)) Phase(i, j, 2) = 1
    endif
  enddo
enddo

HapFreq = 0
FullyPhased = .false.
HapAnis = -99
!nGlobalHaps=nHaps

deallocate(CandGenos)
deallocate(CandHaps)
deallocate(WorkVec)
deallocate(CandPairs)

end subroutine ImputeFromLib

subroutine WriteHapLib(library, currentcore, phase, hapfreq)
  use Global, only: fullfileoutput, windowslinux
  implicit none
  
  type(HapLib), intent(in) :: library
  integer, intent(in) :: currentcore
  integer(kind=1), dimension(:,:,:), intent(in) :: phase
  integer, dimension(:), intent(in) :: hapfreq

  ! This should be refactored out at some point - same as currentcore
  integer :: OutputPoint
  
  integer :: i, j, k, counter, SizeCore, nHaps, nAnisG
  character(len = 300) :: filout
  
  SizeCore = size(phase,2)
  OutputPoint = currentcore
  
  nHaps = library%getSize()
  nAnisG = size(phase,1)

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
    i, HapFreq(i), " ", library%getHap(i)
    write (34) library%getHap(i)
  end do

  print*, "   ", "Final iteration found ", nHaps, "haplotypes"

  counter = 0
  do i = 1, nAnisG
    do j = 1, SizeCore
      do k = 1, 2
	if (Phase(i, j, k) == 0) counter = counter + 1
	if (Phase(i, j, k) == 1) counter = counter + 1
      end do
    end do
  end do
  print*, ""
  write (*, '(a4,a30,f5.2,a1)') "   ", "Final yield for this core was ", 100 * (float(counter)/(2 * nAnisG * SizeCore)), "%"

  write (29, '(i10,f7.2)') CurrentCore, 100 * (float(counter)/(2 * nAnisG * SizeCore))

end subroutine WriteHapLib

!NO WAY THIS SHOULD BE HERE!
subroutine WriteOutResults(phase, allHapAnis)
  use Global, only: CoreIndex, WindowsLinux, nCores, GenotypeID
  implicit none
  
  integer(kind=1), dimension(:,:,:), intent(in) :: phase
  integer, dimension(:,:,:), intent(in) :: allHapAnis

  integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp
  integer, allocatable, dimension(:) :: WorkOut
  double precision, allocatable, dimension(:) :: CoreCount

  nAnisG = size(phase,1)
  nSnp = size(phase,2)
  
  allocate(WorkOut(nCores * 2))
  allocate(CoreCount(nCores * 2))

  if (WindowsLinux == 1) then
    open (unit = 15, file = ".\PhasingResults\FinalPhase.txt", status = "unknown")
    open (unit = 25, file = ".\PhasingResults\CoreIndex.txt", status = "unknown")
    open (unit = 28, file = ".\PhasingResults\SnpPhaseRate.txt", status = "unknown")
    open (unit = 30, file = ".\PhasingResults\IndivPhaseRate.txt", status = "unknown")
    open (unit = 33, file = ".\PhasingResults\FinalHapIndCarry.txt", status = "unknown")
  else
    open (unit = 15, file = "./PhasingResults/FinalPhase.txt", status = "unknown")
    open (unit = 25, file = "./PhasingResults/CoreIndex.txt", status = "unknown")
    open (unit = 28, file = "./PhasingResults/SnpPhaseRate.txt", status = "unknown")
    open (unit = 30, file = "./PhasingResults/IndivPhaseRate.txt", status = "unknown")
    open (unit = 33, file = "./PhasingResults/FinalHapIndCarry.txt", status = "unknown")
  end if

  do i = 1, nAnisG
    write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
    Phase(i,:, 1)
    write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
    Phase(i,:, 2)
  end do

  do i = 1, nCores
    write (25, *) i, CoreIndex(i,:)
  end do


  do i = 1, nSnp
    counter = 0
    do j = 1, nAnisG
      if ((Phase(j, i, 1) == 0).or.(Phase(j, i, 1) == 1)) counter = counter + 1
      if ((Phase(j, i, 2) == 0).or.(Phase(j, i, 2) == 1)) counter = counter + 1
    end do
    write (28, '(i10,f7.2)') i, (100 * (float(counter)/(2 * nAnisG)))
  end do

  do i = 1, nAnisG
    l = 0
    do j = 1, nCores
      CounterP = 0
      CounterM = 0
      do k = CoreIndex(j, 1), CoreIndex(j, 2)
	if ((Phase(i, k, 1) == 0).or.(Phase(i, k, 1) == 1)) counterP = counterP + 1
	if ((Phase(i, k, 2) == 0).or.(Phase(i, k, 2) == 1)) counterM = counterM + 1
      end do
      l = l + 1
      CoreCount(l) = (float(counterP)/((CoreIndex(j, 2) - CoreIndex(j, 1)) + 1)) * 100
      l = l + 1
      CoreCount(l) = (float(counterM)/((CoreIndex(j, 2) - CoreIndex(j, 1)) + 1)) * 100
    end do
    write (30, '(i10,20000f7.2,20000f7.2,20000f7.2,20000f7.2)') i, CoreCount(:)
  end do

  do i = 1, nAnisG
    k = 0
    do j = 1, nCores
      k = k + 2
      WorkOut(k - 1) = AllHapAnis(i, 1, j)
      WorkOut(k) = AllHapAnis(i, 2, j)
    end do
    write (33, '(i10,20000i5,20000i5,20000i5,20000i5,20000i5)') i, WorkOut(:)
  end do

end subroutine WriteOutResults

!subroutine WriteOutPartialResults(phase, allHapAnis, coreID)
!  use Global, only: CoreIndex, WindowsLinux, nCores, GenotypeID
!  implicit none
!  
!  integer(kind=1), dimension(:,:,:), intent(in) :: phase
!  integer, dimension(:,:,:), intent(in) :: allHapAnis
!  integer, intent(in) :: coreID
!
!  integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp
!  integer, allocatable, dimension(:) :: WorkOut
!  double precision, allocatable, dimension(:) :: CoreCount
!  
!  character(:), allocatable :: coreIDtxt
!
!  nAnisG = size(phase,1)
!  nSnp = size(phase,2)
!  
!  allocate(WorkOut(nCores * 2))
!  allocate(CoreCount(nCores * 2))
!  
!  coreIDtxt = itoa(coreID)
!
!  if (WindowsLinux == 1) then
!    open (unit = 15, file = ".\PhasingResults\FinalPhase.txt", status = "unknown")
!    open (unit = 25, file = ".\PhasingResults\CoreIndex.txt", status = "unknown")
!    open (unit = 28, file = ".\PhasingResults\SnpPhaseRate.txt", status = "unknown")
!    open (unit = 30, file = ".\PhasingResults\IndivPhaseRate.txt", status = "unknown")
!    open (unit = 33, file = ".\PhasingResults\FinalHapIndCarry.txt", status = "unknown")
!  else
!    open (unit = 15, file = "./PhasingResults/FinalPhase" // coreIDtxt // ".txt", status = "unknown")
!    open (unit = 25, file = "./PhasingResults/CoreIndex" // coreIDtxt // ".txt", status = "unknown")
!    open (unit = 28, file = "./PhasingResults/SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
!    open (unit = 30, file = "./PhasingResults/IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
!    open (unit = 33, file = "./PhasingResults/FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
!  end if
!
!  do i = 1, nAnisG
!    write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
!    Phase(i,:, 1)
!    write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
!    Phase(i,:, 2)
!  end do
!
!  write (25, *) i, CoreIndex(coreID,:)
!
!
!  do i = 1, nSnp
!    counter = 0
!    do j = 1, nAnisG
!      if ((Phase(j, i, 1) == 0).or.(Phase(j, i, 1) == 1)) counter = counter + 1
!      if ((Phase(j, i, 2) == 0).or.(Phase(j, i, 2) == 1)) counter = counter + 1
!    end do
!    write (28, '(i10,f7.2)') i, (100 * (float(counter)/(2 * nAnisG)))
!  end do
!
!  do i = 1, nAnisG
!    CounterP = 0
!    CounterM = 0
!    !!!! THIS LOOP IS BROKEN FOR PARTIAL OUTPUT !!!!
!    do k = CoreIndex(j, 1), CoreIndex(j, 2)
!    if ((Phase(i, k, 1) == 0).or.(Phase(i, k, 1) == 1)) counterP = counterP + 1
!      if ((Phase(i, k, 2) == 0).or.(Phase(i, k, 2) == 1)) counterM = counterM + 1
!    end do
!    CoreCount(1) = (float(counterP)/((CoreIndex(j, 2) - CoreIndex(j, 1)) + 1)) * 100
!    CoreCount(2) = (float(counterM)/((CoreIndex(j, 2) - CoreIndex(j, 1)) + 1)) * 100
!    end do
!    write (30, '(i10,20000f7.2,20000f7.2,20000f7.2,20000f7.2)') i, CoreCount(:)
!  end do
!
!  do i = 1, nAnisG
!    k = 0
!    do j = 1, nCores
!      k = k + 2
!      WorkOut(k - 1) = AllHapAnis(i, 1, j)
!      WorkOut(k) = AllHapAnis(i, 2, j)
!    end do
!    write (33, '(i10,20000i5,20000i5,20000i5,20000i5,20000i5)') i, WorkOut(:)
!  end do
!
!end subroutine WriteOutPartialResults
!
!function itoa(i) result(res)
!  character(:),allocatable :: res
!  integer,intent(in) :: i
!  character(range(i)+2) :: tmp
!  write(tmp,'(i0)') i
!  res = trim(tmp)
!end function

end module HaplotypeLibrary