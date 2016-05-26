module HaplotypeLibraryDefinition
  implicit none
  public

  type, public :: HaplotypeLibrary
    private
    integer(kind = 1), dimension (:,:), allocatable :: store
    integer, dimension(:), allocatable :: hapFreq
    integer :: size
    integer :: nSnps
    integer :: storeSize, stepSize
    integer, dimension(:), allocatable :: randomOrder
  contains
    private
    procedure, public :: hasHap
    procedure :: addHap
    procedure, public :: getHap
    procedure, public :: getHaps
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
    final :: destroy
  end type HaplotypeLibrary
  
  interface HaplotypeLibrary
    module procedure newHaplotypeLibrary
  end interface HaplotypeLibrary

contains
  function newHaplotypeLibrary(nSnps, storeSize, stepSize) result(library)
    use Random
    type(HaplotypeLibrary) :: library
    integer, intent(in) :: nSnps
    integer, intent(in) :: storeSize
    integer, intent(in) :: stepSize

    integer :: nCount, secs

    library % nSnps = nSnps
    library % size = 0
    library % storeSize = storeSize
    library % stepSize = stepSize
    allocate(library % store(storeSize, nSnps))
    allocate(library % hapFreq(storeSize))
    library % store = 0
    library % hapFreq = 0
    call system_clock(nCount)
    secs = mod(nCount, int(1e6))
    allocate(library % randomOrder(library % nSnps))
    call RandomOrder(library % randomOrder, library % nSnps, 1, -abs(secs))
  end function newHaplotypeLibrary
  
  subroutine destroy(library)
    type(HaplotypeLibrary) :: library
    
    if (allocated(library%store)) then
      deallocate(library%store)
      deallocate(library%hapFreq)
      deallocate(library%randomOrder)
    end if
  end subroutine destroy

  function hasHap(library, haplotype) result(id)
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id
    integer(kind = 1), dimension(:), allocatable :: hap

    allocate(hap(library % nSnps))
    hap = library % store(id,:)
  end function getHap
  
  function getHaps(library, ids) result(haps)
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids
    integer(kind = 1), dimension(:,:), pointer :: haps
    
    integer :: i
    
    allocate(haps(size(ids,1), library%nSnps))
    
    do i = 1, size(ids,1)
      haps(i,:) = library%getHap(ids(i))
    end do
  end function getHaps

  function getPhase(library, id, snp) result(phase)
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id, snp
    integer(kind = 1) :: phase
  
    phase = library % store(id,snp)
  end function getPhase

  function getSize(library) result(size)
    class(HaplotypeLibrary) :: library
    integer :: size

    size = library % size
  end function getSize

  function getHapRel(library) result (rel)
    class(HaplotypeLibrary) :: library
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
    class(HaplotypeLibrary) :: library
    integer :: num

    num = library % nSnps
  end function getNumSnps
  
  subroutine resetHapFreq(library)
    class(HaplotypeLibrary) :: library
    
    library%HapFreq = 0
  end subroutine resetHapFreq
  
  subroutine incrementHapFreq(library, id)
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id
    
    library%HapFreq(id) = library%HapFreq(id) + 1
  end subroutine incrementHapFreq
  
  function getHapFreq(library, id) result (freq)
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id
    integer :: freq
    
    freq = library%hapFreq(id)
  end function getHapFreq
  
  function getCompatHaps(library, genos) result (compatHaps)
    use parameters, only : percgenohaplodisagree
    
    class(HaplotypeLibrary) :: library
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
  end function uniquehaps

end module HaplotypeLibraryDefinition
