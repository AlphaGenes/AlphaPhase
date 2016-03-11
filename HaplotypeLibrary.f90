module HaplotypeLibrary
  implicit none
  private

type, public :: HapLib
  private
  integer(kind=1), dimension (:,:), pointer :: store => null()
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
    
    library%nSnps = nSnps
    library%size = 0
    library%storeSize = storeSize
    library%stepSize = stepSize
    allocate(library%store(storeSize,nSnps))
    library%store = 0
    call system_clock(nCount)
    secs = mod(nCount, int(1e6))
    call RandomOrder(library%randomOrder, library%nSnps, 1, -abs(secs))
  end subroutine initalise  

  function hasHap(library, haplotype) result(id)
    class(HapLib) :: library
    integer(kind=1), dimension(:), intent(in) :: haplotype
    integer :: id
    
    integer :: i, j
    logical :: match
    
    id = 0
    
    do i = 1, library%size
      match = .true.
      do j = 1, library%nSnps
	if (library%store(i, library%RandomOrder(j)) /= haplotype(library%RandomOrder(j))) then
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
    integer(kind=1), dimension(:), intent(in) :: haplotype
    !integer :: id
    
    integer :: newStoreSize
    integer(kind=1), dimension(:,:), allocatable, target :: newStore
    
    if (library%Size == library%storeSize) then
      newStoreSize = library%storeSize + library%stepSize
      allocate(newStore(newStoreSize,library%nSnps))
      newStore = 0
      newStore(1:library%Size,:) = library%Store
      !deallocate(library%Store)
      library%Store => newStore
      library%StoreSize = newStoreSize
    end if
    
    library%Size = library%Size + 1
    library%Store(library%Size,:) = haplotype
    !id = library%Size
    
  !end function addHap
  end subroutine addHap
  
  function matchWithError(library, haplotype, allowedError) result(matches)
    class(HapLib) :: library
    integer(kind=1), dimension(:), intent(in) :: haplotype
    integer, intent(in) :: allowedError
    integer, dimension(:), allocatable :: matches
    
    integer, dimension(:), allocatable :: tempMatches
    integer :: i, j, e, num
    logical :: match
    
    allocate(tempMatches(library%size))
    
    num = 0
    
    do i = 1, library%size
      e = 0
      match = .true.
      do j = 1, library%nSnps
	if (haplotype(library%RandomOrder(j)) /= 9) then
	  if (library%store(i, library%RandomOrder(j)) /= haplotype(library%RandomOrder(j))) then
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
    integer(kind=1), dimension(:), allocatable :: hap

    hap = library%store(id,:)
  end function getHap
  
  function getPhase(library, id, snp) result(phase)
    class(HapLib) :: library
    integer, intent(in) :: id, snp
    integer(kind=1) :: phase
    
    integer(kind=1), dimension(:), allocatable :: hap
    
    hap = library%getHap(id)
    phase = hap(snp)
  end function getPhase
  
  function getSize(library) result(size)
    class(HapLib) :: library
    integer :: size
    
    size = library%size
  end function getSize
  
  function getHapRel(library) result (rel)
    class(HapLib) :: library
    integer, allocatable, dimension (:,:) :: rel
    
    integer :: i, j, k
    
    integer :: counter
    
    allocate(rel(library%size,library%size))
    
    rel = 0
    do i = 1, library%size
      do j = i + 1, library%size
	counter = 0
	!do k = StartCoreSnp, EndCoreSnp
	do k = 1, library%nSnps
	  if (library%store(i, k) == library%store(j, k)) counter = counter + 1
	end do
	rel(i, j) = counter
	rel(j, i) = counter
      end do
    end do
  end function getHapRel
  
  
  function getNumSnps(library) result (num)
    class(HapLib) :: library
    integer :: num
    
    num = library%nSnps    
  end function getNumSnps
    
end module HaplotypeLibrary