module HaplotypeLibraryDefinition
  use Constants
  use HaplotypeModule
  implicit none
  private

  type, public :: HaplotypeLibrary
    private
    type(Haplotype), dimension(:), allocatable :: newstore
    integer, dimension(:), allocatable :: hapFreq
    integer :: size
    integer :: nSnps
    integer :: storeSize, stepSize
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
    procedure, public :: getCompatHapsFreq
    final :: destroy
  end type HaplotypeLibrary
  
  interface HaplotypeLibrary
    module procedure newHaplotypeLibrary
    module procedure haplotypeLibraryFromFile
  end interface HaplotypeLibrary

contains
  function newHaplotypeLibrary(nSnps, storeSize, stepSize) result(library)
    type(HaplotypeLibrary) :: library
    integer, intent(in) :: nSnps
    integer, intent(in) :: storeSize
    integer, intent(in) :: stepSize

    library % nSnps = nSnps
    library % size = 0
    library % storeSize = storeSize
    library % stepSize = stepSize
    allocate(library % hapFreq(storeSize))
    library % hapFreq = 0    
    allocate(library%newstore(library%storeSize))
  end function newHaplotypeLibrary
  
  function haplotypeLibraryFromFile(filename, stepsize) result(library)
    type(HaplotypeLibrary) :: library
    character(*), intent(in) :: filename
    integer, intent(in) :: stepSize
    
    integer(kind=1), dimension(:), allocatable ::holdHap
    integer :: i, j
    
    open (unit=2001,file=trim(filename),status="old",form="unformatted")

    ! Read the number of Hap in the library and how long they are
    read(2001) library%storeSize,library%nSnps
    
    library%size = 0
    library%stepSize = stepSize
    allocate(library % hapFreq(library%storeSize))
    library % hapFreq = 0
    
    allocate(library%newstore(library%storeSize))

    allocate(holdHap(library%nSnps))
    do i=1,library%storeSize
      read(2001) holdHap
      j = library%matchAddHap(holdHap)
    enddo
    close (2001)
    
  end function haplotypeLibraryFromFile
  
  subroutine destroy(library)
    type(HaplotypeLibrary) :: library
    
    if (allocated(library%newstore)) then
      deallocate(library%newstore)
      deallocate(library%hapFreq)
    end if
  end subroutine destroy

  function hasHap(library, hap) result(id)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    type(Haplotype) :: h
    integer :: id

    integer :: i
    
    id = 0
    h = Haplotype(hap)
    do i = 1, library%size
      if (library%newstore(i)%compareHaplotype(h)) then
	id = i
	exit
      end if
    end do
  end function hasHap

  function addHap(library, hap) result(id)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    integer :: id

    integer :: newStoreSize
    type(Haplotype), dimension(:), pointer :: tempNewStore
    integer, dimension(:), allocatable :: tempHapFreq
    
    if (library % Size == library % storeSize) then
      newStoreSize = library % storeSize + library % stepSize
      
      allocate(tempHapFreq(library % storeSize))
      tempHapFreq = library%hapFreq
      deallocate(library%hapFreq)
      allocate(library%hapFreq(newStoreSize))
      library % hapFreq = 0
      library % hapFreq(1:library % Size) = tempHapFreq
      deallocate(tempHapFreq)
      
      allocate(tempNewStore(library % storeSize))
      tempNewStore = library%newStore
      deallocate(library%newStore)
      allocate(library%newStore(newStoreSize))
      library % newStore(1:library % Size) = tempNewStore
      deallocate(tempNewStore)
      
      library % StoreSize = newStoreSize
    end if

    library % Size = library % Size + 1
    library%newStore(library%Size) = Haplotype(hap)
    
    library%hapfreq(library%size) = 1
    id = library%size
  end function addHap
  
  function matchAddHap(library, hap) result (id)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    integer :: id
    
    id = library%hasHap(hap)
    if (id == 0) then
      id = library%addHap(hap)
    else
      library%hapfreq(id) = library%hapfreq(id) + 1
    end if    
  end function matchAddHap

  function matchWithError(library, hap, allowedError) result(matches)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    integer, intent(in) :: allowedError
    integer, dimension(:), pointer :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, e, num, invalid
    type(Haplotype) :: h

    allocate(tempMatches(library % size))
    
    invalid = 0
    do i = 1, size(hap)
      if ((hap(i) /= 0) .and. (hap(i) /= 1) .and. (hap(i) /= MissingPhaseCode)) then
	invalid = invalid + 1
      end if
    end do
    
    num = 0
    
    if (invalid <= allowedError) then
      h = Haplotype(hap)
      do i = 1, library % size
	e = invalid + library%newstore(i)%mismatchesMod(h)
	if (e <= allowedError) then
	  num = num + 1
	  tempMatches(num) = i
	end if
      end do
    end if
    
    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function matchWithError
  
  function limitedMatchWithError(library, hap, allowedError, limit) result(matches)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    integer, intent(in) :: allowedError
    integer, dimension(:), intent(in) :: limit
    integer, dimension(:), pointer :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, k, e, num, invalid    
    type(Haplotype) :: h

    allocate(tempMatches(library % size))

    num = 0
    
    invalid = 0
    do i = 1, size(hap)
      if ((hap(i) /= 0) .and. (hap(i) /= 1) .and. (hap(i) /= MissingPhaseCode)) then
	invalid = invalid + 1
      end if
    end do 
    
    if (invalid <= allowedError) then
      h = Haplotype(hap)

      do k = 1, size(limit)
	i = limit(k)
	e = invalid + library%newstore(i)%mismatchesMod(h)

	if (e <= allowedError) then
	  num = num + 1
	  tempMatches(num) = i
	end if
      end do
    end if
    
    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function limitedMatchWithError
  
  function limitedCompatPairsWithError(library, genos, ErrorAllow, limit, nAnisG) result(pairs)
    use Constants
    use GenotypeModule
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: genos
    integer, intent(in) :: ErrorAllow
    integer, dimension(:), intent(in) :: limit
    integer, intent(in) :: nAnisG
    integer, dimension(:,:), pointer :: pairs
    
    integer, dimension(:,:), pointer :: tempPairs
    integer :: i, j, p, ii, jj

    type(Genotype) :: g
    
    g = Genotype(genos)
    allocate(tempPairs(nAnisG*2,2))
    
    p = 0
    i = 1
    do while ((i <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
      j = i + 1
      do while ((j <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
	ii = limit(i)
	jj = limit(j)
	if (g%compatibleHaplotypes(library%newstore(ii), library%newstore(jj), ErrorAllow)) then
	  p = p + 1
	  tempPairs(p,1) = ii
	  tempPairs(p,2) = jj
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
    hap = library % newstore(id) % toIntegerArray()
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
    
    phase = library % newstore(id) % getPhaseMod(snp)
  end function getPhase

  function getSize(library) result(size)
    class(HaplotypeLibrary) :: library
    integer :: size

    size = library % size
  end function getSize

  function getHapRel(library) result (rel)
    class(HaplotypeLibrary) :: library
    integer, allocatable, dimension (:,:) :: rel

    integer :: i, j

    integer :: counter

    allocate(rel(library % size, library % size))

    rel = 0
    do i = 1, library % size
      do j = i + 1, library % size
	counter = library%newstore(i)%numberSame(library%newstore(j))
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
  
  function getCompatHaps(library, genos, percgenohaplodisagree) result (compatHaps)
    class(HaplotypeLibrary) :: library
    integer(kind=1), dimension(:), intent(in) :: genos
    double precision, intent(in) :: percgenohaplodisagree
    integer, dimension(:), pointer :: compatHaps
    
    compatHaps = getCompatHapsFreq(library, genos, 1, percgenohaplodisagree)    
  end function getCompatHaps
  
  function getCompatHapsFreq(library, genos, freq, percgenohaplodisagree) result (compatHaps)
    use GenotypeModule
    class(HaplotypeLibrary) :: library
    integer(kind=1), dimension(:), intent(in) :: genos
    integer, intent(in) :: freq
    double precision, intent(in) :: percgenohaplodisagree
    integer, dimension(:), pointer :: compatHaps
    
    integer, dimension(:), allocatable :: tempCompatHaps
    integer :: i, numCompatHaps, ErrorAllow
    
    type(Genotype) :: g
    
    g = Genotype(genos)
   
    ErrorAllow = int(PercGenoHaploDisagree * library%nSnps)
    allocate(tempCompatHaps(library%size))
    numCompatHaps = 0
    do i = 1, library%size
      if (library%hapFreq(i) >= freq) then
	if (g%compatibleHaplotype(library%newstore(i), ErrorAllow)) then
	  numCompatHaps = numCompatHaps + 1
	  tempCompatHaps(numCompatHaps) = i
	end if
      end if
    end do
    allocate(compatHaps(numCompatHaps))
    compatHaps = tempCompatHaps(1:numCompatHaps)
    deallocate(tempCompatHaps)
  end function getCompatHapsFreq

end module HaplotypeLibraryDefinition
