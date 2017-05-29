module HaplotypeLibraryModule
  use ConstantModule
  use HaplotypeModule
  implicit none

  type :: HaplotypeLibrary
    type(Haplotype), dimension(:), allocatable :: newstore
    integer, dimension(:), allocatable :: hapFreq
    integer :: size
    integer :: nSnps
    integer :: storeSize, stepSize
  contains
    procedure :: hasHap
    procedure :: addHap
    procedure :: getHap
    procedure :: getHaps
    procedure :: matchWithError
    procedure :: matchWithErrorAndMinOverlap
    procedure :: limitedMatchWithError
    procedure :: limitedMatchWithErrorAndMinOverlap
    procedure :: limitedCompatPairsWithError
    procedure :: getSize
    procedure :: setSize
    procedure :: getHapRel
    procedure :: getNumSnps
    procedure :: setNumSnps
    procedure :: resetHapFreq
    procedure :: incrementHapFreq
    procedure :: getHapFreq
    procedure :: setHapFreq
    procedure :: getCompatHaps
    procedure :: getCompatHapsFreq
    procedure :: allZero
    procedure :: allOne
    procedure :: allZeroOrMissing
    procedure :: allOneOrMissing
    procedure :: allMissing
    procedure :: oneZeroNoOnes
    procedure :: oneOneNoZeros
    procedure :: numberFullyPhased
    procedure :: numberPercentPhased
    procedure :: rationalise
    procedure :: removeHap
    procedure :: updateHap
    final :: destroyHaplotypeLibrary
  end type HaplotypeLibrary

  interface HaplotypeLibrary
    module procedure newHaplotypeLibrary
    module procedure haplotypeLibraryFromFile
  end interface HaplotypeLibrary

contains
  subroutine destroyHaplotypeLibrary(library)
    type(HaplotypeLibrary) :: library
    
    if (allocated(library%newstore)) then
      deallocate(library%newstore)
      deallocate(library%hapfreq)
    end if
  end subroutine destroyHaplotypeLibrary

  function newHaplotypeLibrary(nSnps, storeSize, stepSize) result(library)
    type(HaplotypeLibrary) :: library
    integer, intent(in) :: nSnps
    integer, intent(in) :: storeSize
    integer, intent(in) :: stepSize

    integer :: i
    
    library % nSnps = nSnps
    library % size = 0
    library % storeSize = storeSize
    library % stepSize = stepSize
    allocate(library % hapFreq(storeSize))
    library % hapFreq = 0
    allocate(library%newstore(library%storeSize))
    
    do i = 1, storeSize
      library%newstore(i) = newHaplotypeMissing(nSnps)
    end do
  end function newHaplotypeLibrary

  function haplotypeLibraryFromFile(filename, stepsize, text) result(library)
     ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    
    type(HaplotypeLibrary) :: library
    character(*), intent(in) :: filename
    integer, intent(in) :: stepSize
    logical, optional, intent(in) :: text

    integer(kind=1), dimension(:), allocatable ::holdHap
    integer :: i, j, f
    character(len=10000) :: line
    logical :: binary

    binary = .true.
    if (present(text)) then
      if (text) then
	binary = .false.
      end if
    end if

    if (binary) then
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
	j = library%addHap(Haplotype(holdHap))
      enddo
      close (2001)
    else
      library%storeSize = 0
      open(unit=2001,file=trim(filename),status="old")
      do
	read(2001,'(14X,A10000)', iostat=i) line
	if (i == 0) then
	  library%storeSize = library%storeSize + 1
	else
	  exit
	endif
      enddo

      library%nSnps = len(trim(line))
      library%size = 0
      library%stepSize = stepSize
      allocate(library % hapFreq(library%storeSize))
      library % hapFreq = 0
      allocate(library%newstore(library%storeSize))
      allocate(holdHap(library%nSnps))

      rewind(2001)
      do i = 1, library%storeSize
	read(2001,'(6X,I6,2X,'//itoa(library%nSnps)//'I1)') f, holdHap
	j = library%addHap(Haplotype(holdHap))
	library%hapFreq(i) = f
      end do
      close(2001)

    end if

  end function haplotypeLibraryFromFile

  function hasHap(library, hap) result(id)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    
    class(HaplotypeLibrary) :: library
    type(Haplotype), intent(in) :: hap
    integer :: id

    integer :: i

    id = 0
    do i = 1, library%size
      if (library%newstore(i)%compareHaplotype(hap)) then
	id = i
	exit
      end if
    end do
  end function hasHap

  function addHap(library, hap) result(id)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    type(Haplotype), intent(in) :: hap
    integer :: id

    integer :: newStoreSize
    type(Haplotype), dimension(:), allocatable :: tempNewStore
    integer, dimension(:), allocatable :: tempHapFreq
    
    if (library % Size == library % storeSize) then
      newStoreSize = library % storeSize + library % stepSize

      allocate(tempHapFreq(library % storeSize))
      tempHapFreq = library%hapFreq
      if (allocated(library%hapFreq)) then
        deallocate(library%hapFreq)
      endif
      allocate(library%hapFreq(newStoreSize))
      library % hapFreq = 0
      library % hapFreq(1:library % Size) = tempHapFreq
      allocate(tempNewStore(library % storeSize))
      tempNewStore = library%newStore
      if (allocated(library%newStore)) then
        deallocate(library%newStore)
      endif
      allocate(library%newStore(newStoreSize))
      library % newStore(1:library % Size) = tempNewStore
      library % StoreSize = newStoreSize
      deallocate(tempNewStore)
      deallocate(tempHapFreq)
    end if

    library % Size = library % Size + 1
    library%newStore(library%Size) = hap

    library%hapfreq(library%size) = 1
    id = library%size
  end function addHap
  
  subroutine updateHap(library, i, hap)
      use haplotypemodule
      class(HaplotypeLibrary) :: library
    integer, intent(in) :: i
    type(Haplotype) :: hap
    
    library%newstore(i) = hap    
  end subroutine updateHap
  
  function matchWithError(library, hap, allowedError) result(matches)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    type(Haplotype), intent(in) :: hap
    integer, intent(in) :: allowedError
    integer, dimension(:), pointer :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, e, num, invalid

    allocate(tempMatches(library % size))

    num = 0

    invalid = hap%numberError()

    if (invalid <= allowedError) then

      do i = 1, library%size
	e = invalid + library%newstore(i)%mismatches(hap)

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

  function matchWithErrorAndMinOverlap(library, hap, allowedError, minOverlap) result(matches)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    type(Haplotype), intent(in) :: hap
    integer, intent(in) :: allowedError
    integer, intent(in) :: minOverlap
    integer, dimension(:), allocatable :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, e, num, invalid

    allocate(tempMatches(library % size))

    num = 0

    invalid = hap%numberError()

    if (invalid <= allowedError) then
      if (allowedError == 0) then
	do i = 1, library%size
	  if (library%newstore(i)%noMismatches(hap)) then
	    if (library%newstore(i)%overlap(hap) >= minOverlap) then
	      num = num + 1
	      tempMatches(num) = i
	    end if
	  end if
	end do
      else
	do i = 1, library%size
	  e = invalid + library%newstore(i)%mismatches(hap)

	  if (e <= allowedError) then
	    if (library%newstore(i)%overlap(hap) >= minOverlap) then
	      num = num + 1
	      tempMatches(num) = i
	    end if
	  end if
	end do
      end if
    end if

    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function matchWithErrorAndMinOverlap

  function limitedMatchWithError(library, hap, allowedError, limit) result(matches)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    type(Haplotype), intent(in) :: hap
    integer, intent(in) :: allowedError
    integer, dimension(:), intent(in) :: limit
    integer, dimension(:), pointer :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, k, e, num, invalid

    allocate(tempMatches(library % size))

    num = 0

    invalid = hap%numberError()

    if (invalid <= allowedError) then

      do k = 1, size(limit)
	i = limit(k)
	e = invalid + library%newstore(i)%mismatches(hap)

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

  function limitedMatchWithErrorAndMinOverlap(library, hap, allowedError, minOverlap, limit) result(matches)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    type(Haplotype), intent(in) :: hap
    integer, intent(in) :: allowedError
    integer, dimension(:), intent(in) :: limit
    integer, intent(in) :: minOverlap
    integer, dimension(:), allocatable :: matches

    integer, dimension(:), allocatable :: tempMatches
    integer :: i, k, e, num, invalid

    allocate(tempMatches(library % size))

    num = 0

    invalid = hap%numberError()

    if (invalid <= allowedError) then

      do k = 1, size(limit)
	i = limit(k)
	e = invalid + library%newstore(i)%mismatches(hap)

	if (e <= allowedError) then
	  if (library%newstore(i)%overlap(hap) >= minOverlap) then
	    num = num + 1
	    tempMatches(num) = i
	  end if
	end if
      end do
    end if

    allocate(matches(num))
    matches(:) = tempMatches(1:num)
    deallocate(tempMatches)
  end function limitedMatchWithErrorAndMinOverlap

  function limitedCompatPairsWithError(library, geno, ErrorAllow, limit, nAnisG) result(pairs)
    ! No sensible way to make this funciton work with multi-HD / minOverlap stuff so leaving unchanged
    !! Should probably do something like this at end with only fully phased haps.
    
    use GenotypeModule
    class(HaplotypeLibrary) :: library
    type(Genotype), intent(in) :: geno
    integer, intent(in) :: ErrorAllow
    integer, dimension(:), intent(in) :: limit
    integer, intent(in) :: nAnisG
    integer, dimension(:,:), pointer :: pairs

    integer, dimension(:,:), pointer :: tempPairs
    integer :: i, j, p, ii, jj

    allocate(tempPairs(nAnisG*2,2))

    p = 0
    i = 1
    do while ((i <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
      j = i + 1
      do while ((j <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
	ii = limit(i)
	jj = limit(j)
	if (geno%compatibleHaplotypes(library%newstore(ii), library%newstore(jj), ErrorAllow)) then
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
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id
    type(Haplotype) :: hap


    hap = library % newstore(id)
  end function getHap

  function getHaps(library, ids) result(haps)
    ! Following use statements needed here due to compiler issues (Daniel / 16.0.0)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids
    type(Haplotype), dimension(:), pointer :: haps

    integer :: i

    allocate(haps(size(ids,1)))

    do i = 1, size(ids,1)
      haps(i) = library%getHap(ids(i))
    end do
  end function getHaps

  function getSize(library) result(size)
    class(HaplotypeLibrary) :: library
    integer :: size

    size = library % size
  end function getSize

  subroutine setSize(library, size)

    class(HaplotypeLibrary) :: library
    integer,intent(in) :: size

      library%size = size
  end subroutine setSize

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

    subroutine setNumSnps(library, num)
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: num

    library % nSnps = num
  end subroutine setNumSnps

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


  subroutine setHapFreq(library, id, freq)
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id,freq

    library%hapFreq(id) = freq
  end subroutine setHapFreq

  function getCompatHaps(library, g, errorAllow) result (compatHaps)
    use GenotypeModule
    class(HaplotypeLibrary) :: library
    type(Genotype), intent(in), pointer :: g
    integer, intent(in) :: errorAllow
    integer, dimension(:), pointer :: compatHaps
    
    compatHaps = getCompatHapsFreq(library, g, 1, errorallow)    
  end function getCompatHaps
  
  function getCompatHapsFreq(library, g, freq, errorallow) result (compatHaps)
    use GenotypeModule
    class(HaplotypeLibrary) :: library
    type(Genotype), intent(in), pointer :: g
    integer, intent(in) :: freq, errorallow
    integer, dimension(:), allocatable :: compatHaps

    integer, dimension(:), allocatable :: tempCompatHaps
    integer :: i, numCompatHaps 
    
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

  function allOne(library, ids) result(all)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all

    integer :: i, j, sections
    type(Haplotype) :: hap

    sections = library%nSnps/64+1

    all = 0
    all = NOT(all)

    hap = Haplotype(library%nSnps)

    do i = 1, size(ids)
      hap = library%newstore(ids(i))
      do j = 1, sections
	all(j) = IAND(all(j), IAND(hap%phase(j), NOT(hap%missing(j))))
      end do
    end do

  end function allOne

  function allZero(library, ids) result(all)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all

    integer :: i, j, sections, overhang
    type(Haplotype) :: hap

    sections = library%nSnps/64+1

    all = 0
    all = NOT(all)

    hap = Haplotype(library%nSnps)

    do i = 1, size(ids)
      hap = library%newstore(ids(i))
      do j = 1, sections
	all(j) = IAND(all(j), IAND(NOT(hap%phase(j)), NOT(hap%missing(j))))
      end do
    end do

    overhang = 64 - (library%nSnps - (sections - 1) * 64)
    do i = 64 - overhang + 1, 64
        all(sections) = ibclr(all(sections), i - 1)
    end do

  end function allZero

  function allOneOrMissing(library, ids) result(all)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all

    integer :: i, j, sections
    type(Haplotype) :: hap

    sections = library%nSnps/64+1

    all = 0
    all = NOT(all)

    hap = Haplotype(library%nSnps)

    do i = 1, size(ids)
      hap = library%newstore(ids(i))
      do j = 1, sections
	all(j) = IAND(all(j), IOR(hap%phase(j), hap%missing(j)))
      end do
    end do

  end function allOneOrMissing

  function allZeroOrMissing(library, ids) result(all)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all

    integer :: i, j, sections, overhang
    type(Haplotype) :: hap

    sections = library%nSnps/64+1

    all = 0
    all = NOT(all)

    hap = Haplotype(library%nSnps)

    do i = 1, size(ids)
      hap = library%newstore(ids(i))
      do j = 1, sections
	all(j) = IAND(all(j), IOR(NOT(hap%phase(j)), hap%missing(j)))
      end do
    end do

    overhang = 64 - (library%nSnps - (sections - 1) * 64)
    do i = 64 - overhang + 1, 64
        all(sections) = ibclr(all(sections), i - 1)
    end do

  end function allZeroOrMissing

  function allMissing(library, ids) result(all)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all

    integer :: i, j, sections
    type(Haplotype) :: hap

    sections = library%nSnps/64+1

    all = 0
    all = NOT(all)

    hap = Haplotype(library%nSnps)

    do i = 1, size(ids)
      hap = library%newstore(ids(i))
      do j = 1, sections
	all(j) = IAND(all(j), hap%missing(j))
      end do
    end do

  end function allMissing

  function oneZeroNoOnes(library, ids) result(all)
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all, azm, am

    integer :: i, sections

    sections = library%nSnps/64+1

    azm = library%allZeroOrMissing(ids)
    am = library%allMissing(ids)

    do i = 1, sections
      all(i) = IAND(azm(i), NOT(am(i)))
    end do

  end function oneZeroNoOnes

  function oneOneNoZeros(library, ids) result(all)
    class(HaplotypeLibrary) :: library
    integer, dimension(:), intent(in) :: ids

    integer(kind=8), dimension(library%nSnps/64 + 1) :: all, aom, am

    integer :: i, sections

    sections = library%nSnps/64+1

    aom = library%allOneOrMissing(ids)
    am = library%allMissing(ids)

    do i = 1, sections
      all(i) = IAND(aom(i), NOT(am(i)))
    end do

  end function oneOneNoZeros

  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
  
  function rationalise(library, percThreshold, c) result(newlib)
    use CoreModule
    
    class(HaplotypeLibrary), intent(in) :: library
    double precision, intent(in) :: percThreshold
    type(Core), optional :: c
    type(HaplotypeLibrary) :: newlib
    
    integer :: threshold, i
    integer, dimension(library%size) :: newIDs
    
    newlib = HaplotypeLibrary(library%nSnps,library%size,library%stepsize)
    
    threshold = int(percThreshold * library%nSnps)
    newIDs = MissingHaplotypeCode
    
    do i = 1, library%size
      if (library%newstore(i)%numberNotMissing() >= threshold) then
	newlib%size = newlib%size + 1
	newlib%newstore(newlib%size) = library%newstore(i)
	newlib%hapfreq(newlib%size) = library%hapfreq(i)
	newIDs(i) = newlib%size
      end if
    end do
    
    if (present(c)) then
      do i = 1, size(c%hapAnis,1)
	if (c%hapAnis(i,1) /= MissingHaplotypeCode) then
	  c%hapAnis(i,1) = newIDs(c%hapAnis(i,1))
	end if
	if (c%hapAnis(i,2) /= MissingHaplotypeCode) then
	  c%hapAnis(i,2) = newIDs(c%hapAnis(i,2))
	end if
      end do
    end if
  end function rationalise
  
  subroutine removeHap(library, id)
    class(HaplotypeLibrary) :: library
    integer, intent(in) :: id
    
    library%newstore(id:library%size-1) = library%newstore(id+1:library%size)
    library%size = library%size - 1
  end subroutine removeHap
  
  function numberFullyPhased(library) result(num)
    class(HaplotypeLibrary), intent(in) :: library
    integer :: num
    
    integer :: i
    
    num = 0
    do i = 1, library%size
      if (library%newstore(i)%fullyphased()) then
	num = num + 1
      end if
    end do
  end function numberFullyPhased
  
  function numberPercentPhased(library, percThreshold) result (num)
    class(HaplotypeLibrary), intent(in) :: library
    double precision, intent(in) :: percThreshold
    integer :: num
    
    integer :: threshold
    
    integer :: i
    
    threshold = int(percThreshold * library%nSnps)
    
    num = 0
    do i = 1, library%size
      if (library%newstore(i)%numberNotMissing() >= threshold) then
	num = num + 1
      end if
    end do
  end function numberPercentPhased

end module HaplotypeLibraryModule
