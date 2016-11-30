module HaplotypeLibraryDefinition
!  use Constants
  use HapMod
  implicit none
  private

  type, public :: HaplotypeLibrary
    private
    integer(kind = 1), dimension (:,:), allocatable :: store
    integer(kind = 8), dimension (:,:), allocatable :: bitstore
    integer :: bitoverhang
    integer :: numsections
    type(HaplotypeType), dimension(:), pointer :: newstore
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
    allocate(library % store(storeSize, nSnps))
    allocate(library % hapFreq(storeSize))
    library % store = 0
    library % hapFreq = 0
    
    library%numsections = nSnps / 64 + 1
    library%bitoverhang = 64 - (nSnps - (library%numsections - 1) * 64)
    allocate(library%bitstore(storeSize,library%numsections))
    library%bitstore = 0
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
    allocate(library % store(library%storeSize, library%nSnps))
    allocate(library % hapFreq(library%storeSize))
    library % store = 0
    library % hapFreq = 0
    
    library%numsections = library%nSnps / 64 + 1
    library%bitoverhang = 64 - (library%nSnps - (library%numsections - 1) * 64)
    allocate(library%bitstore(library%storeSize,library%numsections))
    library%bitstore = 0
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
    
    if (allocated(library%store)) then
      deallocate(library%store)
      deallocate(library%bitstore)
!      deallocate(library%newstore)
      deallocate(library%hapFreq)
    end if
  end subroutine destroy

  function hasHap(library, hap) result(id)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    type(HaplotypeType) :: h
    integer(kind = 8), dimension(:), pointer :: bits
    integer :: id

    integer :: i, j
    logical :: match

!    bits => HaplotypeToBits(hap, library%numsections)
!    id = 0
!    do i = 1, library%size
!      match = .true.
!      do j = 1, library % numsections
!	if (library%bitstore(i,j) /= bits(j)) then
!	  match = .false.
!	  exit
!	end if
!      end do
!      if (match) then
!	id = i
!	exit
!      end if
!    end do
!    deallocate(bits)
    
    id = 0
    print *, h%compareHaplotype(h)
!    call newHaplotypeN(hap,h)
!    print *, h%sections
!    do i = 1, library%size
!      if (library%newstore(i) == h) then
!      if (h%compareHaplotype(library%newstore(i))) then
!	id = i
!	exit
!      end if
!    end do
  end function hasHap

  function addHap(library, hap) result(id)
    class(HaplotypeLibrary) :: library
    integer(kind = 1), dimension(:), intent(in) :: hap
    integer :: id

    integer :: newStoreSize
    integer(kind = 1), dimension(:,:), allocatable :: tempStore
    integer(kind = 8), dimension(:,:), allocatable :: tempBitStore
!    type(Haplotype), dimension(:), pointer :: tempNewStore
    integer, dimension(:), allocatable :: tempHapFreq
    
    integer(kind=8), dimension(:), pointer :: bits

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
      
      allocate(tempBitStore(library % storeSize, library % numsections))
      tempBitStore = library%bitStore
      deallocate(library%bitStore)
      allocate(library%bitStore(newStoreSize, library%numsections))
      library % bitStore = 0
      library % bitStore(1:library % Size, :) = tempBitStore
      deallocate(tempBitStore)
      
!      allocate(tempNewStore(library % storeSize))
!      tempNewStore = library%newStore
!      deallocate(library%newStore)
!      allocate(library%newStore(newStoreSize))
!      library % newStore(1:library % Size) = tempNewStore
!      deallocate(tempNewStore)
      
      library % StoreSize = newStoreSize
    end if

    library % Size = library % Size + 1
    library % Store(library % Size,:) = hap
    bits => HaplotypeToBits(hap, library%numsections)
    library % BitStore(library%Size,:) = bits
!    library%newStore(library%Size) = Haplotype(hap)
    deallocate(bits)
    
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
    integer :: i, j, e, num, invalid
    logical :: match
    integer(kind=8), dimension(:), pointer :: bits, present

    allocate(tempMatches(library % size))
    
    invalid = 0
    do i = 1, size(hap)
      if ((hap(i) /= 0) .and. (hap(i) /= 1) .and. (hap(i) /= MissingPhaseCode)) then
	invalid = invalid + 1
      end if
    end do 

    num = 0
 
    if (invalid <= allowedError) then
      bits => HaplotypeToBits(hap, library%numsections)
      present => HaplotypePresent(hap, library%numsections)    

      do i = 1, library % size      
	e = invalid

	do j = 1, library % numsections
	  e = e + POPCNT(IAND( &
		      IEOR(library%bitstore(i,j), bits(j)), &
		      present(j)))
	  if (e > allowedError) then
	    match = .false.
	    exit
	  endif
	end do

	if (e <= allowedError) then
	  num = num + 1
	  tempMatches(num) = i
	end if
      end do
      deallocate(bits)
      deallocate(present)
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
    integer :: i, j, k, e, num, invalid
    logical :: match
    integer(kind=8), dimension(:), pointer :: bits, present

    allocate(tempMatches(library % size))

    num = 0

!    do k = 1, size(limit)
!      i = limit(k)
!      e = 0
!      match = .true.
!      do j = 1, library % nSnps
!	if (haplotype(library % RandomOrder(j)) /= 9) then
!	  if (library % store(i, library % RandomOrder(j)) /= haplotype(library % RandomOrder(j))) then
!	    e = e + 1
!	    if (e > allowedError) then
!	      match = .false.
!	      exit
!	    endif
!	  end if
!	end if
!      end do
!      if (match) then
!	num = num + 1
!	tempMatches(num) = i
!      end if
!    end do
    
    invalid = 0
    do i = 1, size(hap)
      if ((hap(i) /= 0) .and. (hap(i) /= 1) .and. (hap(i) /= MissingPhaseCode)) then
	invalid = invalid + 1
      end if
    end do 
    
    if (invalid <= allowedError) then
      bits => HaplotypeToBits(hap, library%numsections)
      present => HaplotypePresent(hap, library%numsections)    

      do k = 1, size(limit)
	i = limit(k)
	e = invalid

	do j = 1, library % numsections
	  e = e + POPCNT(IAND( &
		      IEOR(library%bitstore(i,j), bits(j)), &
		      present(j)))
	  if (e > allowedError) then
	    match = .false.
	    exit
	  endif
	end do

	if (e <= allowedError) then
	  num = num + 1
	  tempMatches(num) = i
	end if
      end do
      
      deallocate(bits)
      deallocate(present)
    end if
    
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
    integer :: i, j, k, p, e, ii, jj
    logical :: match
    
    integer(kind=8), dimension(:), allocatable :: zeros, ones, twos, present
    integer :: cursection, curpos
    integer(kind=8) :: matchz, matcho, matcht
    
    allocate(zeros(library%numsections),ones(library%numsections),twos(library%numsections),present(library%numsections))
    zeros = 0
    ones = 0
    twos = 0
    present = 0
    
    cursection = 1
    curpos = 1
    do i = 1, size(genos)
      select case (genos(i))
      case (0)
	zeros(cursection) = ibset(zeros(cursection), curpos)
	present(cursection) = ibset(present(cursection), curpos)
      case (1)
	ones(cursection) = ibset(ones(cursection), curpos)
	present(cursection) = ibset(present(cursection), curpos)
      case (2)
	twos(cursection) = ibset(twos(cursection), curpos)
	present(cursection) = ibset(present(cursection), curpos)
      end select
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
    
    
    allocate(tempPairs(nAnisG*2,2))
    
    p = 0
    i = 1
    do while ((i <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
      j = i + 1
      do while ((j <= size(limit)) .and. ((p*p) <= (nAnisG - 1)))
	match = .true.
	e = 0
	ii = limit(i)
	jj = limit(j)
!	do k = 1, library%nSnps
!	  if (Genos(library % randomOrder(k)) /= MissingGenotypeCode)then
!	    if ((library%getPhase(limit(i), library % randomOrder(k)) &
!	    + library%getPhase(limit(j), library % randomOrder(k))) &
!	    /= Genos(library%randomOrder(k))) then
!	      e = e + 1
	do k = 1, library%numsections
	  matchz = IAND(zeros(k), NOT(IOR(library%bitstore(ii,k),library%bitstore(jj,k))))
	  matcho = IAND(ones(k), IEOR(library%bitstore(ii,k),library%bitstore(jj,k)))
	  matcht = IAND(twos(k), IAND(library%bitstore(ii,k),library%bitstore(jj,k)))
	  
	  e = e + POPCNT( IAND( &
		  IAND(NOT(matchz),NOT(matcho)), &
		  IAND(NOT(matcht),present(k)) &
		  ))
	      if (e > ErrorAllow) then
		match = .false.
		exit
	      end if
!	    end if
!	  endif
	end do
	
	if (match) then
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
    deallocate(zeros,ones,twos,present)
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
  
  function getCompatHaps(library, genos, percgenohaplodisagree) result (compatHaps)
    class(HaplotypeLibrary) :: library
    integer(kind=1), dimension(:), intent(in) :: genos
    double precision, intent(in) :: percgenohaplodisagree
    integer, dimension(:), pointer :: compatHaps
    
    compatHaps = getCompatHapsFreq(library, genos, 1, percgenohaplodisagree)    
  end function getCompatHaps
  
  function getCompatHapsFreq(library, genos, freq, percgenohaplodisagree) result (compatHaps)
    class(HaplotypeLibrary) :: library
    integer(kind=1), dimension(:), intent(in) :: genos
    integer, intent(in) :: freq
    double precision, intent(in) :: percgenohaplodisagree
    integer, dimension(:), pointer :: compatHaps
    
    integer, dimension(:), allocatable :: tempCompatHaps
    integer :: j, i, numCompatHaps, disagree, ErrorAllow
    
    integer(kind=8), dimension(:), allocatable :: zeros, twos
    integer :: cursection, curpos
   
    allocate(zeros(library%numsections),twos(library%numsections))
    zeros = 0
    twos = 0
    
    cursection = 1
    curpos = 1
    do i = 1, size(genos)
      select case (genos(i))
      case (0)
	zeros(cursection) = ibset(zeros(cursection), curpos)
      case (2)
	twos(cursection) = ibset(twos(cursection), curpos)
      end select
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
    
    ErrorAllow = int(PercGenoHaploDisagree * library%nSnps)
    allocate(tempCompatHaps(library%size))
    numCompatHaps = 0
    do i = 1, library%size
      if (library%hapFreq(i) >= freq) then
	Disagree = 0
	do j = 1, library % numsections
	  Disagree = Disagree + POPCNT(IOR( &
		    IAND(zeros(j), library%bitstore(i,j)), &
		    IAND(twos(j), NOT(library%bitstore(i,j))) &
		    ))
	end do
	if (Disagree <= ErrorAllow) then
	  numCompatHaps = numCompatHaps + 1
	  tempCompatHaps(numCompatHaps) = i
	end if
      end if
    end do
    allocate(compatHaps(numCompatHaps))
    compatHaps = tempCompatHaps(1:numCompatHaps)
    deallocate(tempCompatHaps)
  end function getCompatHapsFreq

  !!!!!! PRIVATE FUNCTIONS !!!!!!!
    
  function HaplotypeToBits(hap, numsections) result(bits)
    integer(kind=1), dimension(:), intent(in) :: hap
    integer, intent(in) :: numsections
    integer(kind=8), dimension(:), pointer :: bits
    integer :: curpos, cursection, i
    
    allocate(bits(numsections))
    cursection = 1
    curpos = 1
    bits = 0
    do i = 1, size(hap)
      select case (hap(i))
!      case (0)
!	bits(cursection) = ibclr(bits(cursection), curpos)
      case (1)
	bits(cursection) = ibset(bits(cursection), curpos)
      end select
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function  HaplotypeToBits
  
  function HaplotypePresent(hap, numsections) result(bits)
    integer(kind=1), dimension(:), intent(in) :: hap
    integer, intent(in) :: numsections
    integer(kind=8), dimension(:), pointer :: bits
    integer :: curpos, cursection, i
    
    allocate(bits(numsections))
    cursection = 1
    curpos = 1
    bits = 0
    do i = 1, size(hap)
      if ((hap(i) == 0) .or. (hap(i) == 1)) then
	bits(cursection) = ibset(bits(cursection), curpos)
      end if
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function HaplotypePresent
  
  function BitsMissingToHaplotype(bits, missing, overhang) result(hap)
    integer(kind=8), dimension(:), intent(in) :: bits, missing
    integer, intent(in) :: overhang
    integer(kind=1), dimension(:), allocatable :: hap
    
    integer :: i, j, e
    
    allocate(hap(size(bits)*64 - overhang))
    
    do i = 1, size(bits) 
      if (i < size(bits)) then
	e = 64
      else
	e = 64 - overhang
      end if
      do j = 1, e
	if (.not. btest(missing(i),j)) then
	  hap((i-1)*64+j) = 9
	else
	  if (btest(bits(i), j)) then
	    hap((i-1)*64+j) = 1
	  else
	    hap((i-1)*64+j) = 0
	  end if
	end if
      end do
    end do
    
  end function BitsMissingToHaplotype
  
  function BitsToHaplotype(bits, overhang) result(hap)
    integer(kind=8), dimension(:), intent(in) :: bits
    integer, intent(in) :: overhang
    integer(kind=1), dimension(:), allocatable :: hap
    
    integer :: i, j, e
    
    allocate(hap(size(bits)*64 - overhang))
    
    do i = 1, size(bits) 
      if (i < size(bits)) then
	e = 64
      else
	e = 64 - overhang
      end if
      do j = 1, e
	if (btest(bits(i), j)) then
	  hap((i-1)*64+j) = 1
	else
	  hap((i-1)*64+j) = 0
	end if
      end do
    end do
    
  end function BitsToHaplotype

end module HaplotypeLibraryDefinition
