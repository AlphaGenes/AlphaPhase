module ChipsDefinition
  !!!!! NEEDS BIT OPERATIONS !!!!!
  
  
  implicit none
  private

  type, public :: Chips
    private
    logical, dimension(:,:), pointer :: snps
    integer, dimension(:), pointer :: animals
  contains
    private
    procedure, public :: getChipGenotypes
    procedure, public :: getOverlap
    procedure, public :: expandHaplotype
  end type Chips
  
  interface Chips
    module procedure newChips
  end interface Chips

contains
  function newChips(snps, animals) result(c)
    logical, dimension(:,:),intent(in) :: snps
    integer, dimension(:), intent(in) :: animals
    
    type(Chips) :: c
    
    ! This should probably be pointers but it was breaking things!
    allocate(c%snps(size(snps,1),size(snps,2)))
    allocate(c%animals(size(animals)))
    
    c%snps = snps
    c%animals = animals
    
  end function newChips
  
  function getChipGenotypes(c, genos, chipID) result(chip)
    class(Chips), intent(in) :: c
    integer(kind=1), dimension(:,:), intent(in) :: genos
    integer, intent(in) :: chipID
    
    integer(kind=1), dimension(:,:), allocatable :: chip
    
    integer :: i, j, cs, ca
    
    allocate(chip(count(c%animals == chipID), count(c%snps(chipID,:))))
    
    print *, size(c%snps,2), size(genos,2)
    
    ca = 0
    do i = 1, size(c%animals)
      if (c%animals(i) == chipID) then
	ca = ca + 1
	cs = 0
	do j = 1, size(c%snps,2)
	  if (c%snps(chipID,j)) then
	    cs = cs + 1
	    chip(ca,cs) = genos(i,j)
	  end if
	end do
      end if
    end do
  end function getChipGenotypes
  
  function getOverlap(c, chipIDs) result(overlap)
    class(Chips), intent(in) :: c
    integer, dimension(:), intent(in) :: chipIDs
    
    logical, dimension(size(c%snps,2)) :: overlap
    
    integer :: i
    
    overlap = c%snps(chipIDs(1),:)
    do i = 2, size(chipIDs, 1)
      overlap = overlap .and. c%snps(chipIDs(i),:)
    end do
  end function getOverlap
  
  function expandHaplotype(c, original, chipID) result(full)
    use Constants
    
    class(Chips), intent(in) :: c
    integer(kind=1), dimension(:), intent(in) :: original
    integer :: chipID
    
    integer(kind=1), dimension(size(c%snps,2)) :: full
    
    integer :: i, cs
    
    cs = 0
    do i = 1, size(c%snps,2)
      if (c%snps(chipID,i)) then
	cs = cs + 1
	full(i) = original(cs)
      else
	full(i) = MissingPhaseCode
      end if
    end do
  end function expandHaplotype
  
  !! DAMN IT - WHAT TO DO ONE ONE HAP IS ALREADY MERGED??
  
  function matchHaplotypes(c, hap1, chipIDs1, hap2, chipIDs2, allowed) result(match)
    class(Chips) :: c
    integer(kind=1), dimension(:), intent(in) :: hap1, hap2
    integer, dimension(:), intent(in) :: chipIDs1, chipIDs2
    integer, intent(in) :: allowed
    
    logical :: match
    logical, dimension(size(chipIDs1)) :: overlap1, overlap2
    
    integer :: count, i
    
    count = 0
    
    overlap1 = getOverlap(c, chipIDs1)
    overlap2 = getOverlap(c, chipIDs2)
    
    do i = 1, size(c%snps,2)
      if (overlap1(i) .and. overlap2(i) .and. (hap1(i) /= hap2(i))) then
	! Missing?
	
	count = count + 1
	if (count > allowed) then
	  exit
	end if
      end if
    end do
    
    match = (count <= allowed)
  end function matchHaplotypes
  
  function mergeHaplotypes(c, hap1, chipID1, hap2, chipID2) result(merged)
    use Constants
    
    class(Chips) :: c
    integer(kind=1), dimension(:), intent(in) :: hap1, hap2
    integer, intent(in) :: chipID1, chipID2
    
    integer(kind=1), dimension(size(c%snps,2)) :: merged
    
    integer :: i
    
    do i = 1, size(merged)
      if (c%snps(chipID1,i)) then
	if (c%snps(chipID1,i)) then
	  if (hap1(i) == hap2(i)) then
	    merged(i) = hap1(i)
	  else
	    merged(i) = MissingPhaseCode
	  end if
	else
	  merged(i) = hap1(i)
	end if
      else
	if (c%snps(chipID2,i)) then
	  merged(i) = hap2(i)
	else
	  merged(i) = MissingPhaseCode
	end if
      end if
    end do
    
  end function mergeHaplotypes
    
end module ChipsDefinition