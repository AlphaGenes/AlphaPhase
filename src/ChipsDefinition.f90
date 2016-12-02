module ChipsDefinition
  !!!!! NEEDS BIT OPERATIONS !!!!!
  
  
  implicit none
  private

  type, public :: Chips
    private
    logical, dimension(:,:), pointer :: snps
    integer, dimension(:), pointer :: animals
    integer :: numChips
  contains
    private
    procedure, public :: getChipCore
    procedure, public :: getOverlap
    procedure, public :: expandHaplotype
    procedure, public :: getAnimals
    procedure, public :: getSnps
    procedure, public :: getNumChips
    procedure, public :: mergeToAll
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
    
    c%numChips = maxval(animals)
    
  end function newChips

  function getChipCore(ch, orig, chipID) result(c)
    use CoreDefinition
    
    class(Chips), intent(in) :: ch    
    type(Core), intent(in) :: orig
    integer, intent(in) :: chipID
    
    type(Core) :: c
!    
!    if (ch%numChips == 1) then
      c = orig
!    else
!      c = Core(orig,ch%getAnimals(chipID),ch%getSNPs(chipID))
!    end if
  end function getChipCore
  
  function getAnimals(c,chipID) result(animals)
    class(Chips), intent(in) :: c
    integer, intent(in) :: chipID
    
    integer, dimension(:), allocatable :: animals
    
    integer i, ci
    
    allocate(animals(count(c%animals==chipID)))
    ci = 0
    do i = 1, size(c%animals,1)
      if (c%animals(i) == chipID) then
	ci = ci + 1
	animals(ci) = i
      end if
    end do
  end function getAnimals
  
  function getSNPs(c,chipID) result(snps)
    class(Chips), intent(in) :: c
    integer, intent(in) :: chipID
    
    integer, dimension(:), allocatable :: snps
    
    integer i, ci
    
    allocate(snps(count(c%snps(chipID,:))))
    ci = 0
    do i = 1, size(c%snps,2)
      if (c%snps(chipID,i)) then
	ci = ci + 1
	snps(ci) = i
      end if
    end do
  end function getSNPs
  
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
  
  function getNumChips(c) result(num)
    class(Chips), intent(in) :: c
    
    integer :: num
    
    num = c%numChips
  end function getNumChips  
  
  !! This is a hacky way to do this (should really be a subroutine of the first Core) but this keeps as much as possible that is
  !! to do with the temporary HD hack in one place
  subroutine mergeToAll(c, allLib, subLib, allCore, subCore)
    use HaplotypeLibraryDefinition
    use CoreDefinition
    
    class(Chips), intent(in) :: c
    class(HaplotypeLibrary), intent(inout) :: allLib
    class(HaplotypeLibrary), intent(in) :: subLib
    class(Core), intent(inout) :: allCore
    class(Core), intent(in) :: subCore
    
    
    
  end subroutine mergeToAll
  
  !! NOTHING BELOW SHOULD REALLY BE HERE BUT HERE FOR NOW TO KEEP ALL THE MULTI_HD STUFF IN ONE PLACE
  !! Will move somewhere more sensible in future
  
  function matchHaplotypes(c, hap1, hap2, allowed, minOverlap) result(match)
    use Constants
    
    class(Chips) :: c
    integer(kind=1), dimension(:), intent(in) :: hap1, hap2
    integer, intent(in) :: allowed, minOverlap
    
    logical :: match
        
    integer :: count, overlap, i
    
    count = 0
    overlap = 0
    
    do i = 1, size(c%snps,2)
      if ((hap1(i) /= MissingPhaseCode) .and. (hap2(i) /= MissingPhaseCode)) then
	overlap = overlap + 1
	if (hap1(i) /= hap2(i)) then
	  count = count + 1
	  if (count > allowed) then
	    exit
	  end if
	end if
      end if
    end do
    
    match = ((count <= allowed) .and. (overlap >= minOverlap))
  end function matchHaplotypes
  
  function mergeHaplotypes(c, hap1, hap2) result(merged)
    use Constants
    
    class(Chips) :: c
    integer(kind=1), dimension(:), intent(in) :: hap1, hap2
    
    integer(kind=1), dimension(size(c%snps,2)) :: merged
    
    integer :: i
    
    do i = 1, size(merged)
      if (hap1(i) == hap2(i)) then
	merged(i) = hap1(i)
      else
	merged(i) = MissingPhaseCode
      end if
    end do
    
  end function mergeHaplotypes
  
end module ChipsDefinition