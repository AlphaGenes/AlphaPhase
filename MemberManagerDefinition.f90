module MemberManagerDefinition
  use CoreDefinition
  use Parameters, only: itterateType, itterateNumber, numIter
  
  implicit none
  private

  type, public :: MemberManager
    private
    
    type(Core), pointer :: c
    integer :: number
    integer :: curPos
    integer, dimension(:), allocatable :: order
    integer :: numIter
    integer :: completeIter
    logical :: noneLeft
  contains
    private
    procedure, public :: hasNext
    procedure, public :: getNext
    final :: destroy      
  end type MemberManager
  
  interface MemberManager
    module procedure newMemberManager
  end interface MemberManager
  
contains

  function newMemberManager(c) result(manager)
    class(Core), intent(in) :: c
    type(MemberManager) :: manager
    
    if (itterateType .eq. "Off") then
      call createAll(manager, c)
    end if
    if (itterateType .eq. "InputOrder") then
      call createInputOrder(manager, c, itterateNumber, numIter)
    end if
    if (itterateType .eq. "RandomOrder") then
      call createRandomOrder(manager, c, itterateNumber, numIter)
    end if
  end function newMemberManager
  
  subroutine destroy(manager)
    type(MemberManager) :: manager
    
    if (allocated(manager%order)) then
      deallocate(manager%order)
    end if
  end subroutine destroy
  
  subroutine createInputOrder(manager, c, number, numIter)
    class(MemberManager) :: manager
    class(Core), intent(in), target :: c
    integer, intent(in) :: number, numIter
    
    integer :: nAnisG, i

    manager%c => c
    
    nAnisG = c%getNAnisG()

    allocate(manager%order(nAnisG))
    do i = 1, nAnisG
      manager%order(i) = i
    end do
    
    manager%noneleft = .false.
    manager%number = number
    manager%curPos = 1
    manager%numIter = numIter
    manager%completeIter = 0
  end subroutine createInputOrder
  
  subroutine createRandomOrder(manager, c, number, numIter)
    use Random
    
    class(MemberManager) :: manager
    class(Core), intent(in), target :: c
    integer, intent(in) :: number, numIter
    
    integer :: nAnisG, i, secs, nCount

    manager%c => c
    
    call system_clock(nCount)
    secs = mod(nCount, int(1e6))

    nAnisG = c%getNAnisG()
    allocate(manager%order(nAnisG))
    call RandomOrder(manager%order, nAnisG, 1, -abs(secs))
    
    manager%noneleft = .false.
    manager%number = number
    manager%curPos = 1
    manager%numIter = numIter
    manager%completeIter = 0
  end subroutine createRandomOrder
  
  subroutine createAll(manager, c)
    class(MemberManager) :: manager
    class(Core), intent(in), target :: c
    
    manager%noneleft = .false.
    call createInputOrder(manager, c, c%getNAnisG(),1)
  end subroutine createAll
  
  function hasNext(manager) result (has)
    class(MemberManager) :: manager
    logical :: has
    
    has = ((manager%completeIter < manager%numIter) .and. (.not. manager%noneLeft))
  end function hasNext
  
  function getNext(manager) result(members)
    use Random
    class(MemberManager) :: manager
    integer, dimension(:), allocatable :: members
    
    !integer, dimension(:), allocatable :: tempMembers
    integer :: maxNumber, unPhased
    integer :: c, i
    
    integer :: nCount, secs, nAnisG
    
    unPhased = 0
    do i = 1, manager%c%getNAnisG()
      if (.not. manager%c%getBothFullyPhased(i)) then
	unPhased = unPhased + 1
      end if
    end do
    if (unPhased < manager%number) then
      manager%noneLeft = .true.
      maxNumber = unPhased
    else
      maxNumber = manager%number
    end if
    allocate(members(maxNumber))
    c = 0
    !allocate(tempMembers(manager%number))
    !do while ((manager%curPos < size(manager%order,1)) .and. (c < manager%number))
    !do while (c < manager%number)
    do while (c < maxNumber)
      if (.not. manager%c%getBothFullyPhased(manager%curPos)) then
	c  = c + 1
	!tempMembers(c) = manager%curPos
	members(c) = manager%curPos
      end if
      if (manager%curPos == size(manager%order,1)) then
	manager%curPos = 1
	manager%completeIter = manager%completeIter + 1
	if (itterateType .eq. "RandomOrder") then
	  call system_clock(nCount)
	  secs = mod(nCount, int(1e6))
	  call RandomOrder(manager%order, size(manager%order,1), 1, -abs(secs))
	end if
      else
	manager%curPos = manager%curPos + 1
      end if
    end do
    
    !allocate(members(c))
    !members(:) = tempMembers(1:c)
  end function getNext
  
end module MemberManagerDefinition