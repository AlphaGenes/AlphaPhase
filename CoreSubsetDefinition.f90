!module CoreSubsetDefinition
!  use CoreDefinition
!  
!  implicit none
!  private
!
!  type, public :: CoreSubset
!    private
!    !Almost definitely shouldn't be public but for now...
!    type(Core), allocatable :: parent
!    integer(kind = 1), allocatable, dimension(:,:,:), public :: phase
!    integer(kind = 4), allocatable, dimension(:), public :: sireGenotyped
!    integer(kind = 4), allocatable, dimension(:), public :: damGenotyped
!    integer(kind = 4), allocatable, dimension(:), public :: full2sub
!    integer(kind = 4), allocatable, dimension(:), public :: sub2full
!  contains
!    private
!    procedure, public :: create
!  end type CoreSubset
!
!contains
!
!  subroutine create(set, parent, origSireGenotyped, origDamGenotyped, members)
!    implicit none
!    
!    class(CoreSubset) :: set
!    type(Core) :: parent
!    integer(kind = 4), dimension(:), intent(in) :: origSireGenotyped, origDamGenotyped
!    logical, dimension(:), intent(in) :: members
!    
!    integer :: newNumber, origNumber, nSnp
!    integer :: i, c
!    
!    set%parent = parent
!    
!    newNumber = count(members)
!    origNumber = size(origGenos,1)
!    nSnp = size(origGenos,2)
!    
!    ! I need to find a more elegant solution to this!
!    if (allocated(set%phase)) then
!      deallocate(set%phase)
!      deallocate(set%sireGenotyped)
!      deallocate(set%damGenotyped)
!      deallocate(set%full2sub)
!      deallocate(set%sub2full)
!    end if
!    
!    allocate(set%phase(newNumber,nSnp,2))
!    allocate(set%sireGenotyped(newNumber))
!    allocate(set%damGenotyped(newNumber))
!    allocate(set%full2sub(0:origNumber))
!    allocate(set%sub2full(newNumber))
!    
!    c = 0
!    set%full2sub(0) = 0
!    
!    do i = 1, origNumber
!      if (members(i)) then
!	c = c + 1
!	set%full2sub(i) = c
!	set%sub2full(c) = i
!      else
!	set%full2sub(i) = 0
!      end if
!    end do
!    
!    do i = 1, newNumber
!      set%sireGenotyped(i) = set%full2sub(origSireGenotyped(set%sub2full(i)))
!      set%damGenotyped(i) = set%full2sub(origDamGenotyped(set%sub2full(i)))
!    end do
!  end subroutine create
!
!end module CoreSubsetDefinition