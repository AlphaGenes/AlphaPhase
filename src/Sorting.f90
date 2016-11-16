module Sorting
  use Constants
  
  type sortElement
    integer :: element
    integer :: index
  end type
  
  character(:), dimension(:), pointer :: sortArray
  
  contains
  
  subroutine sortWithIndex(array, indexes)
    character(*), dimension(:), intent(inout), target :: array
    integer, dimension(size(array)), intent(out) :: indexes
        
    character(len(array)), dimension(size(array)) :: tempArray
    integer :: i
    integer(8) :: as, es
    
    sortArray => array
    as = size(indexes)
    es = sizeof(indexes(1))
    
    do i = 1, size(array)
      indexes(i) = i
    end do
    
    call qsort(indexes, as, es, sortCompare)
    
    do i = 1, size(array)
      tempArray(i) = array(indexes(i))
    end do
    
    array = tempArray
  end subroutine sortWithIndex
  
  function compare(a, b) result(i)
    character(*), intent(in) :: a, b
    integer(2) :: i
   
    if (a < b) then
      i = -1
    end if
    if (a > b) then
      i = 1
    end if
    if (a == b) then
      i = 0
    end if
!    print *, a, b, i
  end function compare
  
  function sortCompare(a, b) result(i)
    integer, intent(in) :: a, b
    integer(2) :: i
   
    i = compare(sortArray(a), sortArray(b))
  end function sortCompare
end module Sorting
  
