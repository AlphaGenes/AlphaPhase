!!! CUURENTLY UNUSED... BUT WILL BE AT SOME POINT !!!

module PedigreeDefinition
  implicit none
  private

  type, public :: Pedigree
    private
    
    integer(kind = 4), allocatable, dimension (:), public :: sire, dam
    
  contains
    private
    procedure, public :: create
      
  end type Pedigree
  
contains
  
  subroutine create(p, sire, dam)
    implicit none

    class(Pedigree) :: p
    integer(kind = 4), dimension (:), intent(in) :: sire, dam
    
    integer :: nAnisG

    nAnisG = size(sire,1)

    if (allocated(p%sire)) then
      deallocate(p%sire)
      deallocate(p%dam)
    end if

    allocate(p%sire(nAnisG))
    allocate(p%dam(nAnisG))

    p%sire = sire
    p%dam = dam
  end subroutine create
  
end module PedigreeDefinition