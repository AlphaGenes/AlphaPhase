!!! CUURENTLY UNUSED... BUT WILL BE AT SOME POINT !!!

module PedigreeDefinition
  use Constants
  implicit none
  private

  type, public :: Pedigree
    private    
    !integer(kind = 4), dimension (:), pointer :: sire, dam
    integer(kind = 4), dimension (:), allocatable :: sire, dam
    character(lengan), dimension(:), allocatable :: id
    
  contains
    private
    procedure, public :: create
    procedure, public :: getSire
    procedure, public :: getDam
    procedure, public :: getID
    procedure, public :: getNAnis
      
  end type Pedigree
  
contains
  
  subroutine create(p, sire, dam, id)
    use Constants
    
    class(Pedigree) :: p
    integer, dimension(:), intent(in), target :: sire, dam
    character(lengan), dimension(:), intent(in), target :: id(:)
    
    integer :: nAnisG

    allocate(p%sire(size(sire,1)))
    p%sire = sire
    
    allocate(p%dam(size(dam,1)))
    p%dam = dam
    
    allocate(p%id(size(id,1)))
    p%id = id
    
    !p%sire => sire
    !p%dam => dam
    !p%id => id
    
  end subroutine create
  
  function getSire(p,animal) result (sire)
    class(Pedigree) :: p
    integer, intent(in) :: animal
    integer :: sire
    
    sire = p%sire(animal)
  end function getSire
  
  function getDam(p,animal) result (dam)
    class(Pedigree) :: p
    integer, intent(in) :: animal
    integer :: dam
    
    dam = p%dam(animal)
  end function getDam
  
  function getID(p,animal) result (id)
    use Constants
    
    class(Pedigree) :: p
    integer, intent(in) :: animal
    character(lengan) :: id
    
    id = p%id(animal)
  end function getID
  
  function getNAnis(p) result(num)
    class(Pedigree) :: p
    integer :: num
    
    num = size(p%id,1)
  end function getNAnis
  
end module PedigreeDefinition