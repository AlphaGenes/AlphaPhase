module PedigreeDefinition
  use Constants
  implicit none
  private

  type, public :: Pedigree
    private    
    integer(kind = 4), dimension(:), allocatable :: sire, dam
    character(:), dimension(:), pointer :: id
    
  contains
    private
    procedure, public :: getSire
    procedure, public :: getDam
    procedure, public :: getID
    procedure, public :: getNAnis
  end type Pedigree
  
  interface Pedigree
    module procedure newPedigree
  end interface Pedigree
  
contains
  
  function newPedigree(sire, dam, id) result(p)
    use Constants
    
    integer, dimension(:), intent(in) :: sire, dam
    character(*), dimension(:), intent(in):: id(:)
    type(Pedigree) :: p

    allocate(p%sire(size(sire,1)))
    p%sire = sire
    
    allocate(p%dam(size(dam,1)))
    p%dam = dam
    
    allocate(character(len(id)) :: p%id(size(id,1)))
    p%id = id
  end function newPedigree
  
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
    character(len(p%id)) :: id
    
    id = p%id(animal)
  end function getID
  
  function getNAnis(p) result(num)
    class(Pedigree) :: p
    integer :: num
    
    num = size(p%id,1)
  end function getNAnis
  
end module PedigreeDefinition