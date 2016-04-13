module CoreSubsetDefinition
  use CoreDefinition
  
  implicit none
  private

  type, public :: CoreSubset
    private
    !Almost definitely shouldn't be public but for now...
    type(Core), pointer :: parent
    integer(kind = 4), pointer, dimension(:), public :: sireGenotyped
    integer(kind = 4), pointer, dimension(:), public :: damGenotyped
    integer(kind = 4), allocatable, dimension(:), public :: full2sub
    integer(kind = 4), allocatable, dimension(:), public :: sub2full
    integer :: nAnisG
  contains
    private
    procedure, public :: create
    procedure, public :: getCoreAndTailGenos
    procedure, public :: getCoreGenos
    procedure, public :: getSingleCoreAndTailGenos
    procedure, public :: getSingleCoreGenos
    procedure, public :: setPhase
    procedure, public :: getNAnisG
    procedure, public :: getNSnp
    procedure, public :: getNCoreSnp
    procedure, public :: getPhase
    procedure, public :: getPhaseGeno
    procedure, public :: getSire
    procedure, public :: getDam
    procedure, public :: getYield
  end type CoreSubset

contains

  subroutine create(set, parent, sireGenotyped, damGenotyped, members)
    implicit none
    
    class(CoreSubset) :: set
    type(Core), target :: parent
    integer(kind = 4), dimension(:), intent(in), target :: sireGenotyped, damGenotyped
    !logical, dimension(:), intent(in) :: members
    integer, dimension(:), intent(in) :: members
    
    !integer :: newNumber, origNumber, nSnp
    integer :: i
    
    set%parent => parent
    set%sireGenotyped => sireGenotyped
    set%damGenotyped => damGenotyped
    
    !set%nAnisG = count(members)
    set%nAnisG = size(members,1)
    
    ! I need to find a more elegant solution to this!
    if (allocated(set%full2sub)) then
      deallocate(set%full2sub)
      deallocate(set%sub2full)
    end if

    allocate(set%full2sub(0:set%parent%getNAnisG()))
    allocate(set%sub2full(set%nAnisG))
    
    set%full2sub = 0
    
    do i = 1, set%nAnisG
      set%full2sub(members(i)) = i
      set%sub2full(i) = members(i)
    end do
    
!    do i = 1, newNumber
!      set%sireGenotyped(i) = set%full2sub(origSireGenotyped(set%sub2full(i)))
!      set%damGenotyped(i) = set%full2sub(origDamGenotyped(set%sub2full(i)))
!    end do
  end subroutine create
  
  function getNAnisG(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%nAnisG
  end function getNAnisG
  
  function getNSnp(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%parent%getNSnp()
  end function getNSnp
  
  function getNCoreSnp(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%parent%getNCoreSnp()
  end function getNCoreSnp
  
  subroutine setPhase(set, animal, snp, phase, val)
    implicit none
    
    class(CoreSubset) :: set
    integer, intent(in) :: animal, snp, phase
    integer(kind=1) :: val
    
    call set%parent%setPhase(set%sub2full(animal),snp,phase,val)
  end subroutine setPhase
  
  function getPhase(set,animal,snp,phase) result(p)
    implicit none
    class(CoreSubset) :: set
    integer, intent(in) :: animal, snp, phase
    integer(kind=1) :: p
    
    p = set%parent%getPhase(set%sub2full(animal),snp,phase)
  end function getPhase
  
  function getSingleCoreAndTailGenos(set,i) result (ctGenos)
    implicit none
    
    class(CoreSubset), target :: set
    integer, intent(in) :: i
    integer(kind=1), dimension(:), pointer :: ctGenos
    
    ctGenos => set%parent%getSingleCoreAndTailGenos(set%sub2full(i))
    
    return
  end function getSingleCoreAndTailGenos
  
  function getSingleCoreGenos(set, i) result (cGenos)
    implicit none
    
    class(CoreSubset), target :: set
    integer, intent(in) :: i
    integer(kind=1), dimension(:), pointer :: cGenos
    
    cGenos => set%parent%getSingleCoreGenos(set%sub2full(i))
    
    return
  end function getSingleCoreGenos
  
  function getPhaseGeno(set,animal,snp) result (p)
    implicit none
    class(CoreSubset) :: set
    integer, intent(in) :: animal, snp
    integer(kind=1) :: p
    
    p = set%parent%getPhaseGeno(set%sub2full(animal),snp)
  end function getPhaseGeno
  
  function getCoreAndTailGenos(set) result (ctGenos)
    implicit none
    
    class(CoreSubset) :: set
    integer(kind=1), dimension(:,:), pointer :: ctGenos
    integer :: i

    allocate(ctGenos(set%nAnisG,set%parent%getNSnp()))
    
    do i = 1, set%nAnisG
      ctGenos(i,:) = set%parent%getSingleCoreAndTailGenos(set%sub2full(i))
    end do
    
    return
  end function getCoreAndTailGenos
  
  function getCoreGenos(set) result (cGenos)
    implicit none
    
    class(CoreSubset) :: set
    integer(kind=1), dimension(:,:), pointer :: cGenos
    integer :: i

    allocate(cGenos(set%nAnisG,set%parent%getNCoreSnp()))
    
    do i = 1, set%nAnisG
      cGenos(i,:) = set%parent%getSingleCoreGenos(set%sub2full(i))
    end do
    
    return
  end function getCoreGenos
  
  function getSire(set,animal) result(sire)
    class(CoreSubSet) :: set
    integer, intent(in) :: animal
    integer :: sire
    
    sire = set%full2sub(set%SireGenotyped(set%sub2full(animal)))
  end function getSire
  
  function getDam(set,animal) result(dam)
    class(CoreSubSet) :: set
    integer, intent(in) :: animal
    integer :: dam
    
    dam = set%full2sub(set%DamGenotyped(set%sub2full(animal)))
  end function getDam
  
  function getYield(set,phase) result (yield)
    implicit none
    class(CoreSubSet) :: set
    integer, intent(in) :: phase
    double precision :: yield
    integer :: counter, i
    
    counter = 0
    
    do i = 1, set%nAnisG
      counter = count(set%parent%getHaplotype(set%sub2full(i),phase) == 0) + counter
      counter = count(set%parent%getHaplotype(set%sub2full(i),phase) == 1) + counter      
    end do
    yield = (float(counter)/(set%nAnisG * set%parent%getNCoreSnp())) * 100
  end function getYield
end module CoreSubsetDefinition