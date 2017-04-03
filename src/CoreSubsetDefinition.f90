module CoreSubsetDefinition
  use CoreDefinition
  use GenotypeModule
  use PedigreeModule
  
  implicit none

  type :: CoreSubset
    type(Core), pointer :: parentCore !< cor is genotype and phase for every animal, for subset of all the sSnPs.
    type(PedigreeHolder), pointer :: parentPedigree 
    integer(kind = 4), allocatable, dimension(:) :: full2sub !< mapping from fillIndex to subset fillIndex 
    integer(kind = 4), allocatable, dimension(:) :: sub2full !< mapping from subsetIndex to fullIndex
    integer :: nAnisG
  contains
    procedure :: getCoreAndTailGenosCoreSubset
    procedure :: getCoreGenosCoreSubset
    procedure :: getSingleCoreAndTailGenos
    procedure :: getSingleCoreGenos
    procedure :: getNAnisGCoreSubset
    procedure :: getNSnpCoreSubset
    procedure :: getNCoreSnpCoreSubset
    procedure :: getNCoreTailSnpCoreSubset
    procedure :: getSireCoreSubset
    procedure :: getDamCoreSubset
    procedure :: getYieldCoreSubset
    procedure :: getHaplotypeCoreSubset    
    procedure :: setSwappableCoreSubset
    
    final :: destroyCoreSubset
  end type CoreSubset
  
  interface CoreSubSet
    module procedure newCoreSubSet
  end interface CoreSubSet

contains

  function newCoreSubSet(parentCore, parentPedigree, members) result(set)
       
    type(Core), target :: parentCore
    type(PedigreeHolder), target :: parentPedigree
    integer, dimension(:), intent(in) :: members
    type(CoreSubset) :: set
    
    integer :: i
    
    set%parentCore => parentCore
    set%parentPedigree => parentPedigree
    
    set%nAnisG = size(members,1)
    
    allocate(set%full2sub(0:set%parentCore%getNAnisG()))
    allocate(set%sub2full(set%nAnisG))
    
    set%full2sub = 0
    
    do i = 1, set%nAnisG
      set%full2sub(members(i)) = i
      set%sub2full(i) = members(i)
    end do

  end function newCoreSubset
  
  subroutine destroyCoreSubset(set)
    type(CoreSubSet) :: set

    if (allocated(set%full2sub)) then
      deallocate(set%full2sub)
      deallocate(set%sub2full)
    end if
  end subroutine destroyCoreSubset
  
  function getNAnisGCoreSubset(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%nAnisG
  end function getNAnisGCoreSubset
  
  function getNSnpCoreSubset(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%parentCore%getNSnp()
  end function getNSnpCoreSubset
  
  function getNCoreSnpCoreSubset(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%parentCore%getNCoreSnp()
  end function getNCoreSnpCoreSubset
  
  function getNCoreTailSnpCoreSubset(set) result(num)
    class(CoreSubSet) :: set
    integer :: num
    
    num = set%parentCore%getNCoreTailSnp()
  end function getNCoreTailSnpCoreSubset
  
  function getSingleCoreAndTailGenos(set,i) result (ctGenos)
        
    class(CoreSubset), target :: set
    integer, intent(in) :: i
    type(Genotype), pointer :: ctGenos
    
    ctGenos => set%parentCore%coreAndTailGenos(set%sub2full(i))
    
    return
  end function getSingleCoreAndTailGenos
  
  function getSingleCoreGenos(set, i) result (cGenos)
        
    class(CoreSubset), target :: set
    integer, intent(in) :: i
    type(Genotype), pointer :: cGenos
    
    cGenos => set%parentCore%coreGenos(set%sub2full(i))
    
    return
  end function getSingleCoreGenos
  
  function getCoreAndTailGenosCoreSubset(set) result (ctGenos)
        
    class(CoreSubset) :: set
    type(Genotype), dimension(:), pointer :: ctGenos
    integer :: i

    allocate(ctGenos(set%nAnisG))
    
    do i = 1, set%nAnisG
      ctGenos(i) = set%parentCore%coreAndTailGenos(set%sub2full(i))
    end do
    
    return
  end function getCoreAndTailGenosCoreSubset
  
  function getCoreGenosCoreSubset(set) result (cGenos)
        
    class(CoreSubset) :: set
    type(Genotype), dimension(:), pointer :: cGenos
    integer :: i

    allocate(cGenos(set%nAnisG))
    
    do i = 1, set%nAnisG
      cGenos(i) = set%parentCore%coreGenos(set%sub2full(i))
    end do
    
    return
  end function getCoreGenosCoreSubset
  
  function getSireCoreSubset(set,animal) result(sire)
    class(CoreSubSet) :: set
    integer, intent(in) :: animal
    integer :: sire, s, p

    s = set%sub2full(animal)
    p = set%parentPedigree%getSireDamHDIDByIndex(set%parentPedigree%pedigree(set%parentPedigree%hdMap(s)),2)
    sire = set%full2sub(p)
  end function getSireCoreSubset
  
  function getDamCoreSubset(set,animal) result(dam)
    class(CoreSubSet) :: set
    integer, intent(in) :: animal
    integer :: dam
    
    dam = set%full2sub(set%parentPedigree%getSireDamHDIDByIndex(set%parentPedigree%pedigree(set%parentPedigree%hdMap(set%sub2full(animal))),2))
  end function getDamCoreSubset
  
  function getYieldCoreSubset(set,phase) result (yield)
    use HaplotypeModule
    class(CoreSubSet) :: set
    integer, intent(in) :: phase
    double precision :: yield
    integer :: counter, i
    type(Haplotype), pointer :: hap
    
    counter = 0
    
    do i = 1, set%nAnisG
      hap => set%parentCore%phase(set%sub2full(i),phase)
      counter = counter + hap%numberNotMissing()
    end do
    
    yield = (float(counter)/(set%nAnisG * set%parentCore%getNCoreSnp())) * 100
  end function getYieldCoreSubset
  
  subroutine setSwappableCoreSubset(set, animal, val)
    class(CoreSubSet) :: set
    integer, intent(in) :: animal
    integer(kind=1), intent(in) :: val
    
    set%parentCore%swappable(set%sub2full(animal)) = val
  end subroutine setSwappableCoreSubset
  
  function getHaplotypeCoreSubset(set, animal, phase) result(hap)
    use HaplotypeModule
    class(CoreSubSet), intent(in) :: set
    integer, intent(in) :: animal, phase
    
    type(Haplotype), pointer :: hap
    
    hap => set%parentCore%phase(set%sub2full(animal), phase)
  end function getHaplotypeCoreSubset
    
    
end module CoreSubsetDefinition