module CoreDefinition
  use Constants
  use GenotypeModule
  use HaplotypeModule
  implicit none
  private

  type, public :: Core
    type(Genotype), dimension(:), pointer :: coreAndTailGenos
    type(Genotype), dimension(:), pointer :: coreGenos
    type(Haplotype), dimension(:,:), pointer :: phase
    logical, allocatable, dimension(:,:) :: fullyPhased
    integer, dimension(:,:), allocatable, public :: hapAnis
    
    integer(kind=1), dimension(:), allocatable :: swappable
    
    integer :: nCoreAndTailSnps, nCoreSnps
    
  contains
    private
    procedure, public :: getCoreAndTailGenos
    procedure, public :: getNAnisG
    procedure, public :: getNSnp
    procedure, public :: getNCoreSnp
    procedure, public :: getNCoreTailSnp
    procedure, public :: getYield
    procedure, public :: getTotalYield
!    procedure, public :: getHaplotype
    procedure, public :: setHaplotype
    procedure, public :: setHaplotypeToUnphased
    procedure, public :: resetFullyPhased
    procedure, public :: setFullyPhased
    procedure, public :: getFullyPhased
    procedure, public :: getPercentFullyPhased
    procedure, public :: resetHapAnis
    procedure, public :: setHapAnis
    procedure, public :: getHapAnis
    procedure, public :: getBothFullyPhased
    procedure, public :: getCoreGenos
    
    procedure, public :: setSwappable
    procedure, public :: getSwappable
    
    final :: destroy
  end type Core
  
  interface Core
    module procedure newCore
    module procedure newPhaseCore
!    module procedure newHDCore
  end interface Core

contains

  function newCore(genos, startCoreSnp, endCoreSnp, endSurrSnp) result(c)
    
    integer(kind = 1), dimension(:,:), intent(in) :: genos
    integer, intent(in) :: startCoreSnp, endCoreSnp
    type(Core) :: c
    
    integer, intent(in) :: endSurrSnp
    
    integer :: nAnisG, i
    
    nAnisG = size(genos,1)
    c%nCoreAndTailSnps = size(genos,2)
    c%nCoreSnps = endCoreSnp - startCoreSnp + 1
    
    allocate(c%coreAndTailGenos(nAnisG))
    allocate(c%coreGenos(nAnisG))
    allocate(c%phase(nAnisG,2))
    do i = 1, nAnisG
      c%coreGenos(i) = Genotype(genos(i,startCoreSnp:endCoreSnp))
      c%coreAndTailGenos(i) = Genotype(genos(i,1:endSurrSnp))
      c%phase(i,1) = Haplotype(c%nCoreSnps)
      c%phase(i,2) = Haplotype(c%nCoreSnps)
    end do
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapAnis(nAnisG,2))
    
    allocate(c%swappable(nAnisG))
    

    c%fullyPhased = .false.
    c%hapAnis = MissingHaplotypeCode
    
    c%swappable = 0
  end function newCore
  
  function newPhaseCore(phase) result(c)
    
    integer(kind=1), dimension(:,:,:), intent(in) :: phase
    type(Core) :: c
    
    integer :: nAnisG, nSnp, i
    
    nAnisG = size(phase,1)
    nSnp = size(phase,2)
    
    allocate(c%phase(nAnisG,2))
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapAnis(nAnisG,2))
    
    c%nCoreSnps = nSnp
    c%nCoreAndTailSnps = nSnp
    c%fullyPhased = .false.
    do i = 1, size(phase,1)
      c%phase(i,1) = Haplotype(phase(i,:,1))
      c%phase(i,2) = Haplotype(phase(i,:,2))
    end do
    c%hapAnis = MissingHaplotypeCode
  end function newPhaseCore
  
  !! This creates a subcore for HLI while we are using the HD hack.  Once HLI copes with missing data we should be able
  !! to remove this
!  function newHDCore(orig, indivs, snps) result(c)
!    type(Core), intent(in) :: orig
!    integer, dimension(:), intent(in) :: indivs, snps
!    type(Core) :: c
!    
!    integer :: nAnisG, nSnp
!    integer :: i, s, oi, os
!
!    nAnisG = size(indivs,1)
!    nSnp = size(snps,1)
!    
!    allocate(c%genos(nAnisG,nSnp))
!    allocate(c%phase(nAnisG,2))
!    allocate(c%fullyphased(nAnisG,2))
!    allocate(c%hapAnis(nAnisG,2))
!    
!    allocate(c%swappable(nAnisG))
!    
!    c%startCoreSnp = 1
!    c%endCoreSnp = nSnp
!    c%endSurrSnp = 0
!    c%fullyPhased = .false.    
!    c%hapAnis = MissingHaplotypeCode  
!    
!    do i = 1, nAnisG
!      oi = indivs(i)
!      do s = 1, nSnp
!	os = snps(s)
!	c%genos(i,s) = orig%genos(oi,os)
!	c%phase(i,s,:) = orig%phase(oi,os,:)
!      end do
!    end do	
!  end function newHDCore
  
  subroutine destroy(c)
    type(Core) :: c
    
    if (allocated(c%coreGenos)) then
      deallocate(c%fullyPhased)
      deallocate(c%hapAnis)
    end if
  end subroutine destroy
  
  function getCoreAndTailGenos(c) result (ctGenos)
        
    class(Core), target :: c
    type(Genotype), dimension(:), pointer :: ctGenos
    
    ctGenos => c%coreAndTailGenos
    
    return
  end function getCoreAndTailGenos
  
  function getNAnisG(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = size(c%phase,1)
  end function getNAnisG
  
  function getNSnp(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = c%nCoreAndTailSnps
  end function getNSnp
  
  function getNCoreSnp(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = c%nCoreSnps
  end function getNCoreSnp
  
  function getNCoreTailSnp(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = c%nCoreAndTailSnps
  end function getNCoreTailSnp
  
  function getYield(c,phase) result (yield)
    class(Core) :: c
    integer, intent(in) :: phase
    integer :: counter, i
    double precision :: yield
    
    do i = 1, size(c%phase,1)
      counter = counter + c%phase(i,phase)%numberNotMissing()
    end do
    yield = (float(counter)/(size(c%phase,1) * size(c%phase,2))) * 100
  end function getYield
  
  function getTotalYield(c) result(yield)
    class(Core) :: c
    integer :: counter, i
    double precision :: yield
    
    counter = 0
    do i = 1, size(c%phase,1)
      counter = counter + c%phase(i,1)%numberNotMissing()
      counter = counter + c%phase(i,2)%numberNotMissing()
    end do
    yield = (float(counter)/(size(c%phase,1) * c%nCoreSnps * 2)) * 100
  end function getTotalYield
  
  function getHaplotype(c,animal, phase) result(haplotype)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    type(Haplotype), pointer :: haplotype
    
    haplotype => c%phase(animal,phase)
  end function getHaplotype
  
  subroutine setHaplotype(c, animal, phase, hap)
    class(Core) :: c 
    integer, intent(in) :: animal, phase
    type(Haplotype) :: hap
    
    c%phase(animal,phase) = hap
  end subroutine setHaplotype
  
  subroutine setHaplotypeToUnphased(c, animal, phase)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    
    call c%phase(animal,phase)%setUnphased()
  end subroutine setHaplotypeToUnphased
  
  subroutine resetFullyPhased(c)
    class(Core) :: c
    
    c%fullyPhased = .false.
  end subroutine resetFullyPhased
  
  subroutine setFullyPhased(c,animal,phase)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    
    c%fullyPhased(animal,phase) = .true.
  end subroutine setFullyPhased
  
  function getFullyPhased(c,animal,phase) result(fully)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    logical :: fully
    
    fully = c%fullyPhased(animal,phase)
  end function getFullyPhased
  
  function getBothFullyPhased(c,animal) result(fully)
    class(Core) :: c
    integer, intent(in) :: animal
    logical :: fully
    
    fully = c%fullyPhased(animal,1) .and. c%fullyPhased(animal,2)
  end function getBothFullyPhased
  
  function getPercentFullyPhased(c) result (percent)
    class(Core) :: c
    double precision :: percent
    
    percent = 100.0 * float(count(c%fullyPhased)) / (size(c%fullyPhased))
  end function getPercentFullyPhased
    
  
  subroutine resetHapAnis(c)
    class(Core) :: c
    
    c%hapAnis = MissingHaplotypeCode
  end subroutine resetHapAnis
  
  subroutine setHapAnis(c,animal,phase,id)
    class(Core) :: c
    integer, intent(in) :: animal, phase, id
    
    c%hapAnis(animal,phase) = id
  end subroutine setHapAnis
  
  function getHapAnis(c,animal,phase) result(id)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    integer :: id
    
    id = c%hapAnis(animal,phase)
  end function getHapAnis
  
  subroutine setSwappable(c, animal, val)
    class(Core) :: c
    integer, intent(in) :: animal
    integer(kind=1), intent(in) :: val
  
    c%swappable(animal) = val
  end subroutine setSwappable
  
  function getSwappable(c, animal) result(val)
    class(Core) :: c
    integer, intent(in) :: animal
    integer(kind=1) :: val
    
    val = c%swappable(animal)
  end function getSwappable
  
  function getCoreGenos(c,animal) result(g)
    class(Core), intent(in) :: c
    integer, intent(in) :: animal
    
    type(Genotype), pointer :: g
    
    g => c%coreGenos(animal)
  end function getCoreGenos
    
end module CoreDefinition

