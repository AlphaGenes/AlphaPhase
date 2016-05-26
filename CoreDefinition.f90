module CoreDefinition
  implicit none
  private

  type, public :: Core
    private
    !Almost definitely shouldn't be public but for now...
    integer(kind = 1), allocatable, dimension(:,:) :: genos
    integer(kind = 1), allocatable, dimension(:,:,:), public :: phase
    logical, allocatable, dimension(:,:) :: fullyPhased
    integer, dimension(:,:), allocatable, public :: hapAnis
    
    integer :: startCoreSnp, endCoreSnp
    
    !Fudge due to oddness in old version.  Should be removed!!!!!!!
    integer :: endSurrSnp
  contains
    private
    !procedure, public :: create
    procedure, public :: getCoreAndTailGenos
    procedure, public :: getSingleCoreAndTailGenos
    procedure, public :: getSingleCoreGenos
    procedure, public :: setPhase
    procedure, public :: getNAnisG
    procedure, public :: getNSnp
    procedure, public :: getNCoreSnp
    procedure, public :: getNCoreTailSnp
    procedure, public :: getPhase
    procedure, public :: getPhaseGeno
    procedure, public :: getYield
    procedure, public :: getTotalYield
    procedure, public :: getHaplotype
    procedure, public :: setHaplotype
    procedure, public :: setHaplotypeToUnphased
    procedure, public :: resetFullyPhased
    procedure, public :: setFullyPhased
    procedure, public :: getFullyPhased
    procedure, public :: resetHapAnis
    procedure, public :: setHapAnis
    procedure, public :: getBothFullyPhased
    procedure, public :: getCoreGeno
    procedure, public :: numNotMissing
    final :: destroy
  end type Core
  
  interface Core
    module procedure newCore
  end interface Core

contains

  function newCore(genos, startCoreSnp, endCoreSnp, endSurrSnp) result(c)
    implicit none
    
    integer(kind = 1), dimension(:,:), intent(in) :: genos
    integer, intent(in) :: startCoreSnp, endCoreSnp
    type(Core) :: c
    
    integer, intent(in) :: endSurrSnp
    
    integer :: nAnisG, nSnp, nCoreSnp
    
    nAnisG = size(genos,1)
    nSnp = size(genos,2)
    nCoreSnp = endCoreSnp - startCoreSnp + 1
    
    allocate(c%genos(nAnisG,nSnp))
    allocate(c%phase(nAnisG,nCoreSnp,2))
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapAnis(nAnisG,2))
    
    c%genos = genos
    c%startCoreSnp = startCoreSnp
    c%endCoreSnp = endCoreSnp
    c%endSurrSnp = endSurrSnp
    c%fullyPhased = .false.
    c%phase = 9
    c%hapAnis = -99
  end function newCore
  
  subroutine destroy(c)
    type(Core) :: c
    
    if (allocated(c%genos)) then
      deallocate(c%genos)
      deallocate(c%phase)
      deallocate(c%fullyPhased)
      deallocate(c%hapAnis)
    end if
  end subroutine destroy
  
  function getCoreAndTailGenos(c) result (ctGenos)
    implicit none
    
    class(Core), target :: c
    integer(kind=1), dimension(:,:), pointer :: ctGenos
    
    !allocate(ctGenos(size(c%genos,1),size(c%genos,2)))
    allocate(ctGenos(size(c%genos,1),c%endSurrSnp))
    
    !ctGenos = c%genos
    ctGenos = c%genos(:,1:c%endSurrSnp)
    
    return
  end function getCoreAndTailGenos
  
  function getSingleCoreAndTailGenos(c,i) result (ctGenos)
    implicit none
    
    class(Core), target :: c
    integer, intent(in) :: i
    integer(kind=1), dimension(:), pointer :: ctGenos
    
    !allocate(ctGenos(size(c%genos,1),size(c%genos,2)))
    !allocate(ctGenos(c%endSurrSnp))
    
    !ctGenos = c%genos
    ctGenos => c%genos(i,1:c%endSurrSnp)
    
    return
  end function getSingleCoreAndTailGenos
  
  function getSingleCoreGenos(c, i) result (cGenos)
    implicit none
    
    class(Core), target :: c
    integer, intent(in) :: i
    integer(kind=1), dimension(:), pointer :: cGenos
    
    !allocate(cGenos(c%endCoreSnp - c%startCoreSnp+1))
    
    cGenos => c%genos(i,c%startCoreSnp:c%endCoreSnp)
    
    return
  end function getSingleCoreGenos
  
  subroutine setPhase(c, animal, snp, phase, val)
    implicit none
    
    class(Core) :: c
    integer, intent(in) :: animal, snp, phase
    integer(kind=1) :: val
    
    c%phase(animal,snp,phase) = val
  end subroutine setPhase
  
  function getNAnisG(c) result(num)
    implicit none
    class(Core) :: c
    integer :: num
    
    num = size(c%genos,1)
  end function getNAnisG
  
  function getNSnp(c) result(num)
    implicit none
    class(Core) :: c
    integer :: num
    
    num = size(c%genos,2)
  end function getNSnp
  
  function getNCoreSnp(c) result(num)
    implicit none
    class(Core) :: c
    integer :: num
    
    num = c%endCoreSnp-c%startCoreSnp + 1
  end function getNCoreSnp
  
  function getNCoreTailSnp(c) result(num)
    implicit none
    class(Core) :: c
    integer :: num
    
    num = c%endSurrSnp
  end function getNCoreTailSnp
  
  function getPhase(c,animal,snp,phase) result(p)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, snp, phase
    integer(kind=1) :: p
    
    p = c%phase(animal,snp,phase)
  end function getPhase
  
  function getPhaseGeno(c,animal,snp) result (p)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, snp
    integer(kind=1) :: p
    
    p = sum(c%phase(animal,snp,:))
  end function getPhaseGeno
  
  function getYield(c,phase) result (yield)
    implicit none
    class(Core) :: c
    integer, intent(in) :: phase
    integer :: counter
    double precision :: yield
    
    counter = count(c%Phase(:, :, phase) == 0)
    counter = count(c%Phase(:, :, phase) == 1) + counter
    yield = (float(counter)/(size(c%phase,1) * size(c%phase,2))) * 100
  end function getYield
  
  function getTotalYield(c) result(yield)
    implicit none
    class(Core) :: c
    integer :: counter
    double precision :: yield
    
    counter = count(c%Phase(:, :, :) == 0)
    counter = count(c%Phase(:, :, :) == 1) + counter
    yield = (float(counter)/(size(c%phase,1) * size(c%phase,2) * 2)) * 100
  end function getTotalYield
  
  function getHaplotype(c,animal, phase) result(haplotype)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, phase
    integer(kind=1), dimension(:), allocatable :: haplotype
    
    allocate(haplotype(size(c%phase,2)))
    haplotype = c%phase(animal,:,phase)
  end function getHaplotype
  
  subroutine setHaplotype(c, animal, phase, haplotype)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, phase
    integer(kind=1), dimension(:) :: haplotype
    
    c%phase(animal,:,phase) = haplotype
  end subroutine setHaplotype
  
  subroutine setHaplotypeToUnphased(c, animal, phase)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, phase
    
    c%phase(animal,:,phase) = 9
  end subroutine setHaplotypeToUnphased
  
  subroutine resetFullyPhased(c)
    implicit none
    class(Core) :: c
    
    c%fullyPhased = .false.
  end subroutine resetFullyPhased
  
  subroutine setFullyPhased(c,animal,phase)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, phase
    
    c%fullyPhased(animal,phase) = .true.
  end subroutine setFullyPhased
  
  function getFullyPhased(c,animal,phase) result(fully)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, phase
    logical :: fully
    
    fully = c%fullyPhased(animal,phase)
    !fully = all(c%phase(animal,:,phase) /= 9)
  end function getFullyPhased
  
  function getBothFullyPhased(c,animal) result(fully)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal
    logical :: fully
    
    fully = c%fullyPhased(animal,1) .and. c%fullyPhased(animal,2)
    !fully = all(c%phase(animal,:,phase) /= 9)
  end function getBothFullyPhased
    
  
  subroutine resetHapAnis(c)
    implicit none
    class(Core) :: c
    
    c%hapAnis = -99
  end subroutine resetHapAnis
  
  subroutine setHapAnis(c,animal,phase,id)
    implicit none
    class(Core) :: c
    integer, intent(in) :: animal, phase, id
    
    c%hapAnis(animal,phase) = id
  end subroutine setHapAnis
  
  function getCoreGeno(c,animal,snp) result(geno)
    class(Core) :: c
    integer, intent(in) :: animal,snp
    integer(kind = 1) :: geno
    
    geno = c%genos(animal,c%startCoreSnp + snp - 1)
  end function getCoreGeno
  
  function numNotMissing(c, animal) result(num)
    use Constants
    
    class(Core) :: c
    integer, intent(in) :: animal
    integer :: num
    
    num = count(c%genos(animal,c%startCoreSnp:c%endCoreSnp) /= MissingGenotypeCode)
  end function numNotMissing
    
end module CoreDefinition

