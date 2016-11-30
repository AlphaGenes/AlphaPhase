module CoreDefinition
  use Constants
  implicit none
  private

  type, public :: Core
    private
    integer(kind = 1), allocatable, dimension(:,:) :: genos
    integer(kind = 1), allocatable, dimension(:,:,:) :: phase
    logical, allocatable, dimension(:,:) :: fullyPhased
    integer, dimension(:,:), allocatable, public :: hapAnis
    
    integer(kind=1), dimension(:), allocatable :: swappable
    
    integer :: startCoreSnp, endCoreSnp
    
    !Fudge due to oddness in old version.  Should be removed!!!!!!!
    integer :: endSurrSnp
  contains
    private
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
    procedure, public :: getPercentFullyPhased
    procedure, public :: resetHapAnis
    procedure, public :: setHapAnis
    procedure, public :: getHapAnis
    procedure, public :: getBothFullyPhased
    procedure, public :: getCoreGeno
    procedure, public :: numNotMissing
    procedure, public :: hapNumMissing
    procedure, public :: getAllPhase
    
    procedure, public :: setSwappable
    procedure, public :: getSwappable
    
    final :: destroy
  end type Core
  
  interface Core
    module procedure newCore
    module procedure newPhaseCore
    module procedure newHDCore
  end interface Core

contains

  function newCore(genos, startCoreSnp, endCoreSnp, endSurrSnp) result(c)
    
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
    
    allocate(c%swappable(nAnisG))
    
    c%genos = genos
    c%startCoreSnp = startCoreSnp
    c%endCoreSnp = endCoreSnp
    c%endSurrSnp = endSurrSnp
    c%fullyPhased = .false.
    c%phase = MissingPhaseCode
    c%hapAnis = MissingHaplotypeCode
    
    c%swappable = 0
  end function newCore
  
  function newPhaseCore(phase) result(c)
    
    integer(kind = 1), dimension(:,:,:), intent(in) :: phase
    type(Core) :: c
    
    integer :: nAnisG, nSnp
    
    nAnisG = size(phase,1)
    nSnp = size(phase,2)
    
    allocate(c%phase(nAnisG,nSnp,2))
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapAnis(nAnisG,2))
    
    c%startCoreSnp = 1
    c%endCoreSnp = nSnp
    c%endSurrSnp = 0
    c%fullyPhased = .false.
    c%phase = phase
    c%hapAnis = MissingHaplotypeCode
  end function newPhaseCore
  
  !! This creates a subcore for HLI while we are using the HD hack.  Once HLI copes with missing data we should be able
  !! to remove this
  function newHDCore(orig, indivs, snps) result(c)
    type(Core), intent(in) :: orig
    integer, dimension(:), intent(in) :: indivs, snps
    type(Core) :: c
    
    integer :: nAnisG, nSnp
    integer :: i, s, oi, os

    nAnisG = size(indivs,1)
    nSnp = size(snps,1)
    
    allocate(c%genos(nAnisG,nSnp))
    allocate(c%phase(nAnisG,nSnp,2))
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapAnis(nAnisG,2))
    
    allocate(c%swappable(nAnisG))
    
    c%startCoreSnp = 1
    c%endCoreSnp = nSnp
    c%endSurrSnp = 0
    c%fullyPhased = .false.    
    c%hapAnis = MissingHaplotypeCode  
    
    do i = 1, nAnisG
      oi = indivs(i)
      do s = 1, nSnp
	os = snps(s)
	c%genos(i,s) = orig%genos(oi,os)
	c%phase(i,s,:) = orig%phase(oi,os,:)
      end do
    end do	
  end function newHDCore
  
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
        
    class(Core), target :: c
    integer(kind=1), dimension(:,:), pointer :: ctGenos
    
    !allocate(ctGenos(size(c%genos,1),size(c%genos,2)))
    allocate(ctGenos(size(c%genos,1),c%endSurrSnp))
    
    !ctGenos = c%genos
    ctGenos = c%genos(:,1:c%endSurrSnp)
    
    return
  end function getCoreAndTailGenos
  
  function getSingleCoreAndTailGenos(c,i) result (ctGenos)
        
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
        
    class(Core), target :: c
    integer, intent(in) :: i
    integer(kind=1), dimension(:), pointer :: cGenos
    
    !allocate(cGenos(c%endCoreSnp - c%startCoreSnp+1))
    
    cGenos => c%genos(i,c%startCoreSnp:c%endCoreSnp)
    
    return
  end function getSingleCoreGenos
  
  subroutine setPhase(c, animal, snp, phase, val)
        
    class(Core) :: c
    integer, intent(in) :: animal, snp, phase
    integer(kind=1) :: val
    
    c%phase(animal,snp,phase) = val
  end subroutine setPhase
  
  function getNAnisG(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = size(c%phase,1)
  end function getNAnisG
  
  function getNSnp(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = size(c%genos,2)
  end function getNSnp
  
  function getNCoreSnp(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = c%endCoreSnp-c%startCoreSnp + 1
  end function getNCoreSnp
  
  function getNCoreTailSnp(c) result(num)
    class(Core) :: c
    integer :: num
    
    num = c%endSurrSnp
  end function getNCoreTailSnp
  
  function getPhase(c,animal,snp,phase) result(p)
    class(Core) :: c
    integer, intent(in) :: animal, snp, phase
    integer(kind=1) :: p
    
    p = c%phase(animal,snp,phase)
  end function getPhase
  
  function getAllPhase(c) result(phase)
    class(Core) :: c
    integer(kind=1), dimension(:,:,:), allocatable :: phase
    
    !! This is terrible and should be chnaged.
    allocate(phase(size(c%phase,1), size(c%phase,2), size(c%phase,3)))
    
    phase = c%phase
  end function getAllPhase
  
  function getPhaseGeno(c,animal,snp) result (p)
    class(Core) :: c
    integer, intent(in) :: animal, snp
    integer(kind=1) :: p
    
    p = sum(c%phase(animal,snp,:))
  end function getPhaseGeno
  
  function getYield(c,phase) result (yield)
    class(Core) :: c
    integer, intent(in) :: phase
    integer :: counter
    double precision :: yield
    
    counter = count(c%Phase(:, :, phase) == 0)
    counter = count(c%Phase(:, :, phase) == 1) + counter
    yield = (float(counter)/(size(c%phase,1) * size(c%phase,2))) * 100
  end function getYield
  
  function getTotalYield(c) result(yield)
    class(Core) :: c
    integer :: counter
    double precision :: yield
    
    counter = count(c%Phase(:, :, :) == 0)
    counter = count(c%Phase(:, :, :) == 1) + counter
    yield = (float(counter)/(size(c%phase,1) * size(c%phase,2) * 2)) * 100
  end function getTotalYield
  
  function getHaplotype(c,animal, phase) result(haplotype)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    integer(kind=1), dimension(:), allocatable :: haplotype
    
    allocate(haplotype(size(c%phase,2)))
    haplotype = c%phase(animal,:,phase)
  end function getHaplotype
  
  subroutine setHaplotype(c, animal, phase, haplotype)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    integer(kind=1), dimension(:) :: haplotype
    
    c%phase(animal,:,phase) = haplotype
  end subroutine setHaplotype
  
  subroutine setHaplotypeToUnphased(c, animal, phase)
    class(Core) :: c
    integer, intent(in) :: animal, phase
    
    c%phase(animal,:,phase) = 9
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
  
  function hapNumMissing(c, animal, phase) result(num)
    use Constants
    
    class(Core) :: c
    integer, intent(in) :: animal, phase
    integer :: num
    
    num = count(c%phase(animal,:,phase) == MissingPhaseCode)
  end function hapNumMissing
  
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
  
  function getGeno(c,animal,snp) result(g)
    class(Core) :: c
    integer, intent(in) :: animal, snp
    integer(kind=1) :: g
    
    g = c%genos(animal,snp)
  end function getGeno
    
end module CoreDefinition

