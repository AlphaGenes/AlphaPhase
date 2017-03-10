module CoreDefinition
  use ConstantModule
  use GenotypeModule
  use HaplotypeModule
  implicit none

  type :: Core
    type(Genotype), dimension(:), pointer :: coreAndTailGenos
    type(Genotype), dimension(:), pointer :: coreGenos
    type(Haplotype), dimension(:,:), pointer :: phase
    logical, allocatable, dimension(:,:) :: fullyPhased
    integer, dimension(:,:), allocatable :: hapAnis
    
    integer(kind=1), dimension(:), allocatable :: swappable
    
    integer :: nCoreAndTailSnps, nCoreSnps
    
  contains
    procedure :: getCoreAndTailGenos
    procedure :: getNAnisG
    procedure :: getNSnp
    procedure :: getNCoreSnp
    procedure :: getNCoreTailSnp
    procedure :: getYield
    procedure :: getTotalYield
    procedure :: setHaplotype
    procedure :: setHaplotypeToUnphased
    procedure :: resetFullyPhased
    procedure :: setFullyPhased
    procedure :: getFullyPhased
    procedure :: getPercentFullyPhased
    procedure :: resetHapAnis
    procedure :: setHapAnis
    procedure :: getHapAnis
    procedure :: getBothFullyPhased
    procedure :: getCoreGenos
    
    procedure :: setSwappable
    procedure :: getSwappable
    procedure :: flipHaplotypes
    procedure :: writeCore
    final :: destroyCore
  end type Core
  
  interface Core
    module procedure newCore
    module procedure newPhaseCore
  end interface Core

contains

  function newCore(p, startTailSnp, startCoreSnp, endCoreSnp, endTailSnp) result(c)
    use PedigreeModule

    type(PedigreeHolder), intent(in) :: p
    integer, intent(in) :: startCoreSnp, endCoreSnp, startTailSnp, endTailSnp
    type(Core) :: c
    type(Genotype) :: tempFullGeno
    
    integer :: nAnisG, i
    
    
    nAnisG = p%nHd
    c%nCoreAndTailSnps = endTailSnp - startTailSnp + 1
    c%nCoreSnps = endCoreSnp - startCoreSnp + 1
    
    allocate(c%coreAndTailGenos(nAnisG))
    allocate(c%coreGenos(nAnisG))
    allocate(c%phase(nAnisG,2))
    
    do i = 1, nAnisG
      tempFullGeno = p%pedigree(p%hdMap(i))%individualGenotype
      c%coreGenos(i) = tempFullGeno%subset(startCoreSnp, endCoreSnp)
      c%coreAndTailGenos(i) = tempFullGeno%subset(startTailSnp, endTailSnp)
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
  
  function newPhaseCore(phase, startSnp, endSnp) result(c)
    
    type(Haplotype), dimension(:,:), intent(in) :: phase
    integer, intent(in) :: startSnp, endSnp
    type(Core) :: c
    
    integer :: nAnisG, i
    type(Haplotype) :: tempFullHaplotype
    
    nAnisG = size(phase,1)
    
    allocate(c%phase(nAnisG,2))
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapAnis(nAnisG,2))
    
    c%nCoreSnps = endSnp - startSnp + 1
    c%nCoreAndTailSnps = c%nCoreSnps
    c%fullyPhased = .false.
    do i = 1, size(phase,1)
      tempFullHaplotype = Phase(i,1)
      c%phase(i,1) = tempFullHaplotype%subset(startSnp,endSnp)
      tempFullHaplotype = Phase(i,2)
      c%phase(i,2) = tempFullHaplotype%subset(startSnp,endSnp)
    end do
    c%hapAnis = MissingHaplotypeCode
  end function newPhaseCore
  
  subroutine destroyCore(c)
    type(Core) :: c
    
    if (allocated(c%fullyPhased)) then
      deallocate(c%fullyPhased)
      deallocate(c%hapAnis)
    end if
  end subroutine destroyCore
  
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
  
  subroutine flipHaplotypes(c, TruePhase)
    use HaplotypeModule

    class(Core) :: c
    type(Haplotype), dimension(:,:), pointer, intent(in) :: TruePhase

    integer :: i, j, CountAgreeStay1, CountAgreeStay2, CountAgreeSwitch1, CountAgreeSwitch2, truth
    type(Haplotype), pointer :: W1, W2
    type(Haplotype), pointer :: hap1, hap2, truehap1, truehap2
    integer :: HA1, HA2
    logical :: FP1, FP2

    integer :: nAnisG, nSnp

    nAnisG = c%getNAnisG()
    nSnp = c%getNCoreSnp()

    do i = 1, nAnisG
      CountAgreeStay1 = 0
      CountAgreeStay2 = 0
      CountAgreeSwitch1 = 0
      CountAgreeSwitch2 = 0
      truth = 0
      hap1 => c%phase(i,1)
      hap2 => c%phase(i,2)
      truehap1 => TruePhase(i,1)
      truehap2 => TruePhase(i,2)
      !! This could be done with bit ops but is not a priority given it's only used when comparing to known phase
      do j = 1, nSnp
	if (truehap1.getPhaseMod(j) == hap1%getPhaseMod(j)) CountAgreeStay1 = CountAgreeStay1 + 1
	if (truehap2.getPhaseMod(j) == hap1%getPhaseMod(j)) CountAgreeSwitch1 = CountAgreeSwitch1 + 1
	if (truehap1.getPhaseMod(j) == hap2%getPhaseMod(j)) CountAgreeSwitch2 = CountAgreeSwitch2 + 1
	if (truehap2.getPhaseMod(j) == hap2%getPhaseMod(j)) CountAgreeStay2 = CountAgreeStay2 + 1
      end do
      if ((CountAgreeSwitch2 > CountAgreeStay2).and.(CountAgreeStay1 <= CountAgreeSwitch1)) truth = 1
      if ((CountAgreeSwitch1 > CountAgreeStay1).and.(CountAgreeStay2 <= CountAgreeSwitch2)) truth = 1
      if (truth == 1) then
	W1 => c%phase(i,1)
	W2 => c%phase(i,2)
	call c%setHaplotype(i,1,W2)
	call c%setHaplotype(i,2,W1)
	
	HA1 = c%hapAnis(i,1)
	HA2 = c%hapAnis(i,2)
	c%hapAnis(i,1) = HA2
	c%hapAnis(i,2) = HA1
	
	FP1 = c%fullyPhased(i,1)
	FP2 = c%fullyPhased(i,2)
	c%fullyPhased(i,1) = FP2
	c%fullyPhased(i,2) = FP1
      end if
    end do

  end subroutine flipHaplotypes

! 

! TODO check with daniel for better way to do this
  subroutine writeCore(coreIn,fileIn)

    class(Core) :: coreIn
    character(len=*), intent(in) :: fileIn
    integer :: i,h, fileUnit


    open(newUnit=fileUnit, file=fileIn, status="unknown")
    write(fileUnit,*) size(coreIn%coreAndTailGenos) !< write number first
    do i=1, size(coreIn%coreAndTailGenos)
      write(fileUnit,*) coreIn%coreAndTailGenos(i)
    enddo
    write(fileUnit,*) size(coreIn%coreGenos) !< write number first
    do i=1, size(coreIn%coreGenos)
      write(fileUnit,*) coreIn%coreGenos(i)
    enddo                 
    write(fileUnit,*) size(coreIn%phase,1),size(coreIn%phase,1)  !< write number first
    do i=1, size(coreIn%phase,1)
      do h=1, size(coreIn%phase,2)
        write(fileUnit,*) coreIn%phase(i,h)
      enddo
    enddo

    do i=1, size(coreIn%fullyPhased,1)
      do h=1, size(coreIn%fullyPhased,2)
        write(fileUnit,*) coreIn%fullyPhased(i,h)
      enddo
    enddo

    do i=1, size(coreIn%hapAnis,1)
      do h=1, size(coreIn%hapAnis,2)
        write(fileUnit,*) coreIn%hapAnis(i,h)
      enddo
    enddo

  end subroutine writeCore
    
end module CoreDefinition

