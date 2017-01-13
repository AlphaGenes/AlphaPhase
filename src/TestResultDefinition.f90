module TestResultDefinition
  use ConstantModule
  use GenotypeModule
  use HaplotypeModule
  implicit none
!  private
  
  type, public :: TestResults
    private
    
    integer, dimension(:,:,:,:), allocatable :: countA
  contains
    procedure, public :: percentAll
    procedure, public :: counts
    procedure, public :: percent
  end type TestResults
  
  interface TestResults
    module procedure newResults
  end interface TestResults
  
  integer, parameter, public :: ALL_ = 1, HET_ = 2, ERROR_ = 3, MISS_ = 4
  integer, parameter, public :: CORRECT_ = 1, INCORRECT_ = 2, NOTPHASED_ = 3

contains

  function newResults(c, TruePhase) result (results)
    use CoreDefinition
    
    type(Core) :: c
    integer(kind=1), dimension(:,:,:), intent(in) :: TruePhase    
    type(TestResults) :: results
    
    integer :: i, j, k
    integer(kind=1) :: p, g
    logical :: het, miss, error
    
    type(Genotype), pointer :: geno
    type(Haplotype), pointer:: hap1, hap2
    
    print*, " "
    print*, " Checking simulation"
    print*, " "
    
    allocate(results%countA(C%getNAnisG(),2,4,3))
    results%countA = 0
    
    do i = 1, c%getNAnisG()
      geno => c%getCoreGenos(i)
      hap1 => c%phase(i,1)
      hap2 => c%phase(i,1)
      do j = 1, c%getNCoreSnp()
	g = geno%getGenotype(j)
	het = (g == 1)
	miss = (g == MissingGenotypeCode)
	error = ((.not. miss) .and. ( (TruePhase(i,j,1) + TruePhase(i,j,2)) /= g))
	do k = 1, 2
	  if (k == 1) then
	    p = hap1%getPhaseMod(j)
	  else
	    p = hap2%getPhaseMod(j)
	  end if
	  if (p == MissingPhaseCode) then
	    results%countA(i,k,ALL_,NOTPHASED_) = results%countA(i,k,ALL_,NOTPHASED_) + 1
	    if (het) then
	      results%countA(i,k,HET_,NOTPHASED_) = results%countA(i,k,HET_,NOTPHASED_) + 1
	    end if
	    if (miss) then
	      results%countA(i,k,MISS_,NOTPHASED_) = results%countA(i,k,MISS_,NOTPHASED_) + 1
	    end if
	    if (error) then
	      results%countA(i,k,ERROR_,NOTPHASED_) = results%countA(i,k,ERROR_,NOTPHASED_) + 1
	    end if
	  else
	    if (p == TruePhase(i,j,k)) then
	      results%countA(i,k,ALL_,CORRECT_) = results%countA(i,k,ALL_,CORRECT_) + 1
	      if (het) then
		results%countA(i,k,HET_,CORRECT_) = results%countA(i,k,HET_,CORRECT_) + 1
	      end if
	      if (miss) then
		results%countA(i,k,MISS_,CORRECT_) = results%countA(i,k,MISS_,CORRECT_) + 1
	      end if
	      if (error) then
		results%countA(i,k,ERROR_,CORRECT_) = results%countA(i,k,ERROR_,CORRECT_) + 1
	      end if
	    else
	      results%countA(i,k,ALL_,INCORRECT_) = results%countA(i,k,ALL_,INCORRECT_) + 1
	      if (het) then
		results%countA(i,k,HET_,INCORRECT_) = results%countA(i,k,HET_,INCORRECT_) + 1
	      end if
	      if (miss) then
		results%countA(i,k,MISS_,INCORRECT_) = results%countA(i,k,MISS_,INCORRECT_) + 1
	      end if
	      if (error) then
		results%countA(i,k,ERROR_,INCORRECT_) = results%countA(i,k,ERROR_,INCORRECT_) + 1
	      end if
	    end if
	  end if
	end do
      end do
    end do
    
    print *, "   Summary statistics                               ", "Paternal  Maternal"
    write (*, '(a40,a12,2f8.1)') "Percent correctly phased    All Snps", " ", &
      results%percentAll(1,ALL_,CORRECT_), results%percentAll(2,ALL_,CORRECT_)
    write (*, '(a49,a3,2f8.1)') "Percent correctly phased    Heterozygous Snps", " ", &
      results%percentAll(1,HET_,CORRECT_), results%percentAll(2,HET_,CORRECT_)
    write (*, *) " "
    write (*, '(a34,a18,2f8.1)') "Percent not phased    All Snps", " ", &
      results%percentAll(1,ALL_,NOTPHASED_), results%percentAll(2,ALL_,NOTPHASED_)
    write (*, '(a43,a9,2f8.1)') "Percent not phased    Heterozygous Snps", " ", &
      results%percentAll(1,HET_,NOTPHASED_), results%percentAll(2,HET_,NOTPHASED_)
    write (*, *) " "
    write (*, '(a42,a10,2f8.1)') "Percent incorrectly phased    All Snps", " ", &
      results%percentAll(1,ALL_,INCORRECT_), results%percentAll(2,ALL_,INCORRECT_)
    write (*, '(a51,a1,2f8.1)') "Percent incorrectly phased    Heterozygous Snps", " ", &
      results%percentAll(1,HET_,INCORRECT_), results%percentAll(2,HET_,INCORRECT_)
  end function newResults
  
  function percentAll(results, phase, group, state) result(p)
    class(TestResults) :: results
    integer, intent(in) :: phase,group,state
    integer :: s
    double precision :: p

    s = sum(results%countA(:,phase,group,:))
    if (s /= 0) then
      p = 100 * float(sum(results%countA(:,phase,group,state))) / s
    else
      p = 0.0
    end if
  end function percentAll
  
  function percent(results, animal, phase, group, state) result(p)
    class(TestResults) :: results
    integer, intent(in) :: animal,phase,group,state
    integer :: s
    double precision :: p
    
    s = sum(results%countA(animal,phase,group,:))
    if (s /= 0) then
      p = 100 * float(results%countA(animal,phase,group,state)) / s
    else
      p = 0.0
    end if
  end function percent
  
  function counts(results, animal, phase, group, state) result(c)
    class(TestResults) :: results
    integer, intent(in) :: animal,phase,group,state
    integer:: c

    c = results%countA(animal,phase,group,state)
  end function counts 

  subroutine Flipper(c, TruePhase)
    use HaplotypeModule
    !! This needs a new home
    use CoreDefinition

    type(Core) :: c
    integer(kind=1), dimension(:,:,:), intent(in) :: TruePhase

    integer :: i, j, CountAgreeStay1, CountAgreeStay2, CountAgreeSwitch1, CountAgreeSwitch2, truth
    type(Haplotype), pointer :: W1, W2
    type(Haplotype), pointer :: hap1, hap2

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
      hap2 => c%phase(i,1)
      do j = 1, nSnp
	if (TruePhase(i, j, 1) == hap1%getPhaseMod(j)) CountAgreeStay1 = CountAgreeStay1 + 1
	if (TruePhase(i, j, 2) == hap1%getPhaseMod(j)) CountAgreeSwitch1 = CountAgreeSwitch1 + 1
	if (TruePhase(i, j, 1) == hap2%getPhaseMod(j)) CountAgreeSwitch2 = CountAgreeSwitch2 + 1
	if (TruePhase(i, j, 2) == hap2%getPhaseMod(j)) CountAgreeStay2 = CountAgreeStay2 + 1
      end do
      if ((CountAgreeSwitch2 > CountAgreeStay2).and.(CountAgreeStay1 <= CountAgreeSwitch1)) truth = 1
      if ((CountAgreeSwitch1 > CountAgreeStay1).and.(CountAgreeStay2 <= CountAgreeSwitch2)) truth = 1
      if (truth == 1) then
	W1 => c%phase(i,1)
	W2 => c%phase(i,2)
	call c%setHaplotype(i,1,W2)
	call c%setHaplotype(i,2,W1)
      end if
    end do

  end subroutine Flipper

end module TestResultDefinition
