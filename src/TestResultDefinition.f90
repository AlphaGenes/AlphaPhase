module TestResultDefinition
  use Constants
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
    use Constants
    
    type(Core) :: c
    integer(kind=1), dimension(:,:,:), intent(in) :: TruePhase    
    type(TestResults) :: results
    
    integer :: i, j, k
    integer(kind=1) :: p, g
    logical :: het, miss, error
    
    print*, " "
    print*, " Checking simulation"
    print*, " "
    
    allocate(results%countA(C%getNAnisG(),2,4,3))
    results%countA = 0
    
    do i = 1, c%getNAnisG()
      do j = 1, c%getNCoreSnp()
	g = c%getCoreGeno(i,j)
	het = (g == 1)
	miss = (g == MissingGenotypeCode)
	error = ((.not. miss) .and. ( (TruePhase(i,j,1) + TruePhase(i,j,2)) /= g))
	do k = 1, 2
	  p = c%getPhase(i,j,k)
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
    !! This needs a new home
    use CoreDefinition

    type(Core) :: c
    integer(kind=1), dimension(:,:,:), intent(in) :: TruePhase

    integer :: i, j, CountAgreeStay1, CountAgreeStay2, CountAgreeSwitch1, CountAgreeSwitch2, truth
    integer(kind = 1), allocatable, dimension(:) :: W1, W2

    integer :: nAnisG, nSnp

    nAnisG = c%getNAnisG()
    nSnp = c%getNCoreSnp()

    allocate(W1(nSnp))
    allocate(W2(nSnp))

    do i = 1, nAnisG
      CountAgreeStay1 = 0
      CountAgreeStay2 = 0
      CountAgreeSwitch1 = 0
      CountAgreeSwitch2 = 0
      truth = 0
      do j = 1, nSnp
	if (TruePhase(i, j, 1) == c%getPhase(i, j, 1)) CountAgreeStay1 = CountAgreeStay1 + 1
	if (TruePhase(i, j, 2) == c%getPhase(i, j, 1)) CountAgreeSwitch1 = CountAgreeSwitch1 + 1
	if (TruePhase(i, j, 1) == c%getPhase(i, j, 2)) CountAgreeSwitch2 = CountAgreeSwitch2 + 1
	if (TruePhase(i, j, 2) == c%getPhase(i, j, 2)) CountAgreeStay2 = CountAgreeStay2 + 1
      end do
      if ((CountAgreeSwitch2 > CountAgreeStay2).and.(CountAgreeStay1 <= CountAgreeSwitch1)) truth = 1
      if ((CountAgreeSwitch1 > CountAgreeStay1).and.(CountAgreeStay2 <= CountAgreeSwitch2)) truth = 1
      if (truth == 1) then
	W1(:) = c%getHaplotype(i,1)
	W2(:) = c%getHaplotype(i,2)
	call c%setHaplotype(i,1,W2)
	call c%setHaplotype(i,2,W1)
      end if
    end do

  end subroutine Flipper

end module TestResultDefinition
