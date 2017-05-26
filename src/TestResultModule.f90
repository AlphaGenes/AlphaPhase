module TestResultModule
  use ConstantModule
  use GenotypeModule
  use HaplotypeModule
  implicit none
  
  type:: TestResults
    integer, dimension(:,:,:,:), allocatable :: countA
  contains
    procedure :: percentAll
    procedure :: counts
    procedure :: percent
  end type TestResults
  
  interface TestResults
    module procedure newResults
  end interface TestResults
  
  integer, parameter :: ALL_ = 1, HET_ = 2, ERROR_ = 3, MISS_ = 4
  integer, parameter :: CORRECT_ = 1, INCORRECT_ = 2, NOTPHASED_ = 3

contains

  function newResults(c, TruePhase) result (results)
    use CoreModule
    
    type(Core) :: c
    type(Haplotype), dimension(:,:), intent(in) :: TruePhase    
    type(TestResults) :: results
    
    integer :: i, j, k
    integer(kind=1) :: p, t, g
    logical :: het, miss, error
    
    type(Genotype) :: geno
    type(Haplotype) :: hap1, hap2, trueHap1, trueHap2
    
    
    allocate(results%countA(C%getNAnisG(),2,4,3))
    results%countA = 0
    
    do i = 1, c%getNAnisG()
      geno = c%getCoreGenos(i)
      hap1 = c%phase(i,1)
      hap2 = c%phase(i,2)
      trueHap1 = TruePhase(i,1)
      trueHap2 = TruePhase(i,2)
      do j = 1, c%getNCoreSnp()
	g = geno%getGenotype(j)
	het = (g == 1)
	miss = (g == MissingGenotypeCode)
	error = ((.not. miss) .and. ( (trueHap1%getPhase(j) + trueHap2%getPhase(j)) /= g))
	do k = 1, 2
	  if (k == 1) then
	    p = hap1%getPhase(j)
	  else
	    p = hap2%getPhase(j)
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
	    if (k == 1) then
	      t = truehap1%getPhase(j)
	    else
	      t = truehap2%getPhase(j)
	    end if
	    if (p == t) then
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
    
    print '(7x, a18, 33x, a8, 3x, a8)', "Summary statistics", "Paternal", "Maternal"
    print '(7x, a24, 7x, a8, 15x, f5.1, 6x, f5.1)', "Percent correctly phased", "All Snps", &
      results%percentAll(1,ALL_,CORRECT_), results%percentAll(2,ALL_,CORRECT_)
    print '(7x, a24, 7x, a17, 6x, f5.1, 6x, f5.1)', "Percent correctly phased", "Heterozygous Snps", &
      results%percentAll(1,HET_,CORRECT_), results%percentAll(2,HET_,CORRECT_)
    print '(7x, a18, 13x, a8, 15x, f5.1, 6x, f5.1)', "Percent not phased", "All Snps", &
      results%percentAll(1,ALL_,NOTPHASED_), results%percentAll(2,ALL_,NOTPHASED_)
    print '(7x, a18, 13x, a17, 6x, f5.1, 6x, f5.1)', "Percent not phased", "Heterozygous Snps", &
      results%percentAll(1,HET_,NOTPHASED_), results%percentAll(2,HET_,NOTPHASED_)
    print '(7x, a26, 5x, a8, 15x, f5.1, 6x, f5.1)', "Percent incorrectly phased","All Snps", &
      results%percentAll(1,ALL_,INCORRECT_), results%percentAll(2,ALL_,INCORRECT_)
    print '(7x, a26, 5x, a17, 6x, f5.1, 6x, f5.1)', "Percent incorrectly phased","Heterozygous Snps", &
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

end module TestResultModule
