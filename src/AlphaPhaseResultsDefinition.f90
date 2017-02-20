module AlphaPhaseResultsDefinition
  use CoreDefinition
  use HaplotypeLibraryDefinition
  use SurrogateDefinition
  use TestResultDefinition
  use PedigreeModule
  implicit none
  
  
  type :: AlphaPhaseResults
    type(Core), dimension(:), allocatable :: cores
    type(HaplotypeLibrary), dimension(:), allocatable :: libraries
    type(Surrogate), dimension(:), allocatable :: surrogates
    type(TestResults), dimension(:), allocatable :: testResults
    integer, dimension(:), allocatable :: ids
    integer, dimension(:), allocatable :: startIndexes
    integer, dimension(:), allocatable :: endIndexes
    integer :: nCores
    
    contains 
      procedure :: getFullPhase
  end type AlphaPhaseResults
  
  interface AlphaPhaseResults
    module procedure newAlphaPhaseResults
  end interface AlphaPhaseResults
  
contains

  function newAlphaPhaseResults(nCores, containSurrogates, containTestResults) result(results)
    integer, intent(in) :: nCores
    logical, optional, intent(in) :: containSurrogates, containTestResults
    type(AlphaPhaseResults) :: results
    
    allocate(results%cores(nCores))
    allocate(results%libraries(nCores))
    if (present(containSurrogates)) then
      if (containSurrogates) then
	allocate(results%surrogates(nCores))
      end if
    end if
    if (present(containTestResults)) then
      if (containTestResults) then
	allocate(results%testResults(nCores))
      end if
    end if
    
    allocate(results%ids(nCores))
    allocate(results%startIndexes(nCores))
    allocate(results%endIndexes(nCores))
    
    results%nCores = nCores
  end function newAlphaPhaseResults
  
  function getFullPhase(results) result(haps)
    class(AlphaPhaseResults) :: results
    type(Haplotype), dimension(:,:), pointer :: haps
    
    integer :: nAnisG, nSnp
    integer :: i, j, k
    
    nAnisG = results%cores(1)%getNAnisG()
    nSnp = 0
    do i = 1, size(results%cores)
      nSnp = nSnp + results%cores(i)%getNCoreSnp()
    end do
    
    allocate(haps(nAnisG,2))
    
    do i = 1, nAnisG
      do j = 1, 2
	haps(i,j) = Haplotype(nSnp)
	do k = 1, size(results%cores)
	  call haps(i,j)%setSubset(results%cores(k)%phase(i,j), results%startIndexes(k))
	end do
      end do
    end do
  end function getFullPhase
  
end module AlphaPhaseResultsDefinition
    