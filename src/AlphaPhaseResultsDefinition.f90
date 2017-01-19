module AlphaPhaseResultsDefinition
  use CoreDefinition
  use HaplotypeLibraryDefinition
  use SurrogateDefinition
  use TestResultDefinition
  use PedigreeDefinition
  implicit none
  
  private
  
  type, public :: AlphaPhaseResults
    type(Core), dimension(:), allocatable :: cores
    type(HaplotypeLibrary), dimension(:), allocatable :: libraries
    type(Surrogate), dimension(:), allocatable :: surrogates
    type(TestResults), dimension(:), allocatable :: testResults
    integer, dimension(:), allocatable :: ids
    integer, dimension(:), allocatable :: startIndexes
    integer, dimension(:), allocatable :: endIndexes
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
  end function newAlphaPhaseResults
  
end module AlphaPhaseResultsDefinition
    