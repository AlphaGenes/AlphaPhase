module OutputParametersDefinition
  implicit none

  type :: OutputParameters

    logical :: outputFinalPhase
    logical :: outputCoreIndex
    logical :: outputSnpPhaseRate
    logical :: outputIndivPhaseRate
    logical :: outputHapIndex
    logical :: outputSwappable
    logical :: outputHapCommonality
    logical :: outputSurrogates
    logical :: outputSurrogatesSummary
    logical :: outputHaplotypeLibraryText
    logical :: outputHaplotypeLibraryBinary
    logical :: outputPhasingYield
    logical :: outputTimer
    logical :: outputIndivMistakes
    logical :: outputIndivMistakesPercent
    logical :: outputCoreMistakesPercent
    logical :: outputMistakes  
    logical :: outputPerCore
    logical :: outputCombined
    
    character(len=300) :: outputDirectory
    
  end type OutputParameters
  
  interface OutputParameters
    module procedure newOutputParameters
  end interface OutputParameters
  
contains
  function newOutputParameters result(outputParams)
    type(OutputParameters) :: outputParams
    
    outputParams%outputDirectory = "."
    outputParams%outputFinalPhase = .true.
    outputParams%outputCoreIndex = .true.
    outputParams%outputSnpPhaseRate = .true.
    outputParams%outputIndivPhaseRate = .true.
    outputParams%outputHapIndex = .true.
    outputParams%outputSwappable = .true.
    outputParams%outputHapCommonality = .true.
    outputParams%outputSurrogates = .true.
    outputParams%outputSurrogatesSummary = .true.
    outputParams%outputHaplotypeLibraryText = .true.
    outputParams%outputPhasingYield = .true.
    outputParams%outputIndivMistakes = .true.
    outputParams%outputMistakes = .true.
    outputParams%outputCombined = .true.


  end function newOutputParameters

  function newOutputParametersImpute result(outputParams)
    type(OutputParameters) :: outputParams

    outputParams%outputFinalPhase = .true.
    outputParams%outputCoreIndex = .true.
    outputParams%outputSnpPhaseRate = .false.
    outputParams%outputIndivPhaseRate = .false.
    outputParams%outputHapIndex = .false.
    outputParams%outputSwappable = .false.
    outputParams%outputHapCommonality = .false.
    outputParams%outputSurrogates = .false.
    outputParams%outputSurrogatesSummary = .false.
    outputParams%outputHaplotypeLibraryText = .false.
    outputParams%outputHaplotypeLibraryBinary = .true.
    outputParams%outputPhasingYield = .false.
    outputParams%outputTimer = .false.
    outputParams%outputIndivMistakes = .false.
    outputParams%outputIndivMistakesPercent = .false.
    outputParams%outputCoreMistakesPercent = .false.
    outputParams%outputMistakes = .false.

  end function newOutputParametersImpute
end module OutputParametersDefinition
