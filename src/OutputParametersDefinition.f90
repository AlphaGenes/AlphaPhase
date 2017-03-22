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
    
  end type OutputParameters
end module OutputParametersDefinition
