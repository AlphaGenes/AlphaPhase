module OutputParametersModule
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
        logical :: outputPerCore
        logical :: outputCombined
        logical :: outputGlobalCoreMistakesPercent

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
        outputParams%outputCombined = .true.
        outputParams%outputGlobalCoreMistakesPercent = .true.


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
        outputParams%outputCombined = .true.
        outputParams%outputIndivMistakes = .false.
        outputParams%outputIndivMistakesPercent = .false.
        outputParams%outputCoreMistakesPercent = .false.

    end function newOutputParametersImpute
end module OutputParametersModule
