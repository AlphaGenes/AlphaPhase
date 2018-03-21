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

        integer :: outputHapCommonalityThreshold

        character(len=300) :: outputDirectory
        character(len=512) :: resultFolderPath

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
        outputParams%outputHaplotypeLibraryBinary = .true.
        outputParams%outputPhasingYield = .true.
        outputParams%outputTimer = .true.
        outputParams%outputCombined = .true.
        outputParams%outputIndivMistakes = .true.
        outputParams%outputIndivMistakesPercent = .true.
        outputParams%outputCombined = .true.
        outputParams%outputGlobalCoreMistakesPercent = .true.
        outputParams%outputHapCommonalityThreshold = 10000

    end function newOutputParameters

    function newOutputParametersImpute(binary) result(outputParams)
        type(OutputParameters) :: outputParams
        integer, optional :: binary

        outputParams%outputFinalPhase = .true.
        outputParams%outputCoreIndex = .true.
        outputParams%outputSnpPhaseRate = .false.
        outputParams%outputIndivPhaseRate = .false.
        outputParams%outputHapIndex = .false.
        outputParams%outputSwappable = .false.
        outputParams%outputHapCommonality = .false.
        outputParams%outputSurrogates = .false.
        outputParams%outputSurrogatesSummary = .false.

        if (present(binary)) then
            outputParams%outputHaplotypeLibraryText = .false.
            outputParams%outputHaplotypeLibraryBinary = .true.
        else
            outputParams%outputHaplotypeLibraryText = .true.
            outputParams%outputHaplotypeLibraryBinary = .false.
        endif
        outputParams%outputPhasingYield = .false.
        outputParams%outputTimer = .false.
        outputParams%outputCombined = .true.
        outputParams%outputIndivMistakes = .false.
        outputParams%outputIndivMistakesPercent = .false.
        outputParams%outputCoreMistakesPercent = .false.
        outputParams%outputGlobalCoreMistakesPercent = .false.
        outputParams%outputHapCommonalityThreshold = 0

    end function newOutputParametersImpute
end module OutputParametersModule
