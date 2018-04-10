module ProgramParametersModule
    use AlphaPhaseParametersModule
    use OutputParametersModule
    use baseSpecFileModule
    implicit none

    type, extends(baseSpecFile) :: ProgramParameters
        integer :: GenotypeFileFormat
        logical :: Simulation
        character (len = 300) :: TruePhaseFile, Library, CoreFile, PrePhaseFile

        type(OutputParameters) :: outputParams
        type(AlphaPhaseParameters) :: params

    end type ProgramParameters

    interface ProgramParameters
        module procedure newProgramParameters
    end interface ProgramParameters

contains
    function newProgramParameters result(programParams)
        type(ProgramParameters) :: programParams

        programParams%params = AlphaPhaseParameters()
        programParams%outputParams = OutputParameters()
        programParams%library = "None"
        programParams%CoreFile = "None"
        programParams%PrePhaseFile = "None"
        programParams%plinkInputFile = ""
        programParams%resultFolderPath = "."
        programParams%PedigreeFile = "NoPedigree"
        programParams%nSnp = 0
    end function newProgramParameters

end module ProgramParametersModule