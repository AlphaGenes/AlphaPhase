module ProgramParametersModule
    use AlphaPhaseParametersModule
    use OutputParametersModule
    implicit none

    type:: ProgramParameters
        character(len=300) GenotypeFile
        integer :: GenotypeFileFormat
        logical :: Simulation
        character (len = 300) :: PedigreeFile, TruePhaseFile, Library, CoreFile, PrePhaseFile
        integer :: nSnp

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
        programParams%PedigreeFile = "NoPedigree"
        programParams%nSnp = 0
    end function newProgramParameters

end module ProgramParametersModule