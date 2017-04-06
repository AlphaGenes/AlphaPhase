module ProgramParametersDefinition
  use AlphaPhaseParametersDefinition
  use OutputParametersDefinition
  implicit none
  
  type:: ProgramParameters
    character(len=300) GenotypeFile
    integer :: GenotypeFileFormat
    logical :: Simulation
    character (len = 300) :: PedigreeFile, TruePhaseFile, Library, CoreFile
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
  end function newProgramParameters
    
end module ProgramParametersDefinition