module ProgramParametersDefinition
  use ParametersDefinition
  implicit none
  
  type:: ProgramParameters
    character(len=300) GenotypeFile
    integer :: GenotypeFileFormat
    logical :: Simulation
    character (len = 300) :: PedigreeFile, TruePhaseFile, Library
    integer :: nSnp

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
    
    type(Parameters) :: params
    
  end type ProgramParameters

  interface ProgramParameters
    module procedure newProgramParameters
  end interface ProgramParameters
  
contains
  function newProgramParameters result(programParams)
    type(ProgramParameters) :: programParams
    
    programParams%params = Parameters()
    programParams%library = "None"
  end function newProgramParameters
    
end module ProgramParametersDefinition