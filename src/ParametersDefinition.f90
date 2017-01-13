module ParametersDefinition
  implicit none
  private
  
  type, public:: Parameters
    !! Having all these as public and having them accessed directly is probably a bad idea but is easiest for now until we do
    !! something consistent across the AlphaSuite
    character(len=300) GenotypeFile
    integer :: GenotypeFileFormat
    integer :: nSnp   ! Possibly doesn't need to be a parameter  
    integer :: CoreAndTailLength
    integer :: Jump
    logical :: Offset
    integer :: UseSurrsN
    integer :: NumSurrDisagree
    double precision :: PercGenoHaploDisagree
    double precision :: GenotypeMissingErrorPercentage
    double precision :: NrmThresh
    logical :: Simulation
    character (len = 300) :: PedigreeFile, TruePhaseFile, Library ! Used in a really weird way that should probably be refactored

    logical :: readCoreAtTime
    character (len = 300) :: itterateType
    integer :: itterateNumber
    integer :: numIter
    character (len = 10) :: startCoreChar, endCoreChar
    integer :: minHapFreq
    
    integer :: nChips
    character(len = 300) :: ChipsSnps, ChipsAnimals
       
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
    
  end type Parameters

  interface Parameters
    module procedure newParameters
  end interface Parameters
  
contains
  function newParameters result(params)
    type(Parameters) :: params
    params%readCoreAtTime = .false.
    !! Place holder - should probably contain lots of defaults !!
  end function newParameters
    
end module ParametersDefinition