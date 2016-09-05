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
    logical :: FullFileOutput
    logical :: Simulation
    character (len = 300) :: PedigreeFile, TruePhaseFile ! Used in a really weird way that should probably be refactored

    logical :: readCoreAtTime
    character (len = 300) :: itterateType
    integer :: itterateNumber
    integer :: numIter
    character (len = 10) :: startCoreChar, endCoreChar
    integer :: minHapFreq

    logical :: consistent
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