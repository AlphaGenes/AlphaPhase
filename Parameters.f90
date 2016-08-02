module Parameters
  character(len=300) GenotypeFile
  integer :: GenotypeFileFormat
  integer :: nSnp   ! Possibly doesn't need to be a parameter  
  integer :: CoreAndTailLength
  integer :: Jump, Offset
  integer :: UseSurrsN
  integer :: NumSurrDisagree
  double precision :: PercGenoHaploDisagree
  double precision :: GenotypeMissingErrorPercentage
  double precision :: NrmThresh
  integer :: FullFileOutput
  integer :: Simulation
  character (len = 300) :: PedigreeFile, TruePhaseFile ! Used in a really weird way that should probably be refactored
  
  logical :: readCoreAtTime
  character (len = 300) :: itterateType
  integer :: itterateNumber
  integer :: numIter
  character (len = 10) :: startCoreChar, endCoreChar
  integer :: minHapFreq
  
  logical :: consistent  
end module Parameters