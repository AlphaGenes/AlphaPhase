module AlphaPhaseParametersDefinition
  implicit none
  
  type:: AlphaPhaseParameters
    integer :: CoreAndTailLength
    integer :: Jump
    logical :: Offset
    integer :: UseSurrsN
    integer :: NumSurrDisagree
    integer ::  minOverlap
    double precision :: PercGenoHaploDisagree
    double precision :: GenotypeMissingErrorPercentage
    double precision :: percMinToKeep, percMinPresent

    character (len = 300) :: iterateType
    integer :: iterateNumber
    integer :: numIter
    character (len = 10) :: startCoreChar, endCoreChar
    integer :: minHapFreq
    
    integer :: tailLength
  end type AlphaPhaseParameters

  interface AlphaPhaseParameters
    module procedure newParameters
  end interface AlphaPhaseParameters
  
contains
  function newParameters result(params)
    type(AlphaPhaseParameters) :: params
    ! DEFAULT VALUES !
    params%useSurrsN = 10
    params%numSurrDisagree = 10
    params%percGenoHaploDisagree = 0
    params%genotypeMissingErrorPercentage = 0
    params%iterateType = "Off"
    params%iterateNumber = 200
    params%numIter = 1
    params%startCoreChar = "0"
    params%endCoreChar = "0"
    params%minHapFreq = 1
    params%minOverlap = 0
    params%percMinPresent = 1
    params%percMinToKeep = 1
    params%coreAndTailLength = -1
    params%tailLength = -1
  end function newParameters
    
end module AlphaPhaseParametersDefinition
