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
    params%iterateType = "Off"
    params%iterateNumber = 200
    params%numIter = 1
    params%startCoreChar = "0"
    params%endCoreChar = "0"
    params%minHapFreq = 1
    params%minOverlap = 0
    params%percMinPresent = 100
    params%percMinToKeep = 100
    params%coreAndTailLength = -1
    params%tailLength = -1
  end function newParameters
    
end module AlphaPhaseParametersDefinition
