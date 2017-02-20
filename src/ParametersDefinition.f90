module AlphaPhaseParametersDefinition
  implicit none
  
  type:: AlphaPhaseParameters
    integer :: CoreAndTailLength
    integer :: Jump
    logical :: Offset
    integer :: UseSurrsN
    integer :: NumSurrDisagree
    double precision :: PercGenoHaploDisagree
    double precision :: GenotypeMissingErrorPercentage

    character (len = 300) :: itterateType
    integer :: itterateNumber
    integer :: numIter
    character (len = 10) :: startCoreChar, endCoreChar
    integer :: minHapFreq
  end type AlphaPhaseParameters

  interface AlphaPhaseParameters
    module procedure newParameters
  end interface AlphaPhaseParameters
  
contains
  function newParameters result(params)
    type(AlphaPhaseParameters) :: params
    !!! DEFAULT VALUES !!!
    params%itterateType = "Off"
    params%itterateNumber = 200
    params%numIter = 1
    params%startCoreChar = "1"
    params%endCoreChar = "Combine"
    params%minHapFreq = 1
  end function newParameters
    
end module AlphaPhaseParametersDefinition