module CoreDefinition
  implicit none
  private

  type, public :: Core
    private
    !Almost definitely shouldn't be public but for now...
    integer(kind = 1), allocatable, dimension(:,:) :: genos
    integer(kind = 1), allocatable, dimension(:,:,:), public :: phase
    logical, allocatable, dimension(:,:), public :: fullyPhased
    integer(kind = 4), dimension(:), allocatable, public:: hapFreq
    integer, dimension(:,:), allocatable, public :: hapAnis
    
    integer :: startCoreSnp, endCoreSnp
    
    !Fudge due to oddness in old version.  Should be removed!!!!!!!
    integer :: endSurrSnp
  contains
    private
    procedure, public :: create
    procedure, public :: getCoreAndTailGenos
    procedure, public :: getCoreGenos
  end type Core

contains

  subroutine create(c, genos, startCoreSnp, endCoreSnp, endSurrSnp)
    implicit none
    
    class(Core) :: c
    integer(kind = 1), dimension(:,:), intent(in) :: genos
    integer, intent(in) :: startCoreSnp, endCoreSnp
    
    integer, intent(in) :: endSurrSnp
    
    integer :: nAnisG, nSnp, nCoreSnp
    
    nAnisG = size(genos,1)
    nSnp = size(genos,2)
    nCoreSnp = endCoreSnp - startCoreSnp + 1
    
    if (allocated(c%genos)) then
      deallocate(c%genos)
      deallocate(c%phase)
      deallocate(c%fullyPhased)
      deallocate(c%hapFreq)
      deallocate(c%hapAnis)
    end if    
    
    allocate(c%genos(nAnisG,nSnp))
    allocate(c%phase(nAnisG,nCoreSnp,2))
    allocate(c%fullyphased(nAnisG,2))
    allocate(c%hapFreq(nAnisG*2))
    allocate(c%hapAnis(nAnisG,2))
    
    c%genos = genos
    c%startCoreSnp = startCoreSnp
    c%endCoreSnp = endCoreSnp
    c%endSurrSnp = endSurrSnp
    c%phase = 9
    c%hapAnis = -99
  end subroutine create
  
  function getCoreAndTailGenos(c) result (ctGenos)
    implicit none
    
    class(Core) :: c
    integer(kind=1), dimension(:,:), allocatable :: ctGenos
    
    !allocate(ctGenos(size(c%genos,1),size(c%genos,2)))
    allocate(ctGenos(size(c%genos,1),c%endSurrSnp))
    
    !ctGenos = c%genos
    ctGenos = c%genos(:,1:c%endSurrSnp)
    
    return
  end function getCoreAndTailGenos
  
  function getCoreGenos(c) result(cGenos)
    implicit none
    
    class(Core) :: c
    integer(kind=1), dimension(:,:), allocatable :: cGenos
    
    allocate(cGenos(size(c%genos,1),c%endCoreSnp - c%startCoreSnp+1))
    
    cGenos = c%genos(:,c%startCoreSnp:c%endCoreSnp)
    
    return
  end function getCoreGenos
    
end module CoreDefinition

