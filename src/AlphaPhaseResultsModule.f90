
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaPhaseResultsModule.f90
!
! DESCRIPTION:
!> @brief    Module containing return type for alphaphase program
!> @details  Module contains information returned by alphaphase
!
!
!> @date     January 4, 2017
!
!> @version  1.0.0
!
!
!-------------------------------------------------------------------------------


module AlphaPhaseResultsModule
  use CoreModule
  use HaplotypeLibraryModule
  use SurrogateModule
  use TestResultModule
  use PedigreeModule
  implicit none


  type :: AlphaPhaseResults
    type(Core), dimension(:), allocatable :: cores
    type(HaplotypeLibrary), dimension(:), allocatable :: libraries
    type(Surrogate), dimension(:), allocatable :: surrogates
    type(TestResults), dimension(:), allocatable :: testResults
    integer, dimension(:), allocatable :: ids
    integer, dimension(:), allocatable :: startIndexes
    integer, dimension(:), allocatable :: endIndexes
    integer :: nCores

    contains
      procedure :: getFullPhase
      procedure :: getFullPhaseIntArray
  end type AlphaPhaseResults

  type :: AlphaPhaseResultsContainer
    type(AlphaPhaseResults), dimension(:), allocatable :: results
    integer, dimension(:), allocatable :: coreLengths !< core lengths correspond to results
    integer :: nResults !< number of CORELENGTHS AND RESULTS
    !< needs to be 2x the number of core lengths for offset and not offset (shifted, not) runs6
  end type AlphaPhaseResultsContainer



  interface AlphaPhaseResults
    module procedure newAlphaPhaseResults
  end interface AlphaPhaseResults

contains

  function newAlphaPhaseResults(nCores, containSurrogates, containTestResults) result(results)
    integer, intent(in) :: nCores
    logical, optional, intent(in) :: containSurrogates, containTestResults
    type(AlphaPhaseResults) :: results

    allocate(results%cores(nCores))
    allocate(results%libraries(nCores))
    if (present(containSurrogates)) then
      if (containSurrogates) then
	allocate(results%surrogates(nCores))
      end if
    end if
    if (present(containTestResults)) then
      if (containTestResults) then
	allocate(results%testResults(nCores))
      end if
    end if

    allocate(results%ids(nCores))
    allocate(results%startIndexes(nCores))
    allocate(results%endIndexes(nCores))

    results%nCores = nCores
  end function newAlphaPhaseResults


  !---------------------------------------------------------------------------
  !< @brief Returns a  2D array of full phase of haplotype objects
  !< @author  David Wilson david.wilson@roslin.ed.ac.uk
  !< @date    Febuary 26, 2017
  !---------------------------------------------------------------------------
  function getFullPhase(results) result(haps)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(AlphaPhaseResults) :: results
    type(Haplotype), dimension(:,:), pointer :: haps

    integer :: nAnisG, nSnp
    integer :: i, j, k

    nAnisG = results%cores(1)%getNAnisG()
    nSnp = 0
    do i = 1, size(results%cores)
      nSnp = nSnp + results%cores(i)%getNCoreSnp()
    end do

    allocate(haps(nAnisG,2))

    do i = 1, nAnisG
      do j = 1, 2
	haps(i,j) = Haplotype(nSnp)
	do k = 1, size(results%cores)
	  call haps(i,j)%setSubset(results%cores(k)%phase(i,j), results%startIndexes(k))
	end do
      end do
    end do
  end function getFullPhase



  !---------------------------------------------------------------------------
  !< @brief Returns a 3D Integer array of full phase
  !< @author  David Wilson david.wilson@roslin.ed.ac.uk
  !< @date    Febuary 26, 2017
  !---------------------------------------------------------------------------
  function getFullPhaseIntArray(results) result(res)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeModule
    class(AlphaPhaseResults) :: results
    integer, dimension(:,:,:), allocatable :: res !< returns array in format (nAnisG,nSnp,2)
    type(Haplotype) :: tmpHap
    integer :: nAnisG, nSnp
    integer :: i, j, k
    nAnisG = results%cores(1)%getNAnisG()
    nSnp = 0
    do i = 1, size(results%cores)
      nSnp = nSnp + results%cores(i)%getNCoreSnp()
    end do
    allocate(res(nAnisG,nSnp,2))

    do i = 1, nAnisG
      do j = 1, 2
	tmpHap = Haplotype(nSnp)
	do k = 1, size(results%cores)
	  call tmpHap%setSubset(results%cores(k)%phase(i,j), results%startIndexes(k))
	end do
  res(i,:,j) = tmpHap%toIntegerArray()
      end do
    end do
  end function getFullPhaseIntArray

! TODO need to add disk input/output subroutines for this datatype

  ! subroutine writeToFiles(this)
  !   type(AlphaPhaseResults) :: this


  ! end subroutine writeToFiles
end module AlphaPhaseResultsModule
