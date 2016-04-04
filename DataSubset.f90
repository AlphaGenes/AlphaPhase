module DataSubset
  implicit none
  private

  type, public :: Subset
    private
    !Almost definitely shouldn't be public but for now...
    integer(kind = 1), allocatable, dimension(:,:), public :: genos
    integer(kind = 1), allocatable, dimension(:,:,:), public :: phase
    logical, allocatable, dimension(:,:), public :: fullyPhased
    integer(kind = 4), allocatable, dimension(:), public :: sireGenotyped
    integer(kind = 4), allocatable, dimension(:), public :: damGenotyped
    integer(kind = 4), dimension(:), allocatable, public:: hapFreq
    integer, dimension(:,:,:), allocatable, public :: allHapAnis
    integer, dimension(:,:), allocatable, public :: hapAnis
    integer(kind = 4), allocatable, dimension(:), public :: full2sub
    integer(kind = 4), allocatable, dimension(:), public :: sub2full
  contains
    private
    procedure, public :: create
  end type Subset

contains

  subroutine create(set, origGenos, origPhase, origfullyPhased, origSireGenotyped, origDamGenotyped, origHapFreq, origHapAnis, &
      origAllHapAnis, members, startSnp, endSnp)
    implicit none
    
    class(Subset) :: set
    integer(kind = 1), dimension(:,:), intent(in) :: origGenos
    integer(kind = 1), dimension(:,:,:), intent(in) :: origPhase
    logical, dimension(:,:), intent(in) :: origFullyPhased
    integer(kind = 4), dimension(:), intent(in) :: origSireGenotyped, origDamGenotyped
    integer(kind = 4), dimension(:), intent(in) :: origHapFreq
    integer, dimension(:,:), intent(in) :: origHapAnis
    integer, dimension(:,:,:), intent(in) :: origAllHapAnis
    logical, dimension(:), intent(in) :: members
    integer, intent(in) :: startSnp, endSnp
    
    integer :: newNumber, origNumber, nSnp
    integer :: i, c
    
    newNumber = count(members)
    origNumber = size(origGenos,1)
    nSnp = endSnp - startSnp + 1
    
    ! I need to find a more elegant solution to this!
    if (allocated(set%genos)) then
      deallocate(set%genos)
      deallocate(set%phase)
      deallocate(set%fullyPhased)
      deallocate(set%hapFreq)
      deallocate(set%hapAnis)
      deallocate(set%allHapAnis)
      deallocate(set%sireGenotyped)
      deallocate(set%damGenotyped)
      deallocate(set%full2sub)
      deallocate(set%sub2full)
    end if
    
    allocate(set%genos(newNumber,nSnp))
    allocate(set%phase(newNumber,nSnp,2))
    allocate(set%fullyPhased(newNumber,2))
    allocate(set%hapFreq(size(origHapFreq,1)))
    allocate(set%hapAnis(newNumber,size(origHapAnis,2)))
    allocate(set%allHapAnis(newNumber,size(origAllHapAnis,2),size(origAllHapAnis,3)))
    allocate(set%sireGenotyped(newNumber))
    allocate(set%damGenotyped(newNumber))
    allocate(set%full2sub(0:origNumber))
    allocate(set%sub2full(newNumber))
    
    c = 0
    set%full2sub(0) = 0
    
    do i = 1, origNumber
      if (members(i)) then
	c = c + 1
	set%genos(c,:) = origGenos(i,startSnp:endSnp)
	set%phase(c,:,:) = origPhase(i,startSnp:endSnp,:)
	set%fullyPhased(c,:) = origFullyPhased(i,:)
	set%hapFreq(c) = origHapFreq(i)
	set%allHapAnis(c,:,:) = origAllHapAnis(i,:,:)
	set%hapAnis(c,:) = origHapAnis(i,:)
	set%full2sub(i) = c
	set%sub2full(c) = i
      else
	set%full2sub(i) = 0
      end if
    end do
    
    do i = 1, newNumber
      set%sireGenotyped(i) = set%full2sub(origSireGenotyped(set%sub2full(i)))
      set%damGenotyped(i) = set%full2sub(origDamGenotyped(set%sub2full(i)))
    end do
  end subroutine create

end module DataSubset