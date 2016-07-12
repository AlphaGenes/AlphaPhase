module SurrogateDefinition
  implicit none
  private

  type, public :: Surrogate
    private
    !Almost definitely shouldn't be public but for now...
    integer(kind = 2), allocatable, dimension(:,:), public :: numOppose
    integer(kind = 1), allocatable, dimension(:,:), public :: partition
    integer(kind = 1), allocatable, dimension(:), public :: method
    integer, public :: threshold
  contains
    private
    final :: destroy
  end type Surrogate
  
  interface Surrogate
    module procedure newSurrogate
  end interface Surrogate

contains
  function newSurrogate(cs, threshold, useNRM, pseudoNRM) result(definition)
    use Clustering
    use CoreSubSetDefinition
    use Parameters, only : ItterateType
    
    class(CoreSubSet), intent(in) :: cs    
    integer, intent(in) :: threshold
    logical, intent(in) :: useNRM
    integer(kind = 1), dimension (:,:), intent(in) :: PseudoNRM
    type(Surrogate) :: definition
    
    !integer(kind = 1), dimension (:,:), allocatable :: genos
    integer(kind = 1), dimension (:,:), pointer :: genos
    
    integer, allocatable, dimension(:,:) :: passThres
    integer, allocatable, dimension(:) :: numPassThres

    integer :: nAnisG, nSnp
    
    integer :: i, j, k, Counter, truth, PseudoDam, PseudoSire, PseudoSireOld, nSnpCommon
    integer, allocatable, dimension(:) :: SurrogateList, ProgCount
    integer :: CountAgreePat, CountAgreeMat, DumSire, DumDam, temp
    character(len = 300) :: filout

    integer :: ii, jj, kk, l, Noffset, Limit, Switch, CurrentGroup, ChangeThreshold, flag, FirstGroup0
    integer, allocatable :: SumsOld(:), Rank(:), SumDiff(:), SumSame(:), Group(:)
    double precision, allocatable :: Sums(:)
    logical :: NoCommonGroup

    integer :: ai, aj, ak

    integer, parameter :: SortOrMedoid = 0 !if 1 it uses Brians Sort, If Zero it uses k-medoids

    integer :: nSurrogates

    integer :: numsections, overhang, cursection, curpos
    integer(kind = 8), allocatable, dimension(:,:) :: homo, additional

    integer :: thres, pass
    
    
    integer, allocatable, dimension (:,:) :: TempSurrArray
    integer, allocatable, dimension (:) :: TempSurrVector
    integer :: rounds, SurrCounter
    integer, allocatable, dimension (:) :: ClusterMember

    ! integer,allocatable,dimension(:,:) :: nSnpErrorThreshAnims

    logical :: subtract1, subtract2 ! Fudge to produce same results as old version.  Effectively does old buggy logic.  Should be removed
    ! once finished with speed ups.

    
    definition%threshold = threshold
    !nAnisG = size(genos,1)
    !nSnp = size(genos,2)
    nAnisG = cs%getNAnisG()
    nSnp = cs%getNCoreTailSnp()
    !allocate(genos(nAnisG,nSNp))
    genos => cs%getCoreAndTailGenos()

    allocate(SurrogateList(nAnisG))
    allocate(ProgCount(nAnisG))

    numsections = nSnp / 64 + 1
    overhang = 64 - (nSnp - (numsections - 1) * 64)

    allocate(homo(nAnisG, numsections))
    allocate(additional(nAnisG, numsections))
    homo = 0
    additional = 0

    !if (allocated(passThres)) then
    !  deallocate(passThres)
    !  deallocate(numPassThres)
    !end if
    allocate(passThres(nAnisG, nAnisG))
    allocate(numPassThres(nAnisG))
    numPassThres = 0
    passThres = 0

    do i = 1, nAnisG
      cursection = 1
      curpos = 1

      do j = 1, nSnp
	select case (Genos(i, j))
	case (0)
	  homo(i, cursection) = ibset(homo(i, cursection), curpos)
	  !additional(i,cursection) = ibclr(additional(i,cursection),curpos) NOT NEEDED AS INITIALISED TO ZERO
	case (1)
	  !homo(i,cursection) = ibclr(homo(i,cursection),curpos) NOT NEEDED AS INITIALISED TO ZERO
	  !additional(i,cursection) = ibclr(additional(i,cursection),curpos) NOT NEEDED AS INITIALISED TO ZERO
	case (2)
	  homo(i, cursection) = ibset(homo(i, cursection), curpos)
	  additional(i, cursection) = ibset(additional(i, cursection), curpos)
	case default
	  !Asuume missing
	  !homo(i,cursection) = ibclr(homo(i,cursection),curpos) NOT NEEDED AS INITIALISED TO ZERO
	  additional(i, cursection) = ibset(additional(i, cursection), curpos)
	end select
	curpos = curpos + 1
	if (curpos == 65) then
	  curpos = 1
	  cursection = cursection + 1
	end if
      end do
    end do
!
!    if (nClusters /= 2) then
!      print*, "nClusters must equal 2"
!      stop
!    end if

    if (ItterateType .eq. "Off") then
      print*, " "
      print*, " Identifying surrogates"
    end if
    
    if (allocated(definition%partition)) then
      deallocate(definition%partition)
      deallocate(definition%numoppose)
      deallocate(definition%method)
    end if
    allocate(definition%partition(nAnisG,nAnisG))
    allocate(definition%numoppose(nAnisG,nAnisG))
    allocate(definition%method(nAnisG))

    definition%partition = 0
    definition%numoppose = 0
    definition%method = 0

    do i = 1, nAnisG
      pass = 0
      do j = i + 1, nAnisG
	Counter = 0
	nSnpCommon = 0

	Counter = mismatches(homo, additional, i, j, numsections)

	definition%numoppose(i, j) = Counter
	definition%numoppose(j, i) = Counter

	if (Counter <= threshold) then
	  numPassThres(i) = numPassThres(i) + 1
	  passThres(i, numPassThres(i)) = j
	  numPassThres(j) = numPassThres(j) + 1
	  passThres(j, numPassThres(j)) = i
	end if
      end do
      definition%numoppose(i, i) = 0
      if ((mod(i, 400) == 0) .and. (ItterateType .eq. "Off")) then
	print*, "   Surrogate identification done for genotyped individual --- ", i
      end if
    end do
    
    ProgCount = 0
    do i = 1, nAnisG
      if (cs%getSire(i) /= 0) then
	ProgCount(cs%getSire(i)) = ProgCount(cs%getSire(i)) + 1
      end if
      if (cs%getDam(i) /= 0) then
	ProgCount(cs%getDam(i)) = ProgCount(cs%getDam(i)) + 1
      end if
    end do

    if (ItterateType .eq. "Off") then
      print*, " "
      print*, " Partitioning surrogates"
    end if
    do i = 1, nAnisG

      DumSire = 0
      DumDam = 0

      if ((cs%getSire(i) /= 0).and.(cs%getDam(i) /= 0)) then
	!do j = 1, nAnisG
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  truth = 0
	  if ((definition%numoppose((cs%getSire(i)), j) <=  threshold) &
	    .and.(definition%numoppose((cs%getDam(i)), j) > threshold)) then
	    definition%partition(i, j) = 1
	  endif
	  if ((definition%numoppose((cs%getDam(i)), j) <= threshold)&
	    .and.(definition%numoppose((cs%getSire(i)), j) > threshold)) then
	    definition%partition(i, j) = 2
	  endif
	end do
	if (definition%numoppose(i, cs%getSire(i)) <= threshold) then
	  definition%partition(i, cs%getSire(i)) = 1
	end if
	if (definition%numoppose(i, cs%getDam(i)) <= threshold) then
	  definition%partition(i, cs%getDam(i)) = 2
	end if
	definition%method(i) = 1
      end if

      if ((definition%method(i) == 0).and.(cs%getSire(i) /= 0)) then
	definition%partition(i, cs%getSire(i)) = 1
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  if (definition%numoppose(cs%getSire(i), j) <= threshold) then
	    definition%partition(i, j) = 1
	  endif
	enddo
	definition%method(i) = 2
      endif

      if ((definition%method(i) == 0).and.(cs%getDam(i) /= 0)) then
	definition%partition(i, cs%getDam(i)) = 2
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  if (definition%numoppose(cs%getDam(i), j) <= threshold) then
	    definition%partition(i, j) = 2
	  endif
	enddo
	definition%method(i) = 3
      endif
      
      if ((definition%method(i) == 0).and.(cs%getSire(i) == 0).and.(cs%getDam(i) == 0).and.useNRM) then
	subtract1 = .false.
	! Same fudge as below.  Definitely should not be here
	do aj = 1, numPassThres(i) + 1
	  if (.not.subtract1) then
	  j = passThres(i, aj)
	  if ((j > i) .or. (j == 0)) then
	    j = i
	    subtract1 = .true.
	  end if
	  else
	    j = passThres(i, aj - 1)
	  end if
	  if (PseudoNRM(i, j) == 1) then
	    DumSire = j
	    definition%partition(i, j) = 1
	    definition%method(i) = 4
	    exit
	  endif
	end do
	subtract1 = .false.
	do aj = 1, numPassThres(i) + 1
	  if (.not.subtract1) then
	    j = passThres(i, aj)
	    if ((j > i) .or. (j == 0)) then
	      j = i
	      subtract1 = .true.
	    end if
	  else
	    j = passThres(i, aj - 1)
	  end if
	  if (PseudoNRM(i, j) == 2) then
	    DumDam = j
	    definition%partition(i, j) = 2
	    definition%method(i) = 4
	    exit
	  end if
	end do
      endif
      
      if ((definition%method(i) == 0).and.(ProgCount(i) /= 0)) then
	DumSire = 0
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  if (i == cs%getDam(j)) then
	    DumSire = j
	    exit
	  endif
	  if (i == cs%getSire(j)) then
	    DumSire = j
	    exit
	  endif
	end do
	if (DumSire /= 0) then
	  definition%partition(i, DumSire) = 1
	  truth = 0
	  do aj = 1, numPassThres(i)
	    j = passThres(i, aj)
	    if ((i == cs%getSire(j)).or.(i == cs%getDam(j))) then
	      if (definition%numoppose(j, DumSire) > threshold) then
		definition%partition(i, j) = 2
		truth = 1
		exit
	      endif
	    end if
	  end do
	end if
	definition%method(i) = 5
      end if

      if (definition%method(i) > 1) then
	subtract1 = .false.
	do aj = 1, numPassThres(i) + 1
	  ! Fudge to get same (buggy) result as before.  Should really be just:
	  ! j = passThres(i,aj)
	  if (.not.subtract1) then
	    j = passThres(i, aj)
	    if ((j > i) .or. (j == 0)) then
	      j = i
	      subtract1 = .true.
	    end if
	  else
	    j = passThres(i, aj - 1)
	  end if
	  CountAgreePat = 0
	  CountAgreeMat = 0
	  subtract2 = .false.
	  do ak = 1, numPassThres(j) + 1
!	  do k = 1, nAnisG
	    !More fudges!
	    k = passThres(j, ak)
	    if (.not.subtract2) then
	      k = passThres(j, ak)
	      if ((k > j) .or. (k == 0)) then
		k = j
		subtract2 = .true.
	      end if
	    else
	      k = passThres(j, ak - 1)
	    end if
!	    if (definition%numoppose(j,k) == 0) then
	      !if (i == 20) print *, j, k
	      if (definition%partition(i, k) == 1) then
		CountAgreePat = CountAgreePat + 1
		exit !here
	      endif
	      if (definition%partition(i, k) == 2) then
		CountAgreeMat = CountAgreeMat + 1
		exit !here
	      endif
!	    endif
	  end do
	  if ((CountAgreePat /= 0).and.(CountAgreeMat == 0)) then
	    definition%partition(i, j) = 1
	  end if
	  if ((CountAgreePat == 0).and.(CountAgreeMat /= 0)) then
	    definition%partition(i, j) = 2
	  end if
	end do
      end if
      
      if (definition%method(i) == 0) then
	SurrCounter = 0
	do j = 1, nAnisG
	  if ((definition%numoppose(i, j) <= threshold).and.(i /= j)) then
	    SurrCounter = SurrCounter + 1
	  endif
	end do
	if (SurrCounter > 0) then
	  allocate(TempSurrArray(SurrCounter, SurrCounter))
	  allocate(TempSurrVector(SurrCounter))
	  SurrCounter = 0
	  do j = 1, nAnisG
	    if ((definition%numoppose(i, j) <= threshold).and.(i /= j)) then
	      SurrCounter = SurrCounter + 1
	      TempSurrVector(SurrCounter) = j
	    endif
	  end do
	  TempSurrArray = 0
	  do j = 1, SurrCounter
	    do k = 1, SurrCounter
	      if (definition%numoppose(TempSurrVector(j), TempSurrVector(k)) <= threshold) then
		TempSurrArray(j, k) = 1
	      end if
	    end do
	    TempSurrArray(j, j) = 1
	  end do

	  allocate(ClusterMember(SurrCounter))
	  ClusterMember(1) = 1
	  do j = 1, SurrCounter
	    if (TempSurrArray(1, j) == 0) then
	      ClusterMember(j) = 2
	    else
	      ClusterMember(j) = 1
	    endif
	  end do
	  rounds = cluster(TempSurrArray, ClusterMember, 2, SurrCounter, .true.)
	  if (rounds <= SurrCounter) then
	    do j = 1, SurrCounter
	      definition%partition(i, TempSurrVector(j)) = ClusterMember(j)
	    enddo
	    definition%method(i) = 7
	  end if

	  deallocate(ClusterMember)
	  deallocate(TempSurrArray)
	  deallocate(TempSurrVector)
	endif
	definition%method(i) = 6
      endif
      if ((mod(i, 400) == 0) .and. (ItterateType .eq. "Off")) then
	print*, "   Partitioning done for genotyped individual --- ", i
      end if
      definition%partition(i, i) = 0
    end do
  end function newSurrogate
  
  function mismatches(homo, additional, first, second, numsections) result(c)
    integer(kind = 8), dimension(:,:), intent(in) :: homo, additional
    integer, intent(in) :: first, second, numsections
    integer :: c, i

    c = 0
    do i = 1, numsections
	c = c + POPCNT(IAND(IAND(homo(first, i), homo(second, i)), &
	IEOR(additional(first, i), additional(second, i))))
    end do
  end function mismatches
  
  subroutine destroy(definition)
    type(Surrogate) :: definition
    
    if (allocated(definition%partition)) then
      deallocate(definition%partition)
      deallocate(definition%numoppose)
      deallocate(definition%method)
    end if
    
  end subroutine destroy
    
end module SurrogateDefinition