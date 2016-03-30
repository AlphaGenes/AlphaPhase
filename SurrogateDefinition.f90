module SurrogateDefinition
  implicit none
  private

  type, public :: SurrDef
    private
    !Almost definitely shouldn't be public but for now...
    integer(kind = 2), allocatable, dimension(:,:), public :: numOppose
    integer(kind = 1), allocatable, dimension(:,:), public :: partition
    integer(kind = 1), allocatable, dimension(:), public :: method
  contains
    private
    procedure, public :: calculate
  end type SurrDef

contains
  subroutine calculate(definition, genos, SireGenotyped, DamGenotyped, threshold)
    !Don't like this but for now!!
    use GlobalClustering
    use Global, only: pseudoNRM, genotypeID
    
    class(SurrDef) :: definition
    integer(kind = 1), dimension (:,:), intent(in) :: genos
    integer(kind = 4), allocatable, dimension (:), intent(in) :: SireGenotyped, DamGenotyped
    integer :: threshold
    
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

    ! integer,allocatable,dimension(:,:) :: nSnpErrorThreshAnims

    logical :: subtract1, subtract2 ! Fudge to produce same results as old version.  Effectively does old buggy logic.  Should be removed
    ! once finished with speed ups.

    nAnisG = size(genos,1)
    nSnp = size(genos,2)

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

    if (nClusters /= 2) then
      print*, "nClusters must equal 2"
      stop
    end if

    print*, " "
    print*, " Identifying surrogates"
    
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
      if (mod(i, 400) == 0) print*, "   Surrogate identification done for genotyped individual --- ", i
    end do
    
    ProgCount = 0
    do i = 1, nAnisG
      if (SireGenotyped(i) /= 0) then
	ProgCount(SireGenotyped(i)) = ProgCount(SireGenotyped(i)) + 1
      end if
      if (DamGenotyped(i) /= 0) then
	ProgCount(DamGenotyped(i)) = ProgCount(DamGenotyped(i)) + 1
      end if
    end do

    print*, " "
    print*, " Partitioning surrogates"
    do i = 1, nAnisG

      DumSire = 0
      DumDam = 0

      if ((SireGenotyped(i) /= 0).and.(DamGenotyped(i) /= 0)) then
	!do j = 1, nAnisG
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  truth = 0
	  if ((definition%numoppose((SireGenotyped(i)), j) <=  threshold) &
	    .and.(definition%numoppose((DamGenotyped(i)), j) > threshold)) then
	    definition%partition(i, j) = 1
	  endif
	  if ((definition%numoppose((DamGenotyped(i)), j) <= threshold)&
	    .and.(definition%numoppose((SireGenotyped(i)), j) > threshold)) then
	    definition%partition(i, j) = 2
	  endif
	end do
	if (definition%numoppose(i, SireGenotyped(i)) <= threshold) then
	  definition%partition(i, SireGenotyped(i)) = 1
	end if
	if (definition%numoppose(i, DamGenotyped(i)) <= threshold) then
	  definition%partition(i, DamGenotyped(i)) = 2
	end if
	definition%method(i) = 1
      end if

      if ((definition%method(i) == 0).and.(SireGenotyped(i) /= 0)) then
	definition%partition(i, SireGenotyped(i)) = 1
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  if (definition%numoppose(SireGenotyped(i), j) <= threshold) then
	    definition%partition(i, j) = 1
	  endif
	enddo
	definition%method(i) = 2
      endif

      if ((definition%method(i) == 0).and.(DamGenotyped(i) /= 0)) then
	definition%partition(i, DamGenotyped(i)) = 2
	do aj = 1, numPassThres(i)
	  j = passThres(i, aj)
	  if (definition%numoppose(DamGenotyped(i), j) <= threshold) then
	    definition%partition(i, j) = 2
	  endif
	enddo
	definition%method(i) = 3
      endif
      
      if ((definition%method(i) == 0).and.(SireGenotyped(i) == 0).and.(DamGenotyped(i) == 0)) then
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
	  if (i == DamGenotyped(j)) then
	    DumSire = j
	    exit
	  endif
	  if (i == SireGenotyped(j)) then
	    DumSire = j
	    exit
	  endif
	end do
	if (DumSire /= 0) then
	  definition%partition(i, DumSire) = 1
	  truth = 0
	  do aj = 1, numPassThres(i)
	    j = passThres(i, aj)
	    if ((i == SireGenotyped(j)).or.(i == DamGenotyped(j))) then
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

	  allocate(Medoids(nClusters, SurrCounter))
	  allocate(ClusterMember(SurrCounter))
	  allocate(MinClust(SurrCounter))
	  ClusterMember(1) = 1
	  do j = 1, SurrCounter
	    if (TempSurrArray(1, j) == 0) then
	      ClusterMember(j) = 2
	    else
	      ClusterMember(j) = 1
	    endif
	  end do
	  call EvaluateMedoids
	  Change = 0
	  MinClust = 1
	  rounds = 1
	  call RePartition
	  do j = 1, SurrCounter
	    call EvaluateMedoids
	    Change = 0
	    call RePartition
	    if (Change == 0) exit
	  enddo
	  if (rounds <= SurrCounter) then
	    do j = 1, SurrCounter
	      definition%partition(i, TempSurrVector(j)) = ClusterMember(j)
	    enddo
	    definition%method(i) = 7
	  end if

	  deallocate(Medoids)
	  deallocate(ClusterMember)
	  deallocate(MinClust)
	  deallocate(TempSurrArray)
	  deallocate(TempSurrVector)
	endif
	definition%method(i) = 6
      endif
      if (mod(i, 400) == 0) print*, "   Partitioning done for genotyped individual --- ", i
      definition%partition(i, i) = 0
    end do
  end subroutine calculate
  
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
end module SurrogateDefinition