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
    
    class(SurrDef) :: definition
    integer(kind = 1), dimension (:,:), intent(in) :: genos
    integer(kind = 4), allocatable, dimension (:), intent(in) :: SireGenotyped, DamGenotyped
    integer :: threshold


    integer :: i, j, k, Counter, truth, PseudoDam, PseudoSire, PseudoSireOld
    integer, allocatable, dimension(:) :: SurrogateList, ProgCount
    integer :: CountAgreePat, CountAgreeMat, DumSire, DumDam, temp

    integer :: ii, jj, kk, l, Noffset, Limit, Switch, CurrentGroup, ChangeThreshold, flag, FirstGroup0
    integer, allocatable :: SumsOld(:), Rank(:), SumDiff(:), SumSame(:), Group(:)
    double precision, allocatable :: Sums(:)
    logical :: NoCommonGroup

    integer :: GetnSnpErrorThreshAnims
    
    ! integer,allocatable,dimension(:,:) :: nSnpErrorThreshAnims
    
    
    !!OLD GLOBAL!!
    integer :: nAnisG, nSnps
   
    nAnisG = size(genos,1)
    nSnps = size(genos,2)
    
    if (allocated(definition%numOppose)) then
      deallocate(definition%numOppose)
      deallocate(definition%partition)
      deallocate(definition%method)
    end if
    allocate(definition%numOppose(nAnisG,nAnisG))
    allocate(definition%partition(nAnisG,nAnisG))
    allocate(definition%method(nAnisG))
    
    definition%partition = 0

    allocate(SurrogateList(nAnisG))
    allocate(ProgCount(nAnisG))

    print*, " "
    print*, " Identifying surrogates"

    do i = 1, nAnisG
      do j = i + 1, nAnisG
	Counter = 0
	do k = 1, nSnps
	  if ((Genos(i, k) == 0).and.(Genos(j, k) == 2)) Counter = Counter + 1
	  if ((Genos(i, k) == 2).and.(Genos(j, k) == 0)) Counter = Counter + 1
	end do
	definition%numOppose(i, j) = Counter
	definition%numOppose(j, i) = Counter
      end do
      definition%numOppose(i, i) = 0
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
    definition%method = 0
    do i = 1, nAnisG
   
      DumSire = 0
      DumDam = 0
      if ((SireGenotyped(i) /= 0).and.(DamGenotyped(i) /= 0)) then
	do j = 1, nAnisG
	  if ((definition%numOppose(i, j) <= threshold).and.(definition%numOppose((SireGenotyped(i)), j) <= threshold) &
	  .and.(definition%numOppose((DamGenotyped(i)), j) > threshold)) then
	    definition%partition(i, j) = 1
	  endif
	  if ((definition%numOppose(i, j) <= threshold).and.(definition%numOppose((DamGenotyped(i)), j) <= threshold)&
	  .and.(definition%numOppose((SireGenotyped(i)), j) > threshold)) then
	    definition%partition(i, j) = 2
	  endif
	end do
	if (definition%numOppose(i, SireGenotyped(i)) <= threshold) then
	  definition%partition(i, SireGenotyped(i)) = 1
	end if
	if (definition%numOppose(i, DamGenotyped(i)) <= threshold) then
	  definition%partition(i, DamGenotyped(i)) = 2
	end if
	definition%method(i) = 1
      end if
      
      if ((definition%method(i) == 0).and.(SireGenotyped(i) /= 0)) then
	definition%partition(i, SireGenotyped(i)) = 1
	do j = 1, nAnisG
	  if ((i /= j).and.(definition%numOppose(SireGenotyped(i), j) <= threshold).and.&
	    (definition%numOppose(i, j) <= threshold)) then
	    definition%partition(i, j) = 1
	  endif
	enddo
        definition%method(i) = 2
      endif

      if ((definition%method(i) == 0).and.(DamGenotyped(i) /= 0)) then
	definition%partition(i, DamGenotyped(i)) = 2
	do j = 1, nAnisG
	  if ((i /= j).and.(definition%numOppose(DamGenotyped(i), j) <= threshold).and.&
	    (definition%numOppose(i, j) <= threshold)) then
	    definition%partition(i, j) = 2
	  endif
	enddo
	definition%method(i) = 3
      endif
      
      if ((definition%method(i) == 0).and.(ProgCount(i) /= 0)) then
	DumSire = 0
	do j = 1, nAnisG
	  if ((i == DamGenotyped(j)).and.(definition%numOppose(i, j) <= threshold)) then
	    DumSire = j
	    exit
	  endif
	  if ((i == SireGenotyped(j)).and.(definition%numOppose(i, j) <= threshold)) then
	    DumSire = j
	    exit
	  endif
	end do
	if (DumSire /= 0) then
	  definition%partition(i, DumSire) = 1
	  truth = 0
	  do j = 1, nAnisG
	    if ((i == SireGenotyped(j)).or.(i == DamGenotyped(j))) then
	      if (definition%numOppose(i, j) <= threshold) then
		if (definition%numOppose(j, DumSire) > threshold) then
		  definition%partition(i, j) = 2
		  truth = 1
		  exit
		endif
	      endif
	    end if
	  end do
	end if
	definition%method(i) = 5
      end if
      
      if (definition%method(i) == 2.or.definition%method(i) == 3.or.definition%method(i) == 4.or.definition%method(i) == 5) then
	do j = 1, nAnisG
	  CountAgreePat = 0
	  CountAgreeMat = 0
	  do k = 1, nAnisG
	    if ((definition%partition(i, k) == 1).and.(definition%numOppose(i, j) <= threshold)&
	      .and.(definition%numOppose(k, j) <= threshold)) then
	      CountAgreePat = CountAgreePat + 1
	      exit !here
	    endif
	    if ((definition%partition(i, k) == 2).and.(definition%numOppose(i, j) <= threshold)&
	      .and.(definition%numOppose(k, j) <= threshold)) then
	      CountAgreeMat = CountAgreeMat + 1
	      exit !here
	    endif
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
	  if ((definition%numOppose(i, j) <= threshold).and.(i /= j)) then
	    SurrCounter = SurrCounter + 1
	  endif
	end do
	if (SurrCounter > 0) then
	  allocate(TempSurrArray(SurrCounter, SurrCounter))
	  allocate(TempSurrVector(SurrCounter))
	  SurrCounter = 0
	  do j = 1, nAnisG
	    if ((definition%numOppose(i, j) <= threshold).and.(i /= j)) then
	      SurrCounter = SurrCounter + 1
	      TempSurrVector(SurrCounter) = j
	    endif
	  end do
	  TempSurrArray = 0
	  do j = 1, SurrCounter
	    do k = 1, SurrCounter
	      if (definition%numOppose(TempSurrVector(j), TempSurrVector(k)) <= threshold) then
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

    deallocate(SurrogateList)
    deallocate(ProgCount)
  end subroutine calculate
end module SurrogateDefinition