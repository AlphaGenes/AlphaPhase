module SurrogateDefinition
  implicit none
  private

  type, public :: SurrDef
    private
    !Almost definitely shouldn't be public but for now...
    integer(kind = 2), allocatable, dimension(:,:), public :: numOppose
    integer(kind = 1), allocatable, dimension(:,:), public :: partition
  contains
    private
    procedure, public :: calculate
  end type SurrDef

contains
  subroutine calculate(definition, genos, SireGenotyped, DamGenotyped, threshold)
    !Don't like this but for now!!
    use GlobalClustering
    !Here for input / output crap!!
    use Global, only : FullFileOutput, OutputPoint, WindowsLinux, GenotypeId
    
    class(SurrDef) :: definition
    integer(kind = 1), dimension (:,:), intent(in) :: genos
    integer(kind = 4), allocatable, dimension (:), intent(in) :: SireGenotyped, DamGenotyped
    integer :: threshold


    integer :: i, j, k, Counter, truth, PseudoDam, PseudoSire, PseudoSireOld
    integer, allocatable, dimension(:) :: Partitioned, SurrogateList, ProgCount
    integer :: CountAgreePat, CountAgreeMat, DumSire, DumDam, temp
    character(len = 300) :: filout

    integer :: ii, jj, kk, l, Noffset, Limit, Switch, CurrentGroup, ChangeThreshold, flag, FirstGroup0
    integer, allocatable :: SumsOld(:), Rank(:), SumDiff(:), SumSame(:), Group(:)
    double precision, allocatable :: Sums(:)
    logical :: NoCommonGroup

    integer :: GetnSnpErrorThreshAnims, nSurrogates
    
    ! integer,allocatable,dimension(:,:) :: nSnpErrorThreshAnims
    
    
    !!OLD GLOBAL!!
    integer :: nAnisG, nSnps
   
    nAnisG = size(genos,1)
    nSnps = size(genos,2)
    
    if (allocated(definition%numOppose)) then
      deallocate(definition%numOppose)
      deallocate(definition%partition)
    end if
    allocate(definition%numOppose(nAnisG,nAnisG))
    allocate(definition%partition(nAnisG,nAnisG))
    
    definition%partition = 0


    allocate(Partitioned(nAnisG))
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
    Partitioned = 0
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
	Partitioned(i) = 1
      end if
      
      if ((Partitioned(i) == 0).and.(SireGenotyped(i) /= 0)) then
	definition%partition(i, SireGenotyped(i)) = 1
	do j = 1, nAnisG
	  if ((i /= j).and.(definition%numOppose(SireGenotyped(i), j) <= threshold).and.&
	    (definition%numOppose(i, j) <= threshold)) then
	    definition%partition(i, j) = 1
	  endif
	enddo
        Partitioned(i) = 2
      endif

      if ((Partitioned(i) == 0).and.(DamGenotyped(i) /= 0)) then
	definition%partition(i, DamGenotyped(i)) = 2
	do j = 1, nAnisG
	  if ((i /= j).and.(definition%numOppose(DamGenotyped(i), j) <= threshold).and.&
	    (definition%numOppose(i, j) <= threshold)) then
	    definition%partition(i, j) = 2
	  endif
	enddo
	Partitioned(i) = 3
      endif
      
      if ((Partitioned(i) == 0).and.(ProgCount(i) /= 0)) then
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
	Partitioned(i) = 5
      end if
      
      if (Partitioned(i) == 2.or.Partitioned(i) == 3.or.Partitioned(i) == 4.or.Partitioned(i) == 5) then
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
      
      if (Partitioned(i) == 0) then
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
	    Partitioned(i) = 7
	  end if

	  deallocate(Medoids)
	  deallocate(ClusterMember)
	  deallocate(MinClust)
	  deallocate(TempSurrArray)
	  deallocate(TempSurrVector)
	endif
	Partitioned(i) = 6
      endif
      
      if (mod(i, 400) == 0) print*, "   Partitioning done for genotyped individual --- ", i
      definition%partition(i, i) = 0
    end do

    if (FullFileOutput == 1) then
      if (WindowsLinux == 1) then
	write (filout, '(".\Miscellaneous\Surrogates",i0,".txt")') OutputPoint
	open (unit = 13, FILE = filout, status = 'unknown')
	write (filout, '(".\Miscellaneous\SurrogatesSummary",i0,".txt")') OutputPoint
	open (unit = 19, FILE = filout, status = 'unknown')

      else
	write (filout, '("./Miscellaneous/Surrogates",i0,".txt")') OutputPoint
	open (unit = 13, FILE = filout, status = 'unknown')
	write (filout, '("./Miscellaneous/SurrogatesSummary",i0,".txt")') OutputPoint
	open (unit = 19, FILE = filout, status = 'unknown')
      endif
      do i = 1, nAnisG
	nSurrogates = 0
	if (nAnisG < 20000) then
	  write (13, '(a20,20000i6,20000i6,20000i6,20000i6)') GenotypeId(i), definition%partition(i,:)
	else
	  write (13, *) GenotypeId(i), definition%partition(i,:)
	end if
	do j = i, nAnisG
	  if (definition%numOppose(i, j) <= threshold) nSurrogates = nSurrogates + 1
	enddo
	write (19, '(a20,20000i6,20000i6,20000i6,20000i6)') &
	GenotypeId(i), count(definition%partition(i,:) == 1), count(definition%partition(i,:) == 2)&
	, count(definition%partition(i,:) == 3), nSurrogates, Partitioned(i)
      enddo
    end if

    deallocate(Partitioned)
    deallocate(SurrogateList)
    deallocate(ProgCount)
  end subroutine calculate
end module SurrogateDefinition