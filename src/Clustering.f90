module Clustering
  implicit none

  interface cluster
    module procedure clusterD
    module procedure cluster1
  end interface cluster
  
contains
  function clusterD(input, assign, numClusters, maxIter, multi) result (rounds)
    
    integer, dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    integer, intent(in) :: numClusters, maxIter
    logical :: multi
    
    integer :: j
    
    integer :: rounds
    logical :: change
    
    rounds = 1
    do j = 1, maxIter + 1
      change = RePartitionD(numClusters, input, assign, multi)
      rounds = rounds + 1
      if ((.not. change) .and. (j /= 1)) exit
    enddo

  end function clusterD
  
  function RePartitionD(numClusters, input, assign, multi) result(change)
    
    integer, intent(in) :: numClusters
    integer, dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    logical, intent(in) :: multi
    
    double precision, allocatable, dimension (:,:) :: Medoids
    logical change

    integer :: i, j, k
    
    double precision :: Dist(numClusters)
    double precision minclust, d
    
    
    allocate(Medoids(numClusters, size(input,2)))
    
    Medoids = 0
    do i = 1, size(input,1)
	  Medoids(assign(i), :) = Medoids(assign(i), :) + input(i, :)
    end do	
    do i = 1, numClusters
      Medoids(i,:) = Medoids(i,:)/count(assign==i)
    end do
    
    change = .false.
    i = 1
    !do i = 1, size(input,1)
    do while ((i <= size(input,1)) .and. (multi .or. (.not. Change)))
      Dist = 0
      do j = 1, size(input,2)
	do k = 1, numClusters
	  Dist(k) = Dist(k) + abs(input(i, j) - Medoids(k, j))
	end do
      end do
      Dist = Dist/size(input,2)
      MinClust = Dist(assign(i))
      k = 1
      !do k = 1, numClusters
      do while ((k <= numClusters) .and. (multi .or. (.not. Change)))
	d = dist(k)
	if ((d <= MinClust).and.(assign(i) /= k)) then
	  MinClust = d
	  assign(i) = k
	  Change = .true.
	end if
	k = k + 1
      end do
      i = i + 1
    end do

    deallocate(Medoids)
  
  end function RePartitionD
  
  function cluster1(input, assign, numClusters, maxIter, multi) result (rounds)
    
    integer(kind=1), dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    integer, intent(in) :: numClusters, maxIter
    logical :: multi
    
    integer :: j
    
    integer :: rounds
    logical :: change
    
    rounds = 1
    do j = 1, maxIter + 1
      change = RePartition1(numClusters, input, assign, multi)
      rounds = rounds + 1
      if ((.not. change) .and. (j /= 1)) exit
    enddo

  end function cluster1
  
  function RePartition1(numClusters, input, assign, multi) result(change)
    
    integer, intent(in) :: numClusters
    integer(kind=1), dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    logical, intent(in) :: multi
    
    double precision, allocatable, dimension (:,:) :: Medoids
    logical change

    integer :: i, j, k
    
    double precision :: Dist(numClusters)
    double precision minclust, d
    
    
    allocate(Medoids(numClusters, size(input,2)))
    
    Medoids = 0
    do i = 1, size(input,1)
	  Medoids(assign(i), :) = Medoids(assign(i), :) + input(i, :)
    end do	
    do i = 1, numClusters
      Medoids(i,:) = Medoids(i,:)/count(assign==i)
    end do
    
    change = .false.
    i = 1
    !do i = 1, size(input,1)
    do while ((i <= size(input,1)) .and. (multi .or. (.not. Change)))
      Dist = 0
      do j = 1, size(input,2)
	do k = 1, numClusters
	  Dist(k) = Dist(k) + abs(input(i, j) - Medoids(k, j))
	end do
      end do
      Dist = Dist/size(input,2)
      MinClust = Dist(assign(i))
      k = 1
      !do k = 1, numClusters
      do while ((k <= numClusters) .and. (multi .or. (.not. Change)))
	d = dist(k)
	if ((d <= MinClust).and.(assign(i) /= k)) then
	  MinClust = d
	  assign(i) = k
	  Change = .true.
	end if
	k = k + 1
      end do
      i = i + 1
    end do

    deallocate(Medoids)
  
  end function RePartition1
  
  function initialAssignmentAlternate(number) result(assign)
    integer, intent(in) :: number
    integer, dimension(:), pointer :: assign
    
    integer :: i
    
    allocate(assign(number))
    do i = 1, number
      if (mod(i, 2) == 0) then
	assign(i) = 1
      else
	assign(i) = 2
      endif
    end do
  end function initialAssignmentAlternate
  
end module Clustering