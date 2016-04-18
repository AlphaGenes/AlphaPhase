module Clustering
  implicit none

contains
  function cluster(input, assign, numClusters, maxIter, multi) result (converged)
    integer, dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    integer, intent(in) :: numClusters, maxIter
    logical :: multi
    
    logical:: converged
    
    double precision, dimension(:,:), allocatable :: mediods
    integer :: i, j, iter
    logical :: swapped
    double precision :: curdist, d    
    
    allocate(mediods(numClusters,size(input,2)))
    
    iter = 1
    converged = .false.
    
    do while ((iter + 1 < maxIter) .and. (.not. converged))
      mediods = 0
      do i = 1, size(input,1)
	mediods(assign(i), :) = mediods(assign(i), :) + input(i, :)
      end do
      do i = 1, numClusters
	mediods(i,:) = mediods(i,:) / count(assign == i)
      end do
      
      i = 1
      swapped = .false.
      do while ((i <= size(input,1)) .and. ((.not. swapped) .or. multi))
	curdist = distD(mediods(assign(i),:), input(i,:))
	j = 1
	do while ((j <= numClusters) .and. ((.not. swapped) .or. multi))
	  d = distD(mediods(j, :), input(i,:))
	  if ((assign(i) /= j) .and. (d <= curdist)) then
	    assign(i) = j
	    curdist = d
	    swapped = .true.
	  end if
	  j = j + 1
	end do
	i = i + 1
      end do

      if (.not. swapped) then
	converged = .true.
      end if
      
     
      iter = iter + 1
    end do
    
    !print *, iter, maxIter, converged
    
    deallocate(mediods)
    
  end function cluster
  
  function clusterA(input, assign, numClusters, maxIter, multi) result (rounds)
    
    integer, dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    integer, intent(in) :: numClusters, maxIter
    logical :: multi
    
    integer :: j
    
    integer :: rounds
    double precision, allocatable, dimension (:,:) :: Medoids
    logical :: change
    
    allocate(Medoids(2, size(input,2)))
    call EvaluateMedoids(numClusters, input, assign, medoids)
    change = RePartition(numClusters, input, assign, medoids)
    rounds = 2
    do j = 1, maxIter
      call EvaluateMedoids(numClusters, input, assign, medoids)
      change = RePartition(numClusters, input, assign, medoids)
      rounds = rounds + 1
      if (.not. change) exit
    enddo
    deallocate(Medoids)
  end function clusterA
  
  subroutine EvaluateMedoids(numClusters, input, assign, medoids)
    
    integer, intent(in) :: numClusters
    integer, dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    double precision, allocatable, dimension (:,:) :: Medoids
    
    integer :: i, j, k, l
    
    double precision :: CountCluster(numClusters)

    Medoids = 0
    CountCluster = 0
    do i = 1, numClusters
      do k = 1, size(input,1)
	if (assign(k) == i) then
	  do l = 1, size(input,2)
	    Medoids(i, l) = Medoids(i, l) + input(k, l)
	  end do
	  CountCluster(i) = CountCluster(i) + 1
	end if
      end do	
      Medoids(i,:) = Medoids(i,:)/CountCluster(i)
    end do

  end subroutine EvaluateMedoids
  
  function RePartition(numClusters, input, assign, Medoids) result(change)
    implicit none
    
    integer, intent(in) :: numClusters
    integer, dimension(:,:), intent(in) :: input
    integer, dimension(:), intent(inout) :: assign
    double precision, allocatable, dimension (:,:) :: Medoids
    logical change

    integer :: i, j, k, l
    
    double precision :: Dist(numClusters)
    double precision, allocatable, dimension (:) :: MinClust
    
    allocate(MinClust(size(input,1)))

    change = .false.
    do i = 1, size(input,1)
      Dist = 0
      do j = 1, size(input,2)
	do k = 1, numClusters
	  Dist(k) = Dist(k) + abs(input(i, j) - Medoids(k, j))
	end do
      end do
      Dist = Dist/size(input,2)
      MinClust(i) = Dist(assign(i))
      do k = 1, numClusters
	if ((Dist(k) <= MinClust(i)).and.(assign(i) /= k)) then
	  MinClust(i) = Dist(k)
	  assign(i) = k
	  Change = .true.
	end if
      end do
    end do

    deallocate(MinClust)
  
  end function RePartition

  
  function distD(p1, p2) result(d)
    double precision, dimension(:), intent(in) :: p1
    integer, dimension(:), intent(in) :: p2
    double precision :: d
    
    d = sum(abs(p1 - p2)) / size(p1,1)
  end function distD
  
  function distg(p1, p2) result(d)
    double precision, dimension(:), intent(in) :: p1
    integer, dimension(:), intent(in) :: p2
    double precision :: d
    integer i
    
    !d = sum(abs(p1 - p2)) / size(p1,1)
    do i = 1, size(p1,1)
      d = d + abs(p1(i)-p2(i))
    end do
    d = d / size(p1,1)
  end function distg
  
end module Clustering