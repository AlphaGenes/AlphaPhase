module Phasing
  implicit none
  
  contains
  
  subroutine Erdos(surrogates, threshold, c)
  use SurrogateDefinition
  use CoreSubsetDefinition
  use Parameters, only: useSurrsN, consistent
  implicit none

  type(SurrDef), intent(in) :: surrogates
  type(CoreSubSet) :: c
  integer(kind=1), dimension(:,:), pointer :: genos
  integer, intent(in) :: threshold
    
  integer(kind = 1), allocatable, dimension (:) :: Visited
  integer, allocatable, dimension (:) :: SurrAveDiff
  integer :: AlleleCount(2)
  integer :: ErdosNumber, HighestErdos
  
  integer :: nAnisG, nCoreSnp
  
  integer :: i, j
  integer :: counter, SizeCore
  double precision :: value

  nAnisG = c%getNAnisG()
  nCoreSnp = c%getNCoreSnp()
  
  allocate(Visited(nAnisG))
  allocate(SurrAveDiff(nAnisG))
  
  genos => c%getCoreGenos()

  do i = 1, nAnisG
    value = 0
    counter = 0
    do j = 1, nAnisG
      if (surrogates%numOppose(i, j) > threshold) then
	value = value+surrogates%numOppose(i, j)
	counter = counter + 1
      endif
    end do
    SurrAveDiff(i) = int(value/counter)
  end do

  HighestErdos = 1
  print*, " "
  print*, " Phasing genotyped individuals for Paternal allele"

  do i = 1, nAnisG
    if (mod(i, 400) == 0) print*, "   Phasing done for genotyped individual --- ", i
    do j = 1, nCoreSnp
      AlleleCount = 0
      Visited = 0
      Visited(i) = 1
      if (Genos(i, j) == 0) then
	AlleleCount(1) = AlleleCount(1) + 1
      elseif (Genos(i, j) == 2) then
	AlleleCount(2) = AlleleCount(2) + 1
      endif
      if (sum(AlleleCount(:)) == UseSurrsN) then
	if (Genos(i, j) == 0) then
	  call c%setPhase(i, j, 1, 0)
	elseif (Genos(i, j) == 2) then
	  call c%setPhase(i, j, 1, 1)
	endif
      else
	call c%setPhase(i, j, 1, IterAllele(i, j, 1, surrogates, threshold, visited, allelecount, erdosnumber, genos, surravediff))
	if (ErdosNumber > HighestErdos) HighestErdos = ErdosNumber
      end if
    end do
  end do

  print*, " "
  print*, " ", HighestErdos, " was the highest Erdos used on Paternal Side"


  HighestErdos = 1
  print*, " "
  print*, " Phasing genotyped individuals for Maternal allele"

  do i = 1, nAnisG
    if (mod(i, 400) == 0) print*, "   Phasing done for genotyped individual --- ", i
    do j = 1, nCoreSnp
      AlleleCount = 0
      Visited = 0
      Visited(i) = 1
      if (Genos(i, j) == 0) then
	AlleleCount(1) = AlleleCount(1) + 1
      elseif (Genos(i, j) == 2) then
	AlleleCount(2) = AlleleCount(2) + 1
      endif
      if (sum(AlleleCount(:)) == UseSurrsN) then
	if (Genos(i, j) == 0) then
	  call c%setPhase(i, j, 2, 0)
	elseif (Genos(i, j) == 2) then
	  call c%setPhase(i, j, 2, 1)
	endif
      else
	call c%setPhase(i, j, 2, IterAllele(i, j, 2, surrogates, threshold, visited, allelecount, erdosnumber, genos, surravediff))
	if (ErdosNumber > HighestErdos) HighestErdos = ErdosNumber
      end if
    end do
  end do
  
  if (.not. consistent) then
    do i = 1, nAnisG
      do j = 1, nCoreSnp
	if (genos(i,j) == 0) then
	  call c%setPhase(i,j,1,0)
	  call c%setPhase(i,j,2,0)
	end if
	if (genos(i,j) == 2) then
	  call c%setPhase(i,j,1,1)
	  call c%setPhase(i,j,2,1)
	end if
      end do
    end do
  end if
  
  print*, " "
  print*, " ", HighestErdos, " was the highest Erdos used on Maternal Side"

  deallocate(Visited)
  deallocate(SurrAveDiff)

end subroutine Erdos

!#################################################################################################################################################################

function IterAllele(animal, snp, SideOn, surrogates, threshold, visited, allelecount, erdosnumber, genos, surravediff) result (iAllele)
  use SurrogateDefinition
  use Parameters, only: UseSurrsN, NumSurrDisagree
  use Constants
  implicit none
  
  type(SurrDef), intent(in) :: surrogates
  integer :: threshold
  integer(kind = 1), allocatable, dimension (:) :: Visited
  integer :: AlleleCount(2)
  integer :: erdosnumber
  integer(kind=1), dimension(:,:), intent(in) :: genos
  integer,dimension(:), intent(in) :: SurrAveDiff
  
  integer(kind = 1) :: iAllele

  integer :: i, j, animal, snp, SideOn
  integer(kind = 1), allocatable, dimension(:) :: ErdosNowVec, ErdosNextVec

  integer :: nAnisG

  nAnisG = size(genos,1)
  allocate(ErdosNowVec(nAnisG))
  allocate(ErdosNextVec(nAnisG))

  ErdosNumber = 1
  do i = 1, nAnisG
    if ((surrogates%numOppose(animal, i) <= threshold).and.(surrogates%partition(animal, i) /= SideOn)) then
      Visited(i) = 1
    endif
    if ((surrogates%numOppose(animal, i) > threshold).and.(surrogates%numOppose(animal, i) <= SurrAveDiff(i))) then
      Visited(i) = 1
    end if
  enddo
  
  ErdosNowVec = 0
  ErdosNextVec = 0
  do i = 1, nAnisG
    if ((surrogates%partition(animal, i) == SideOn).and.(Visited(i) /= 1)) then
      if (genos(i, snp) == 0) then
	Visited(i) = 1
	AlleleCount(1) = AlleleCount(1) + 1
      end if
      if (genos(i, snp) == 2) then
	Visited(i) = 1
	AlleleCount(2) = AlleleCount(2) + 1
      end if
      if (sum(AlleleCount(:)) == UseSurrsN) then
	if ((AlleleCount(1) > 0).and.(AlleleCount(2) <= NumSurrDisagree)) then
	  iAllele = 0
	  deallocate(ErdosNowVec)
	  deallocate(ErdosNextVec)
	  return
	elseif ((AlleleCount(2) > 0).and.(AlleleCount(1) <= NumSurrDisagree)) then
	  iAllele = 1
	  deallocate(ErdosNowVec)
	  deallocate(ErdosNextVec)
	  return
	elseif ((AlleleCount(2) > NumSurrDisagree).and.(AlleleCount(1) > NumSurrDisagree)) then
	  iAllele = 9
	  deallocate(ErdosNowVec)
	  deallocate(ErdosNextVec)
	  return
	end if
      end if
      if ((genos(i, snp) == MissingGenotypeCode).or.(genos(i, snp) == 1)) then
	Visited(i) = 1
	ErdosNextVec(i) = 1
      endif
    endif
  enddo
  
  do while (sum(ErdosNextVec) /= 0)
    ErdosNumber = ErdosNumber + 1
    ErdosNowVec = ErdosNextVec
    ErdosNextVec = 0
    do i = 1, nAnisG
      if (ErdosNowVec(i) /= 0) then
	do j = 1, nAnisG
	  if (surrogates%numOppose(i, j) > threshold) then
	    if (surrogates%numOppose(i, j) <= SurrAveDiff(j)) Visited(j) = 1
	  else
	    if (Visited(j) /= 1) then
	    Visited(j) = 1
	    if (genos(j, snp) == 0) then
	      if (mod(ErdosNumber, 2) == 0) then
		AlleleCount(2) = AlleleCount(2) + 1
	      end if
	      if (mod(ErdosNumber, 2) /= 0) then
		AlleleCount(1) = AlleleCount(1) + 1
	      end if
	    end if
	    if (genos(j, snp) == 2) then
	      if (mod(ErdosNumber, 2) == 0) then
		AlleleCount(1) = AlleleCount(1) + 1
	      end if
	      if (mod(ErdosNumber, 2) /= 0) then
		AlleleCount(2) = AlleleCount(2) + 1
	      end if
	    end if
	    if (sum(AlleleCount(:)) == UseSurrsN) then
	      if ((AlleleCount(1) > 0).and.(AlleleCount(2) <= NumSurrDisagree)) then
		iAllele = 0
		deallocate(ErdosNowVec)
		deallocate(ErdosNextVec)
		return
	      elseif ((AlleleCount(2) > 0).and.(AlleleCount(1) <= NumSurrDisagree)) then
		iAllele = 1
		deallocate(ErdosNowVec)
		deallocate(ErdosNextVec)
		return
	      elseif ((AlleleCount(2) > NumSurrDisagree).and.(AlleleCount(1) > NumSurrDisagree)) then
		iAllele = 9
		deallocate(ErdosNowVec)
		deallocate(ErdosNextVec)
		return
	      end if
	    end if
	    if ((genos(j, snp) == MissingGenotypeCode).or.(genos(j, snp) == 1)) then
	      ErdosNextVec(j) = 1
	    endif
	    end if
	  endif
	end do
      end if
    end do
  end do
  iAllele = 9
  if (sum(AlleleCount(:)) < UseSurrsN) then
    if ((AlleleCount(2) <= NumSurrDisagree).and.(AlleleCount(1) > AlleleCount(2))) iAllele = 0
    if ((AlleleCount(1) <= NumSurrDisagree).and.(AlleleCount(2) > AlleleCount(1))) iAllele = 1
  endif
  deallocate(ErdosNowVec)
  deallocate(ErdosNextVec)
  return

end function IterAllele

subroutine newErdos(surrogates, threshold, c)
  use SurrogateDefinition
  use CoreSubsetDefinition
  use Parameters, only: useSurrsN, consistent, numsurrdisagree
  implicit none

  type(SurrDef), intent(in) :: surrogates
  type(CoreSubSet) :: c
  integer, intent(in) :: threshold
  
  integer(kind=1), dimension(:,:), pointer :: genos
  integer, dimension(:,:), allocatable :: surrList
  integer, dimension(:), allocatable :: nSurrList
  integer :: nAnisG, nSnp
  integer :: i, j, k, side
  logical, dimension(:), allocatable :: visited
  integer, dimension(:), allocatable :: toVisit
  integer, dimension(:), allocatable :: depth
  integer :: toVisitPos, visiting, next
  integer :: found
  integer, dimension(2) :: alleleCount
  integer :: total, counter
  integer, dimension(:), allocatable :: surrAveDiff
  integer(kind=1) :: iAllele
  integer :: current
  integer :: highestErdos
  
  nAnisG = c%getNAnisG()
  nSnp = c%getNCoreSnp()
  genos => c%getCoreGenos()
  
  
  allocate(surrList(nAnisG,nAnisG))
  allocate(nSurrList(nAnisG))
  allocate(visited(nAnisG))
  allocate(toVisit(nAnisG))
  allocate(depth(nAnisG))
  allocate(surrAveDiff(nAnisG))
  
  nSurrList = 0
  
  do i = 1, nAnisG
    do j = 1, nAnisG
      if (surrogates%numOppose(i,j) <= threshold) then
	if (i /= j) then
	  nSurrList(i) = nSurrList(i) + 1
	  surrList(i,nSurrList(i)) = j
	end if
      end if
    end do
  end do
  
  do i = 1, nAnisG
    total = 0
    counter = 0
    do j = 1, nAnisG
      if (surrogates%numOppose(i, j) > threshold) then
	total = total+surrogates%numOppose(i, j)
	counter = counter + 1
      endif
    end do
    SurrAveDiff(i) = int(total/counter)
  end do
  
  do side = 1, 2
    print*, " "
    select case (side)
      case (1)
	print*, " Phasing genotyped individuals for Paternal allele"
      case (2)
	print*, " Phasing genotyped individuals for Maternal allele"
    end select
    highestErdos = 0
    do i = 1, nAnisG
      if (mod(i, 400) == 0) then
	print*, "   Phasing done for genotyped individual --- ", i
      end if
      do j = 1, nSnp
      	visited = .false.
	toVisit = 0
	toVisit(1) = i
	Visited(i) = .true.
	toVisitPos = 1
	found = 0
	depth = 0
	allelecount = 0
	visiting = 1
	
	!!! FUDGES !!!
	if (genos(i,j) == 0) then
	  AlleleCount(1) = AlleleCount(1) + 1
	  found = found + 1
	end if
	if (genos(i,j) == 2) then
	  AlleleCount(2) = AlleleCount(2) + 1
	  found = found + 1
	end if
	!!! END FUDGES !!!
	
	do while ((toVisit(visiting) /= 0) .and.(found < useSurrsN))
	  current = toVisit(visiting)
	  k = 1
	  do while ((k <= nSurrList(current)) .and. (found < useSurrsN))
	    next = surrList(current,k)
	    if (goodToVisit(next, depth, side, surrogates, visited,  &
	      toVisit, visiting, threshold, SurrAveDiff)) then
	      select case (Genos(next, j))
		case (0)
		  if (mod(depth(Visiting) + 1, 2) == 0) then
		    AlleleCount(2) = AlleleCount(2) + 1
		  end if
		  if (mod(depth(Visiting) + 1, 2) /= 0) then
		    AlleleCount(1) = AlleleCount(1) + 1
		  end if
		  found = found + 1
		case (2)
		  if (mod(depth(Visiting) + 1, 2) == 0) then
		    AlleleCount(1) = AlleleCount(1) + 1
		  end if
		  if (mod(depth(Visiting) + 1, 2) /= 0) then
		    AlleleCount(2) = AlleleCount(2) + 1
		  end if
		  found = found + 1
		case default		      
		  toVisitPos = toVisitPos + 1
		  toVisit(toVisitPos) = next
		  depth(toVisitPos) = depth(visiting) + 1
	      end select
	    end if
	    visited(next) = .true.
	    k = k + 1
	  end do
	  visiting = visiting + 1
	  if (depth(visiting) /= depth(visiting-1)) then
	    call phaseSort(toVisit, visiting, toVisitPos)
	  end if
	end do
	
	iAllele = 9
	if (sum(AlleleCount(:)) < UseSurrsN) then
    	  if ((AlleleCount(2) <= NumSurrDisagree).and.(AlleleCount(1) > AlleleCount(2))) iAllele = 0
	  if ((AlleleCount(1) <= NumSurrDisagree).and.(AlleleCount(2) > AlleleCount(1))) iAllele = 1
	else
	  if ((AlleleCount(1) > 0).and.(AlleleCount(2) <= NumSurrDisagree)) then
	    iAllele = 0
	  elseif ((AlleleCount(2) > 0).and.(AlleleCount(1) <= NumSurrDisagree)) then
	    iAllele = 1
	  end if
	endif
	call c%setPhase(i, j, side, iAllele)
	! +1 as we use the genotypes at one depth lower than we're currently visiting (slightly odd but maintains consistency
	! with old version).  See discussion about odd depth first search above (when it's written!)
	if (depth(visiting - 1) + 1 > highestErdos) then
	  highestErdos = depth(visiting - 1) + 1
	end if
      end do
    end do
    print*, " "
    select case (side)
      case (1)
	print*, " ", highestErdos, " was the highest Erdos used on Paternal Side"
      case (2)
	print*, " ", highestErdos, " was the highest Erdos used on Maternal Side"
    end select
  end do
	
end subroutine newErdos

function goodToVisit(next, depths, sideon, surrogates, visited, &
      tovisit, visiting, threshold, surravediff) result(good)
  use SurrogateDefinition
  
  integer, intent(in) :: next, sideon, threshold, visiting
  type(SurrDef), intent(in) :: surrogates
  logical, dimension(:), intent(in) :: visited
  integer, dimension(:), intent(in) :: surravediff
  integer, dimension(:), intent(in) :: tovisit, depths
  logical :: good
  
  integer :: depth, current
  
  integer :: i
  
  depth = depths(visiting)
  current = tovisit(visiting)
  
  good = .not. visited(next)
  
  if (good .and. (depth == 0)) then
    good = (surrogates%partition(current, next) == sideon)
  end if
  
  if (good) then
    do i = 1, visiting - 1
      if ( &
	(surrogates%numOppose(next, tovisit(i)) > threshold) .and. &
	(surrogates%numOppose(next, tovisit(i)) <= SurrAveDiff(next))) then
        good = .false.
      end if
    end do
  end if
  
end function goodToVisit

subroutine phaseSort(array, p1, p2)
  integer, dimension(:), intent(inout) :: array
  integer, intent(in) :: p1, p2
  
  integer(8) :: length, size
  
  length = p2-p1+1
  size = 4
  
  call qsort(array(p1:p2), length, size, cmp_function)
end subroutine phaseSort

function cmp_function(a1, a2) result (cmp)
    integer(2), intent(in) :: a1, a2
    integer(2) :: cmp
    
    cmp = a1-a2
end function

end module Phasing
