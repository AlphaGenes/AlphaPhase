module LongRangePhasing
  use ConstantModule
  implicit none
  
  abstract interface
    function gtvTest(surrogates, surrAveDiff, toVisit, next) result (gtv)
      use SurrogateModule

      type(Surrogate), intent(in) :: surrogates
      integer, dimension(:), intent(in) :: surravediff
      integer, intent(in) :: toVisit, next
      
      logical :: gtv
    end function gtvTest
  end interface
contains
  
  subroutine ErdosWithOverlap(surrogates, c, numsurrdisagree, useSurrsN)
    use SurrogateModule
    use CoreSubsetModule
    use GenotypeModule
    use HaplotypeModule
    
    type(Surrogate), intent(in) :: surrogates
    type(CoreSubSet) :: c
    integer, intent(in) :: numsurrdisagree, useSurrsN

    type(Genotype), dimension(:), pointer :: genos
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
    integer :: current
    integer :: highestErdos
    type(Haplotype), pointer :: hap
   
    nAnisG = c%getNAnisGCoreSubset()
    nSnp = c%getNCoreSnpCoreSubset()
    genos => c%getCoreGenosCoreSubset()


    allocate(surrList(nAnisG,nAnisG))
    allocate(nSurrList(nAnisG))
    allocate(visited(nAnisG))
    allocate(toVisit(nAnisG + 1))
    allocate(depth(nAnisG))
    allocate(surrAveDiff(nAnisG))
    
    nSurrList = 0

    do i = 1, nAnisG
      do j = 1, nAnisG
	if ((surrogates%numOppose(i,j) <= surrogates%threshold) .and. &
	surrogates%enoughIncommon(i,j)) then
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
	if ((surrogates%numOppose(i, j) > surrogates%threshold) .and. &
	surrogates%enoughInCommon(i,j)) then
	  total = total+surrogates%numOppose(i, j)
	  counter = counter + 1
	endif
      end do
      if (counter /= 0) then
	SurrAveDiff(i) = int(total/counter)
      else
	! All animals are surrogates of each other.  Can't see this happening in anything other than simulated data
	! Should probably think about more and check with John and possibly others.
	SurrAveDiff(i) = 0
      end if
    end do

    do side = 1, 2
      highestErdos = 0
      do i = 1, nAnisG
	hap => c%getHaplotypeCoreSubset(i, side)
	do j = 1, nSnp
	  if (hap%isMissing(j)) then
	    visited = .false.
	    toVisit = 0
	    toVisit(1) = i
	    ! Bit of a hack in case all animals are surrogates
	    toVisit(nAnisG + 1) = 0
	    Visited(i) = .true.
	    toVisitPos = 1
	    found = 0
	    depth = 0
	    allelecount = 0
	    visiting = 1

	    if (genos(i)%isZero(j)) then
	      AlleleCount(1) = AlleleCount(1) + 1
	      found = found + 1
	    end if
	    if (genos(i)%isTwo(j)) then
	      AlleleCount(2) = AlleleCount(2) + 1
	      found = found + 1
	    end if

	    do while ((toVisit(visiting) /= 0) .and.(found < useSurrsN))
	      current = toVisit(visiting)
	      k = 1
	      do while ((k <= nSurrList(current)) .and. (found < useSurrsN))
		next = surrList(current,k)
		if (goodToVisitWithOverlap(next, depth, side, surrogates, visited,  &
		  toVisit, visiting, SurrAveDiff)) then
		  if (Genos(next)%isHomo(j)) then
		    if (Genos(next)%isZero(j)) then
		      if (mod(depth(Visiting) + 1, 2) == 0) then
			AlleleCount(2) = AlleleCount(2) + 1
		      end if
		      if (mod(depth(Visiting) + 1, 2) /= 0) then
			AlleleCount(1) = AlleleCount(1) + 1
		      end if
		      found = found + 1
		    else
		      if (mod(depth(Visiting) + 1, 2) == 0) then
			AlleleCount(1) = AlleleCount(1) + 1
		      end if
		      if (mod(depth(Visiting) + 1, 2) /= 0) then
			AlleleCount(2) = AlleleCount(2) + 1
		      end if
		      found = found + 1
		    end if
		  else
		    if (.not. Genos(next)%isMissing(j)) then
		      toVisitPos = toVisitPos + 1
		      toVisit(toVisitPos) = next
		      depth(toVisitPos) = depth(visiting) + 1
		    end if
		  end if
		end if
		visited(next) = .true.
		k = k + 1
	      end do
	      visiting = visiting + 1
	      !Bit of a hack in case all animals are surrogates
	      if (visiting <= nAnisG) then
		if (depth(visiting) /= depth(visiting-1)) then
		  call phaseSort(toVisit, visiting, toVisitPos)
		end if
	      end if
	    end do

	    if (sum(AlleleCount(:)) < UseSurrsN) then
	      if ((AlleleCount(2) <= NumSurrDisagree).and.(AlleleCount(1) > AlleleCount(2))) then
		call hap%setZero(j)
	      end if
	      if ((AlleleCount(1) <= NumSurrDisagree).and.(AlleleCount(2) > AlleleCount(1))) then
		call hap%setOne(j)
	      end if
	    else
	      if ((AlleleCount(1) > 0).and.(AlleleCount(2) <= NumSurrDisagree)) then
		call hap%setZero(j)
	      elseif ((AlleleCount(2) > 0).and.(AlleleCount(1) <= NumSurrDisagree)) then
		call hap%setOne(j)
	      end if
	    endif
	  end if
	end do
      end do
    end do
    
    deallocate(surrList)
    deallocate(nSurrList)
    deallocate(visited)
    deallocate(toVisit)
    deallocate(depth)
    deallocate(surrAveDiff)
    deallocate(genos)

  end subroutine ErdosWithOverlap
  
  subroutine ErdosWithoutOverlap(surrogates, c, numsurrdisagree, useSurrsN)
    use SurrogateModule
    use CoreSubsetModule
    use GenotypeModule
    use HaplotypeModule
    
    type(Surrogate), intent(in) :: surrogates
    type(CoreSubSet) :: c
    integer, intent(in) :: numsurrdisagree, useSurrsN

    type(Genotype), dimension(:), pointer :: genos
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
    integer :: current
    integer :: highestErdos
    type(Haplotype), pointer :: hap
   
    nAnisG = c%getNAnisGCoreSubset()
    nSnp = c%getNCoreSnpCoreSubset()
    genos => c%getCoreGenosCoreSubset()


    allocate(surrList(nAnisG,nAnisG))
    allocate(nSurrList(nAnisG))
    allocate(visited(nAnisG))
    allocate(toVisit(nAnisG + 1))
    allocate(depth(nAnisG))
    allocate(surrAveDiff(nAnisG))
    
    nSurrList = 0

    do i = 1, nAnisG
      do j = 1, nAnisG
	if (surrogates%numOppose(i,j) <= surrogates%threshold) then
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
	if (surrogates%numOppose(i, j) > surrogates%threshold) then
	  total = total+surrogates%numOppose(i, j)
	  counter = counter + 1
	endif
      end do
      if (counter /= 0) then
	SurrAveDiff(i) = int(total/counter)
      else
	! All animals are surrogates of each other.  Can't see this happening in anything other than simulated data
	! Should probably think about more and check with John and possibly others.
	SurrAveDiff(i) = 0
      end if
    end do

    do side = 1, 2
      highestErdos = 0
      do i = 1, nAnisG
	hap => c%getHaplotypeCoreSubset(i, side)
	do j = 1, nSnp
	  if (hap%isMissing(j)) then
	    visited = .false.
	    toVisit = 0
	    toVisit(1) = i
	    ! Bit of a hack in case all animals are surrogates
	    toVisit(nAnisG + 1) = 0
	    Visited(i) = .true.
	    toVisitPos = 1
	    found = 0
	    depth = 0
	    allelecount = 0
	    visiting = 1

	    if (genos(i)%isZero(j)) then
	      AlleleCount(1) = AlleleCount(1) + 1
	      found = found + 1
	    end if
	    if (genos(i)%isTwo(j)) then
	      AlleleCount(2) = AlleleCount(2) + 1
	      found = found + 1
	    end if

	    do while ((toVisit(visiting) /= 0) .and.(found < useSurrsN))
	      current = toVisit(visiting)
	      k = 1
	      do while ((k <= nSurrList(current)) .and. (found < useSurrsN))
		next = surrList(current,k)
		if (goodToVisitWithoutOverlap(next, depth, side, surrogates, visited,  &
		  toVisit, visiting, SurrAveDiff)) then
		  if (Genos(next)%isHomo(j)) then
		    if (Genos(next)%isZero(j)) then
		      if (mod(depth(Visiting) + 1, 2) == 0) then
			AlleleCount(2) = AlleleCount(2) + 1
		      end if
		      if (mod(depth(Visiting) + 1, 2) /= 0) then
			AlleleCount(1) = AlleleCount(1) + 1
		      end if
		      found = found + 1
		    else
		      if (mod(depth(Visiting) + 1, 2) == 0) then
			AlleleCount(1) = AlleleCount(1) + 1
		      end if
		      if (mod(depth(Visiting) + 1, 2) /= 0) then
			AlleleCount(2) = AlleleCount(2) + 1
		      end if
		      found = found + 1
		    end if
		  else
		    if (.not. Genos(next)%isMissing(j)) then
		      toVisitPos = toVisitPos + 1
		      toVisit(toVisitPos) = next
		      depth(toVisitPos) = depth(visiting) + 1
		    end if
		  end if
		end if
		visited(next) = .true.
		k = k + 1
	      end do
	      visiting = visiting + 1
	      !Bit of a hack in case all animals are surrogates
	      if (visiting <= nAnisG) then
		if (depth(visiting) /= depth(visiting-1)) then
		  call phaseSort(toVisit, visiting, toVisitPos)
		end if
	      end if
	    end do

	    if (sum(AlleleCount(:)) < UseSurrsN) then
	      if ((AlleleCount(2) <= NumSurrDisagree).and.(AlleleCount(1) > AlleleCount(2))) then
		call hap%setZero(j)
	      end if
	      if ((AlleleCount(1) <= NumSurrDisagree).and.(AlleleCount(2) > AlleleCount(1))) then
		call hap%setOne(j)
	      end if
	    else
	      if ((AlleleCount(1) > 0).and.(AlleleCount(2) <= NumSurrDisagree)) then
		call hap%setZero(j)
	      elseif ((AlleleCount(2) > 0).and.(AlleleCount(1) <= NumSurrDisagree)) then
		call hap%setOne(j)
	      end if
	    endif
	  end if
	end do
      end do
    end do
    
    deallocate(surrList)
    deallocate(nSurrList)
    deallocate(visited)
    deallocate(toVisit)
    deallocate(depth)
    deallocate(surrAveDiff)
    deallocate(genos)

  end subroutine ErdosWithoutOverlap

  function goodToVisitWithOverlap(next, depths, sideon, surrogates, visited, &
    	tovisit, visiting, surravediff) result(good)
    use SurrogateModule

    integer, intent(in) :: next, sideon, visiting
    type(Surrogate), intent(in) :: surrogates
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
	  ((surrogates%numOppose(next, tovisit(i)) > surrogates%threshold) .and. &
	  (surrogates%numOppose(next, tovisit(i)) <= SurrAveDiff(next))) .or. &
	  .not. surrogates%enoughIncommon(next,toVisit(i))) then
	  good = .false.
	end if
      end do
    end if

  end function goodToVisitWithOverlap
  
  function goodToVisitWithoutOverlap(next, depths, sideon, surrogates, visited, &
    	tovisit, visiting, surravediff) result(good)
    use SurrogateModule

    integer, intent(in) :: next, sideon, visiting
    type(Surrogate), intent(in) :: surrogates
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
	  ((surrogates%numOppose(next, tovisit(i)) > surrogates%threshold) .and. &
	  (surrogates%numOppose(next, tovisit(i)) <= SurrAveDiff(next)))) then
	  good = .false.
	end if
      end do
    end if

  end function goodToVisitWithoutOverlap
  
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

  subroutine CheckCompatHapGeno(c, PercGenoHaploDisagree)
    use CoreSubsetModule
    use GenotypeModule
    use HaplotypeModule
    use BitUtilities
    
    class(CoreSubset) :: c
    double precision, intent(in) :: PercGenoHaploDisagree

    type(Genotype), pointer :: genos
    type(Haplotype), pointer :: hap1, hap2

    integer :: i, CountError, ErrorAllow, counterMissing, nAnisG, nCoreSnp
    integer(kind=8), dimension(:), allocatable :: error
    
    nAnisG  = c%getNAnisGCoreSubset()
    nCoreSnp = c%getNCoreSnpCoreSubset()

    ErrorAllow = int(PercGenoHaploDisagree * nCoreSnp)

    do i = 1, nAnisG
      CountError = 0
      counterMissing = 0
      hap1 => c%getHaplotypeCoreSubset(i,1)
      hap2 => c%getHaplotypeCoreSubset(i,2)
      genos => c%getSingleCoreGenos(i)
      error = genos%getErrors(hap1,hap2)
      countError = bitCount(error)
      counterMissing = hap1%numberBothNotMissing(hap2)
      ErrorAllow = int(PercGenoHaploDisagree * counterMissing)
      if (CountError >= ErrorAllow) then
	call genos%setHaplotypeFromGenotypeIfError(hap1,error)
	call genos%setHaplotypeFromGenotypeIfError(hap2,error)
      endif
    end do
  end subroutine CheckCompatHapGeno

end module LongRangePhasing
