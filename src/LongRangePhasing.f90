module LongRangePhasing
  use Constants
  implicit none
  
contains
  
  subroutine Erdos(surrogates, c, numsurrdisagree, useSurrsN, printProgress)
    use SurrogateDefinition
    use CoreSubsetDefinition
    use GenotypeModule
    
    type(Surrogate), intent(in) :: surrogates
    type(CoreSubSet) :: c
    integer, intent(in) :: numsurrdisagree, useSurrsN
    logical, intent(in) :: printProgress

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
    integer(kind=1) :: iAllele
    integer :: current
    integer :: highestErdos

    nAnisG = c%getNAnisG()
    nSnp = c%getNCoreSnp()
    genos => c%getCoreGenos()


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
	(surrogates%numIncommon(i,j) >= surrogates%incommonThreshold)) then
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
!      if (ItterateType .eq. "Off") then
      if (printProgress) then
	print*, " "
	select case (side)
	  case (1)
	    print*, " Phasing genotyped individuals for Paternal allele"
	  case (2)
	    print*, " Phasing genotyped individuals for Maternal allele"
	end select
      end if
      highestErdos = 0
      do i = 1, nAnisG
	if ((mod(i, 400) == 0) .and. printProgress) then
	  print*, "   Phasing done for genotyped individual --- ", i
	end if
	do j = 1, nSnp
	  if (c%getPhase(i, j, side) == MissingPhaseCode) then
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

	    !!! FUDGES !!!
	    if (genos(i)%getGenotype(j) == 0) then
	      AlleleCount(1) = AlleleCount(1) + 1
	      found = found + 1
	    end if
	    if (genos(i)%getGenotype(j) == 2) then
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
		  toVisit, visiting, SurrAveDiff)) then
		  select case (Genos(next)%getGenotype(j))
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
	      !Bit of a hack in case all animals are surrogates
	      if (visiting <= nAnisG) then
		if (depth(visiting) /= depth(visiting-1)) then
		  call phaseSort(toVisit, visiting, toVisitPos)
		end if
	      end if
	    end do

	    iAllele = MissingPhaseCode
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
	  end if
	end do
      end do
      if (printProgress) then
	print*, " "
	select case (side)
	  case (1)
	    print*, " ", highestErdos, " was the highest Erdos used on Paternal Side"
	  case (2)
	    print*, " ", highestErdos, " was the highest Erdos used on Maternal Side"
	end select
      end if
    end do
    
    deallocate(surrList)
    deallocate(nSurrList)
    deallocate(visited)
    deallocate(toVisit)
    deallocate(depth)
    deallocate(surrAveDiff)
    deallocate(genos)

  end subroutine Erdos

  function goodToVisit(next, depths, sideon, surrogates, visited, &
	tovisit, visiting, surravediff) result(good)
    use SurrogateDefinition

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
	  (surrogates%numIncommon(next,toVisit(i)) < surrogates%incommonThreshold)) then
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

  subroutine CheckCompatHapGeno(c, PercGenoHaploDisagree, printProgress)
    use Constants
    use CoreSubsetDefinition
    use GenotypeModule
    
    class(CoreSubset) :: c
    double precision, intent(in) :: PercGenoHaploDisagree
    logical, intent(in) :: printProgress

    type(Genotype), dimension(:), pointer :: genos

    integer :: i, j, CountError, ErrorAllow, counterMissing, nAnisG, nCoreSnp

    nAnisG  = c%getNAnisG()
    nCoreSnp = c%getNCoreSnp()

    genos => c%getCoreGenos()

    ErrorAllow = int(PercGenoHaploDisagree * nCoreSnp)

    do i = 1, nAnisG
      CountError = 0
      counterMissing = 0
      do j = 1, nCoreSnp
	if ((c%getPhase(i, j, 1) /= MissingPhaseCode).and.(c%getPhase(i, j, 2) /= MissingPhaseCode)) then
	  counterMissing = counterMissing + 1
	  if ((Genos(i)%getGenotype(j) /= MissingGenotypeCode).and.(c%getPhaseGeno(i,j)  /= Genos(i)%getGenotype(j))) CountError = CountError + 1
	end if
      end do
      ErrorAllow = int(PercGenoHaploDisagree * counterMissing)
      if (CountError >= ErrorAllow) then
	do j = 1, nCoreSnp
	  if (Genos(i)%getGenotype(j) /= MissingGenotypeCode) then
	    if ((c%getPhase(i, j, 1) /= MissingPhaseCode).and.(c%getPhase(i, j, 2) /= MissingPhaseCode).and. &
	      (c%getPhaseGeno(i, j) /= Genos(i)%getGenotype(j))) then
	      if (Genos(i)%getGenotype(j) == 1) call c%setPhase(i, j, 2, MissingPhaseCode)
	      if (Genos(i)%getGenotype(j) == MissingGenotypeCode) call c%setPhase(i, j, 2, MissingPhaseCode)
	      if (Genos(i)%getGenotype(j) == 0) then
		call c%setPhase(i, j, 1, 0)
		call c%setPhase(i, j, 2, 0)
	      end if
	      if (Genos(i)%getGenotype(j) == 2) then
		call c%setPhase(i, j, 1, 1)
		call c%setPhase(i, j, 2, 1)
	      end if
	    endif
	  endif
	enddo
      endif
    end do

    
    !! THIS REALLY SHOULDN'T BE HERE!!
    if (printProgress) then
      print*, " "
      write (*, '(a3,f6.2,a45)') "  ", c%getYield(1), "% was the Paternal allele yield for this core"
      write (*, '(a3,f6.2,a45)') "  ", c%getYield(2), "% was the Maternal allele yield for this core"
    end if
    
    deallocate(genos)

  end subroutine CheckCompatHapGeno

end module LongRangePhasing
