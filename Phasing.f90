module Phasing
  contains
  
  subroutine Erdos(surrogates, threshold, genos, phase, startCoreSnp, endCoreSnp)
!  use Global
  use SurrogateDefinition
  use global, only: useSurrsN
  implicit none

  type(SurrDef), intent(in) :: surrogates
  integer(kind=1), dimension(:,:), intent(in) :: genos
  integer, intent(in) :: threshold
  integer(kind=1), dimension(:,:,:), intent(inout) :: phase
  integer, intent(in) :: startCoreSnp, EndCoreSnp
    
  integer(kind = 1), allocatable, dimension (:) :: Visited
  integer, allocatable, dimension (:) :: SurrAveDiff
  integer :: AlleleCount(2)
  integer :: ErdosNumber, HighestErdos
  
  integer :: nAnisG
  
  integer :: i, j
  !integer :: counter, IterAllele, SizeCore
  integer :: counter, SizeCore
  double precision :: value

  integer :: GetnSnpErrorThreshAnims

  nAnisG = size(genos,1)
  
  allocate(Visited(nAnisG))
  allocate(SurrAveDiff(nAnisG))

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

  Phase(:, StartCoreSnp:EndCoreSnp,:) = 9
  HighestErdos = 1
  print*, " "
  print*, " Phasing genotyped individuals for Paternal allele"

  do i = 1, nAnisG
    if (mod(i, 400) == 0) print*, "   Phasing done for genotyped individual --- ", i
    do j = StartCoreSnp, EndCoreSnp
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
	  Phase(i, j, 1) = 0
	elseif (Genos(i, j) == 2) then
	  Phase(i, j, 1) = 1
	endif
      else
	Phase(i, j, 1) = IterAllele(i, j, 1, surrogates, threshold, visited, allelecount, erdosnumber, genos, surravediff)
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
    do j = StartCoreSnp, EndCoreSnp
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
	  Phase(i, j, 2) = 0
	elseif (Genos(i, j) == 2) then
	  Phase(i, j, 2) = 1
	endif
      else
	Phase(i, j, 2) = IterAllele(i, j, 2, surrogates, threshold, visited, allelecount, erdosnumber, genos, surravediff)
	if (ErdosNumber > HighestErdos) HighestErdos = ErdosNumber
      end if
    end do
  end do

  print*, " "
  print*, " ", HighestErdos, " was the highest Erdos used on Maternal Side"

  deallocate(Visited)
  deallocate(SurrAveDiff)

end subroutine Erdos

!#################################################################################################################################################################

function IterAllele(animal, snp, SideOn, surrogates, threshold, visited, allelecount, erdosnumber, genos, surravediff) result (iAllele)
  !use Global
  use SurrogateDefinition
  use Global, only: MissingGenotypeCode, UseSurrsN, NumSurrDisagree
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

  integer :: GetnSnpErrorThreshAnims
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
      ! if ((Surrogates(animal,i,1)>GetnSnpErrorThreshAnims(animal,i)).and.(Surrogates(animal,i,1)<=(GetnSnpErrorThreshAnims(animal,i)+15))) then
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

  99 ErdosNumber = ErdosNumber + 1
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
	    if (mod(ErdosNumber, 2) == 0) AlleleCount(2) = AlleleCount(2) + 1
	    if (mod(ErdosNumber, 2) /= 0) AlleleCount(1) = AlleleCount(1) + 1
	  end if
	  if (genos(j, snp) == 2) then
	    if (mod(ErdosNumber, 2) == 0) AlleleCount(1) = AlleleCount(1) + 1
	    if (mod(ErdosNumber, 2) /= 0) AlleleCount(2) = AlleleCount(2) + 1
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
  if (sum(ErdosNextVec(:)) /= 0) goto 99
  iAllele = 9
  if (sum(AlleleCount(:)) < UseSurrsN) then
    if ((AlleleCount(2) <= NumSurrDisagree).and.(AlleleCount(1) > AlleleCount(2))) iAllele = 0
    if ((AlleleCount(1) <= NumSurrDisagree).and.(AlleleCount(2) > AlleleCount(1))) iAllele = 1
  endif
  deallocate(ErdosNowVec)
  deallocate(ErdosNextVec)
  return

end function IterAllele
end module Phasing
