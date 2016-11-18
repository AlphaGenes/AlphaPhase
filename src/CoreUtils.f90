module CoreUtils
  
contains
  function calculateCores(nSnp, Jump, offset, consistent) result(CoreIndex)
    implicit none

    integer, intent(in) :: nSnp, Jump
    logical, intent(in) :: offset, consistent
    integer, dimension(:,:), pointer :: CoreIndex

    double precision :: corelength
    integer :: i, nCores, left

    if (consistent) then
      if (.not. Offset) then
	nCores = int(nSnp)/Jump
	allocate(CoreIndex(nCores, 2))

	CoreIndex(1, 1) = 1
	CoreIndex(1, 2) = 1 + Jump - 1
	do i = 2, nCores
	  CoreIndex(i, 1) = CoreIndex(i - 1, 1) + Jump
	  CoreIndex(i, 2) = CoreIndex(i - 1, 2) + Jump
	end do
	CoreIndex(nCores, 2) = nSnp
      endif

      if (Offset) then
	nCores = (int(nSnp)/Jump) + 1
	allocate(CoreIndex(nCores, 2))

	CoreIndex(1, 1) = 1
	CoreIndex(1, 2) = int(Jump/2)
	do i = 2, nCores
	  CoreIndex(i, 1) = CoreIndex(i - 1, 2) + 1
	  CoreIndex(i, 2) = CoreIndex(i - 1, 2) + Jump
	end do
	CoreIndex(nCores, 2) = nSnp
      endif
    else
      nCores = nSnp / Jump
      corelength = nSnp / nCores
      left = nSnp - nCores * corelength

      if (.not. Offset) then
	allocate(CoreIndex(nCores, 2))
	CoreIndex(1, 1) = 1
	if (left /= 0) then
	  CoreIndex(1, 2) = 1 + corelength
	else
	  CoreIndex(1, 2) = corelength
	end if
      else
	nCores = nCores + 1
	allocate(CoreIndex(nCores, 2))
	CoreIndex(1, 1) = 1
	if (left /= 0) then
	  CoreIndex(1, 2) = 1 + floor(dble(corelength) / 2.0)
	else
	  CoreIndex(1, 2) = floor(dble(corelength) / 2.0)
	end if
      end if

      do i = 2, nCores
	CoreIndex(i,1) = CoreIndex(i - 1, 2) + 1
	if (i <= left) then
	  CoreIndex(i, 2) = CoreIndex(i - 1, 2) + corelength + 1
	else
	  CoreIndex(i, 2) = CoreIndex(i - 1, 2) + corelength
	end if
      end do

      if (Offset) then
	CoreIndex(nCores,2) = nSnp
      end if
    endif
  end function CalculateCores

  function CalculateTails(CoreIndex, nSnp, Jump, CoreAndTailLength, offset, consistent) result(TailIndex)
    implicit none

    integer, dimension(:,:), intent(in) :: CoreIndex
    integer, intent(in) :: nSnp, Jump, CoreAndTailLength
    logical, intent(in) :: offset, consistent
    integer, dimension(:,:), pointer :: TailIndex

    integer :: resid
    integer :: ltail, rtail, nCores
    integer :: i

    nCores = size(CoreIndex,1)
    allocate(TailIndex(nCores,2))

    if (consistent) then
      if (.not. Offset) then
	allocate(TailIndex(nCores, 2))

	resid = int((CoreAndTailLength - Jump)/2)
	TailIndex(1, 1) = 1
	TailIndex(1, 2) = 1 + CoreAndTailLength - 1
	do i = 2, nCores
	  TailIndex(i, 1) = CoreIndex(i, 1) - resid
	  TailIndex(i, 2) = CoreIndex(i, 2) + resid
	  if (TailIndex(i, 1) < 1) TailIndex(i, 1) = 1
	  if (TailIndex(i, 2) > nSnp) TailIndex(i, 2) = nSnp
	end do
      endif

      if (Offset) then
	resid = int((CoreAndTailLength - Jump)/2)

	allocate(TailIndex(nCores, 2))

	TailIndex(1, 1) = 1
	TailIndex(1, 2) = nSnp
	do i = 2, nCores
	  TailIndex(i, 1) = CoreIndex(i, 1) - resid
	  TailIndex(i, 2) = CoreIndex(i, 2) + resid
	  if (TailIndex(i, 1) < 1) TailIndex(i, 1) = 1
	  if (TailIndex(i, 2) > nSnp) TailIndex(i, 2) = nSnp
	end do
	TailIndex(nCores, 1) = 1   
	TailIndex(nCores, 2) = nSnp
      endif
    else
      ltail = floor(dble(CoreAndTailLength - Jump) / 2.0)
      rtail = ceiling(dble(CoreAndTailLength - Jump) / 2.0)


      do i = 1, nCores
	TailIndex(i,1) = max(1,CoreIndex(i,1) - ltail)
	TailIndex(i,2) = min(nSnp,CoreIndex(i,2) + rtail)
      end do
    endif
  end function CalculateTails
  
  function MinorAlleleFrequency(genos) result(maf)
    use Constants
    
    integer(kind=1), dimension(:,:), intent(in) :: genos
    double precision, dimension(:), pointer :: maf
    
    integer :: i, j
    integer :: sum, count
    
    allocate(maf(size(genos,2)))
    
    do i = 12, 14 !size(genos,2)
      sum = 0
      count = 0
      
      do j = 1, size(genos,1)
	if (genos(j,i) /= MissingGenotypeCode) then
	  sum = sum + genos(j,i)
	  count = count + 1
	end if
      end do
      
      maf(i) = float(sum) / (2*float(count))
      !print '(i6, 10x, f7.5)', sum, maf(i)
      !print '(f7.5)', maf(i)
    end do
    
  end function MinorAlleleFrequency
  
  function NumNeeded(maf, threshold) result(num)
    double precision, dimension(:), intent(in) :: maf
    double precision, intent(in) :: threshold
    integer, dimension(size(maf)) :: num
    integer :: i, j
    
    double precision, dimension(size(maf)) :: cumulative
    double precision :: poppose, m, t
    
    num = size(maf) + 1
    cumulative = 1.0
    t = 1 - threshold
    
    do i = 1, size(maf)
      m = maf(i)
      poppose = 2*m*m*(1-m)*(1-m)
      
      do j = 1, i
	if (num(j) == size(maf) + 1) then
	  cumulative(j) = cumulative(j) * (1-poppose)
	  if (cumulative(j) < t) then
	    num(j) = i - j + 1
	  end if
	end if
      end do
    end do

    print *, minval(num, num /= size(maf)+1), sum(num, num /= size(maf) + 1) / count(num /= size(maf) + 1),  &
      maxval(num, num /= size(maf) + 1)
  end function NumNeeded
  
  function TransitionP(genos) result(same)
    use Constants
    
    integer(kind=1), dimension(:,:), intent(in) :: genos
    double precision, dimension(:,:), pointer :: same
    
    integer :: count0, sum0, count1, sum1
    integer :: i, j
    
    allocate(same(size(genos)-1,0:1))
    
    do i = 12, 14 !size(genos,2) - 1
      count0 = 0
      sum0 = 0
      count1 = 0
      sum1 = 0
      
      do j = 1, size(genos,1)
	if (genos(j,i) == 0) then
	  if (genos(j,i + 1) == 0) then
	    sum0 = sum0 + 2
	  end if
	  if (genos(j,i + 1) == 1) then
	    sum0 = sum0 + 0
	  end if
          count0 = count0 + 2
	end if
	if (genos(j,i) == 2) then
	  if (genos(j,i + 1) == 2) then
	    sum1 = sum1 + 2
	  end if
	  if (genos(j,i + 1) == 1) then
	    sum1 = sum1 + 1
	  end if
	  count1 = count1 + 2
	end if
	
	if (genos(j,i) == 1) then
	  if (genos(j,i + 1) == 0) then
	    sum0 = sum0 + 1
	    count0 = count0 + 1
	    count1 = count1 + 1
	  end if
	  if (genos(j,i + 1) == 2) then
	    sum1 = sum1 + 1
	    count0 = count0 + 1
	    count1 = count1 + 1
	  end if
	end if
      end do
      
      same(i,0) = float(sum0) / float(count0)
      same(i,1) = float(sum1) / float(count1)
    end do
  end function TransitionP
      
   
end module CoreUtils

