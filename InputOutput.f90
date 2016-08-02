module InputOutput
  
contains

  subroutine WriteOutResults(phase, allHapAnis, coreIndex, p)
    use Constants
    use PedigreeDefinition
    implicit none

    integer(kind=1), dimension(:,:,:), intent(in) :: phase
    integer, dimension(:,:,:), intent(in) :: allHapAnis
    integer, dimension(:,:) :: coreIndex
    type(Pedigree) :: p

    integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp, nCores
    integer, allocatable, dimension(:) :: WorkOut
    double precision, allocatable, dimension(:) :: CoreCount

    nAnisG = size(phase,1)
    nSnp = size(phase,2)
    nCores = size(allHapAnis,3)

    allocate(WorkOut(nCores * 2))
    allocate(CoreCount(nCores * 2))

    if (WindowsLinux == 1) then
      open (unit = 15, file = ".\PhasingResults\FinalPhase.txt", status = "unknown")
      open (unit = 25, file = ".\PhasingResults\CoreIndex.txt", status = "unknown")
      open (unit = 28, file = ".\PhasingResults\SnpPhaseRate.txt", status = "unknown")
      open (unit = 30, file = ".\PhasingResults\IndivPhaseRate.txt", status = "unknown")
      open (unit = 33, file = ".\PhasingResults\FinalHapIndCarry.txt", status = "unknown")
    else
      open (unit = 15, file = "./PhasingResults/FinalPhase.txt", status = "unknown")
      open (unit = 25, file = "./PhasingResults/CoreIndex.txt", status = "unknown")
      open (unit = 28, file = "./PhasingResults/SnpPhaseRate.txt", status = "unknown")
      open (unit = 30, file = "./PhasingResults/IndivPhaseRate.txt", status = "unknown")
      open (unit = 33, file = "./PhasingResults/FinalHapIndCarry.txt", status = "unknown")
    end if

    do i = 1, nAnisG
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') p%getId(i), &
      Phase(i,:, 1)
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') p%getId(i), &
      Phase(i,:, 2)
    end do

    do i = 1, nCores
      write (25, *) i, CoreIndex(i,:)
    end do


    do i = 1, nSnp
      counter = 0
      do j = 1, nAnisG
	if ((Phase(j, i, 1) == 0).or.(Phase(j, i, 1) == 1)) counter = counter + 1
	if ((Phase(j, i, 2) == 0).or.(Phase(j, i, 2) == 1)) counter = counter + 1
      end do
      write (28, '(i10,f7.2)') i, (100 * (float(counter)/(2 * nAnisG)))
    end do

    do i = 1, nAnisG
      l = 0
      do j = 1, nCores
	CounterP = 0
	CounterM = 0
	do k = CoreIndex(j, 1), CoreIndex(j, 2)
	  if ((Phase(i, k, 1) == 0).or.(Phase(i, k, 1) == 1)) counterP = counterP + 1
	  if ((Phase(i, k, 2) == 0).or.(Phase(i, k, 2) == 1)) counterM = counterM + 1
	end do
	l = l + 1
	CoreCount(l) = (float(counterP)/((CoreIndex(j, 2) - CoreIndex(j, 1)) + 1)) * 100
	l = l + 1
	CoreCount(l) = (float(counterM)/((CoreIndex(j, 2) - CoreIndex(j, 1)) + 1)) * 100
      end do
      write (30, '(i10,20000f7.2,20000f7.2,20000f7.2,20000f7.2)') i, CoreCount(:)
    end do

    do i = 1, nAnisG
      k = 0
      do j = 1, nCores
	k = k + 2
	WorkOut(k - 1) = AllHapAnis(i, 1, j)
	WorkOut(k) = AllHapAnis(i, 2, j)
      end do
      write (33, '(a20,20000i8,20000i8,20000i8,20000i8,20000i8)') p%getID(i), WorkOut(:)
      !write (33, '(i10,20000i5,20000i5,20000i5,20000i5,20000i5)') i, WorkOut(:)
    end do
    
    close(15)
    close(25)
    close(28)
    close(30)
    close(33)

  end subroutine WriteOutResults

  subroutine writeOutCore(phase, hapAnis, coreID, coreStart, p)
    use Constants
    use PedigreeDefinition
    implicit none
    
    integer(kind=1), dimension(:,:,:), intent(in) :: phase
    integer, dimension(:,:), intent(in) :: hapAnis
    integer, intent(in) :: coreID
    integer, intent(in) :: coreStart
    type(Pedigree), intent(in) :: p
  
    integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp
    integer, allocatable, dimension(:) :: WorkOut
    double precision, allocatable, dimension(:) :: CoreCount
    
    character(:), allocatable :: coreIDtxt
  
    nAnisG = size(phase,1)
    nSnp = size(phase,2)
    
    allocate(WorkOut(2))
    allocate(CoreCount(2))
    
    coreIDtxt = itoa(coreID)
  
    if (WindowsLinux == 1) then
      open (unit = 15, file = ".\PhasingResults\FinalPhase" // coreIDtxt // ".txt", status = "unknown")
      open (unit = 28, file = ".\PhasingResults\SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      open (unit = 30, file = ".\PhasingResults\IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      open (unit = 33, file = ".\PhasingResults\FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
    else
      open (unit = 15, file = "./PhasingResults/FinalPhase" // coreIDtxt // ".txt", status = "unknown")
      open (unit = 28, file = "./PhasingResults/SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      open (unit = 30, file = "./PhasingResults/IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      open (unit = 33, file = "./PhasingResults/FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
    end if
  
    do i = 1, nAnisG
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') p%getID(i), &
      Phase(i,:, 1)
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') p%getID(i), &
      Phase(i,:, 2)
    end do
  
    do i = 1, nSnp
      counter = 0
      do j = 1, nAnisG
        if ((Phase(j, i, 1) == 0).or.(Phase(j, i, 1) == 1)) counter = counter + 1
        if ((Phase(j, i, 2) == 0).or.(Phase(j, i, 2) == 1)) counter = counter + 1
      end do
      write (28, '(i10,f7.2)') i + CoreStart - 1, (100 * (float(counter)/(2 * nAnisG)))
    end do
  
    do i = 1, nAnisG
      CounterP = 0
      CounterM = 0
      do k = 1, nSnp
	if ((Phase(i, k, 1) == 0).or.(Phase(i, k, 1) == 1)) counterP = counterP + 1
        if ((Phase(i, k, 2) == 0).or.(Phase(i, k, 2) == 1)) counterM = counterM + 1
      end do
      CoreCount(1) = (float(counterP)/(nSnp) * 100)
      CoreCount(2) = (float(counterM)/(nSnp) * 100)
      write (30, '(i10,2f7.2)') i, CoreCount(:)
    end do
  
    do i = 1, nAnisG
      WorkOut(1) = hapAnis(i, 1)
      WorkOut(2) = hapAnis(i, 2)
      write (33, '(i10,2i5)') i, WorkOut(:)
    end do
    
    close(15)
    close(28)
    close(30)
    close(33)
  end subroutine writeOutCore
  
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
  
  subroutine CombineResults(nAnisG, CoreIndex, p)
    use Constants
    use PedigreeDefinition
    implicit none    
        
    integer, intent(in) :: nAnisG
    integer, dimension(:,:), intent(in) :: CoreIndex
    type(Pedigree) :: p
    
    integer :: nCores
    
    integer, dimension(:), allocatable :: inUnits
    integer :: i, j, coreLength, inUnit
    character(:), allocatable :: coreIDtxt
    integer(kind=1), dimension(:), allocatable :: tempPhase
    integer, dimension(2) :: tempHapInd
    double precision, dimension(2) :: tempIndivPhase
    double precision :: tempSnpPhase
    character(len=20) :: id
    
    nCores = size(CoreIndex,1)

    if (WindowsLinux == 1) then
      open (unit = 15, file = ".\PhasingResults\FinalPhase.txt", status = "unknown")
      open (unit = 25, file = ".\PhasingResults\CoreIndex.txt", status = "unknown")
      open (unit = 28, file = ".\PhasingResults\SnpPhaseRate.txt", status = "unknown")
      open (unit = 30, file = ".\PhasingResults\IndivPhaseRate.txt", status = "unknown")
      open (unit = 33, file = ".\PhasingResults\FinalHapIndCarry.txt", status = "unknown")
    else
      open (unit = 15, file = "./PhasingResults/FinalPhase.txt", status = "unknown")
      open (unit = 25, file = "./PhasingResults/CoreIndex.txt", status = "unknown")
      open (unit = 28, file = "./PhasingResults/SnpPhaseRate.txt", status = "unknown")
      open (unit = 30, file = "./PhasingResults/IndivPhaseRate.txt", status = "unknown")
      open (unit = 33, file = "./PhasingResults/FinalHapIndCarry.txt", status = "unknown")
    end if

    do i = 1, nCores
      write (25, *) i, CoreIndex(i,:)
    end do
    
    !!! FINAL PHASE !!!
    allocate(inUnits(nCores))
    do i = 1, nCores
      coreIDtxt = itoa(i)
      if (WindowsLinux == 1) then
	open (newunit = inUnits(i), file = ".\PhasingResults\FinalPhase" // coreIDtxt // ".txt", status = "unknown")
      else
	open (newunit = inUnits(i), file = "./PhasingResults/FinalPhase" // coreIDtxt // ".txt", status = "old")
      end if      
    end do
    
    do i = 1, nAnisG * 2
      do j = 1, nCores
	coreLength  = CoreIndex(j,2) - CoreIndex(j,1) + 1
	allocate(tempPhase(coreLength))
	read(inUnits(j),'(a20,' // itoA(coreLength) // 'i2)') id, tempPhase
	if (j == 1) then
	  write(15, '(a20)', advance = 'no') id
	end if
	if (j == nCores) then
	  write(15, '(' // itoA(coreLength) // 'i2)', advance='yes') tempPhase
	else
	  write(15, '(' // itoA(coreLength) // 'i2)', advance='no') tempPhase
	end if
	deallocate(tempPhase)
      end do
    end do
    
    do i = 1, nCores
      close(inUnits(i))
    end do
    
    deallocate(inUnits)
    
    !!! HAPINDCARRY !!!
    allocate(inUnits(nCores))
    do i = 1, nCores
      coreIDtxt = itoa(i)
      if (WindowsLinux == 1) then
	open (newunit = inUnits(i), file = ".\PhasingResults\FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
      else
	open (newunit = inUnits(i), file = "./PhasingResults/FinalHapIndCarry" // coreIDtxt // ".txt", status = "old")
      end if      
    end do
    
    do i = 1, nAnisG
      do j = 1, nCores
	read(inUnits(j),'(a10,2i5)') id, tempHapInd
	if (j == 1) then
	  write(33, '(a10)', advance = 'no') id
	end if
	if (j == nCores) then
	  write(33, '(2i5)', advance='yes') tempHapInd
	else
	  write(33, '(2i5)', advance='no') tempHapInd
	end if
      end do
    end do
    
    do i = 1, nCores
      close(inUnits(i))
    end do
    
    deallocate(inUnits)
    
    !!! INDIVPHASE !!!
    allocate(inUnits(nCores))
    do i = 1, nCores
      coreIDtxt = itoa(i)
      if (WindowsLinux == 1) then
	open (newunit = inUnits(i), file = ".\PhasingResults\IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      else
	open (newunit = inUnits(i), file = "./PhasingResults/IndivPhaseRate" // coreIDtxt // ".txt", status = "old")
      end if      
    end do
    
    do i = 1, nAnisG
      do j = 1, nCores
	read(inUnits(j),'(a10,2f7.2)') id, tempIndivPhase
	if (j == 1) then
	  write(30, '(a10)', advance = 'no') id
	end if
	if (j == nCores) then
	  write(30, '(2f7.2)', advance='yes') tempIndivPhase
	else
	  write(30, '(2f7.2)', advance='no') tempIndivPhase
	end if
      end do
    end do
    
    do i = 1, nCores
      close(inUnits(i))
    end do
    
    deallocate(inUnits)
    
    !!! SNPPHASE !!!
    do i = 1, nCores
      coreLength  = CoreIndex(i,2) - CoreIndex(i,1) + 1
      coreIDtxt = itoa(i)
      if (WindowsLinux == 1) then
	open (newunit = inUnit, file = ".\PhasingResults\SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      else
	open (newunit = inUnit, file = "./PhasingResults/SnpPhaseRate" // coreIDtxt // ".txt", status = "old")
      end if
      do j = 1, coreLength
	read(inUnit,'(a10,f7.2)') id, tempSnpPhase
	write(28,'(a10,f7.2)') id, tempSnpPhase
      end do
      close(inUnit)
    end do
    
    close(15)
    close(25)
    close(28)
    close(30)
    close(33)

  end subroutine CombineResults
  
  function ParsePedigreeData() result(p)
    use Parameters, only : PedigreeFile, GenotypeFile, nSnp, GenotypeFileFormat
    use Constants
    use PedigreeDefinition
    use NRMCode
    use Sorting
    implicit none
    
    type(Pedigree) :: p

    integer :: i, j, k, SumPseudoNrmS, SumPseudoNrmD, truth, counter, CountMissingGenotype, SireGen, DamGen, nAnisG
    real(kind = 4) :: value, valueS, valueD, SumNrm, SumDiag
    integer, allocatable, dimension (:) :: GenoInPed, WorkVec, ReadingVector
    
    integer, allocatable, dimension (:) :: DanRecode, DanPos !, DanDamGenotyped, DanSireGenotyped
    character(lengan), allocatable, dimension(:) :: DanArray
    integer spos, dpos

    integer :: nAnisRawPedigree
    
    ! Removing Pedigree global variable as first step to moving to seperate subroutine
    character(lengan), allocatable :: ped(:,:)
    !character(lengan), allocatable :: Id(:)
    
    integer(kind = 4), allocatable, dimension (:), target :: SireGenotyped, DamGenotyped
    character(lengan), dimension(:), allocatable :: GenotypeId
    
    call CountInData(nAnisRawPedigree, nAnisG)

    allocate(GenotypeId(nAnisG))
    allocate(GenoInPed(nAnisG))
!    allocate(RecodeGenotypeId(nAnisG))
    allocate(Ped(nAnisRawPedigree, 3))
    allocate(WorkVec(nSnp * 2))
    allocate(ReadingVector(nSnp))
    
    if (trim(PedigreeFile) /= "NoPedigree") then
      open (unit = 2, file = trim(PedigreeFile), status = "old")
      do i = 1, nAnisRawPedigree
	read(2, *) ped(i,:)
      enddo
      close(2)
    else
      open (unit = 3, file = trim(GenotypeFile), status = "old")
      do i = 1, nAnisRawPedigree
	ped(i, 2:3) = "0"
	read (3, *) ped(i, 1)
      enddo
      close (3)
    endif

    GenoInPed = 0

    open (unit = 3, file = trim(GenotypeFile), status = "old")
    
    allocate(DanArray(size(ped, 1)))
    allocate(DanPos(size(ped,1)))
    DanArray = adjustr(ped(:,1))
    call SortWithIndex(DanArray,DanPos)
    
    do i = 1, nAnisG
      truth = 0
      if (GenotypeFileFormat == 1) then
	read (3, *) GenotypeId(i), ReadingVector(:)
      end if
      if (GenotypeFileFormat == 2) then
	!read (3, *) GenotypeId(i), Phase(i,:, 1)
	!read (3, *) GenotypeId(i), Phase(i,:, 2)
	read (3, *) GenotypeId(i), ReadingVector(:)
	read (3, *) GenotypeId(i), ReadingVector(:)
      end if
      if (GenotypeFileFormat == 3) then
	read (3, *) GenotypeId(i), WorkVec(:)
      endif
!      do j = 1, nAnisRawPedigree
!	if (GenotypeId(i) == ped(j, 1)) then
!	  truth = 1
!	  exit
!	endif
!      enddo
      truth = BinarySearch(DanArray,adjustr(GenotypeID(i)))
      if (truth == 0) GenoInPed(i) = 1
    enddo
    deallocate(Ped)
    deallocate(DanArray)
    deallocate(DanPos)
    
    close(3)

    if (trim(PedigreeFile) /= "NoPedigree") then
      nAnisRawPedigree = nAnisRawPedigree + count(GenoInPed(:) == 1)
    else
      nAnisRawPedigree = nAnisG
    endif
    allocate(Ped(nAnisRawPedigree, 3))
    if (trim(PedigreeFile) /= "NoPedigree") then
      open (unit = 2, file = trim(PedigreeFile), status = "old")
      do i = 1, nAnisRawPedigree - count(GenoInPed(:) == 1)
	read(2, *) ped(i,:)
      end do
      counter = nAnisRawPedigree - count(GenoInPed(:) == 1)
      do i = 1, nAnisG
	if (GenoInPed(i) == 1) then
	  counter = counter + 1
	  ped(counter, 1) = GenotypeId(i)
	  ped(counter, 2) = "0"
	  ped(counter, 3) = "0"
	endif
      enddo
      close(2)
    else
      do i = 1, nAnisG
	ped(i, 1) = GenotypeId(i)
	ped(i, 2) = "0"
	ped(i, 3) = "0"
      enddo

    endif

    allocate(SireGenotyped(nAnisG))
    allocate(DamGenotyped(nAnisG))    
    
!    allocate(DanRecode(size(ped,1)))
!    DanRecode = 0
!    do i = 1, nAnisG
!      do j = 1, size(ped,1)
!	if (ped(j,1) == GenotypeId(i)) then
!	  DanRecode(j) = i
!	  exit
!	end if
!      end do
!    end do
!    
!    allocate(SireGenotyped(nAnisG))
!    allocate(DamGenotyped(nAnisG))
    
!    SireGenotyped = 0
!    DamGenotyped = 0
!    do i = 1, size(ped,1)
!      if (DanRecode(i) /= 0) then
!	do j = 1, nAnisG
!	  if (GenotypeID(j) .eq. ped(i,2)) then
!	    SireGenotyped(DanRecode(i)) = j
!	  end if
!	  if (GenotypeID(j) .eq. ped(i,3)) then
!	    DamGenotyped(DanRecode(i)) = j
!	  end if
!	end do
!      end if
!    end do
!
!    deallocate(DanRecode)
    
    allocate(DanArray(size(GenotypeID)))
    allocate(DanPos(size(GenotypeID)))
    DanArray = adjustr(GenotypeID)
    call SortWithIndex(DanArray,DanPos)
    
    allocate(DanRecode(size(ped,1)))
    do i = 1, size(ped,1)
      dpos = BinarySearch(DanArray, adjustr(ped(i,1)))
      if (dpos > 0) then
	DanRecode(i) = DanPos(dpos)
      else
	DanRecode(i) = 0
      endif
    end do
    
    SireGenotyped = 0
    DamGenotyped = 0
    do i = 1, size(ped,1)
      if (DanRecode(i) /= 0) then
	spos = BinarySearch(DanArray, adjustr(ped(i,2)))
	if (spos > 0) then
	  SireGenotyped(DanRecode(i)) = DanPos(spos)
	endif
	dpos = BinarySearch(DanArray, adjustr(ped(i,3)))
	if (dpos > 0) then
	  DamGenotyped(DanRecode(i)) = DanPos(dpos)
	endif
      end if
    end do

    deallocate(DanArray)
    deallocate(DanRecode)
    deallocate(DanPos)

    allocate (nrmped(size(ped,1),size(ped,2)))
    nrmped = ped
    
    p = Pedigree(sireGenotyped, damGenotyped, genotypeId)
    
    deallocate(sireGenotyped, damGenotyped, genotypeID)
  end function ParsePedigreeData
  
  function ParseGenotypeData(startSnp, endSnp, nAnisG) result(Genos)
    use Parameters, only: GenotypeFile, nSnp, GenotypeFileFormat
    use Constants
    implicit none

    integer, intent(in) :: startSnp, endSnp
    integer, intent(in) :: nAnisG
    integer(kind=1), allocatable, dimension(:,:) :: Genos

    integer :: i, j, k
    integer, allocatable, dimension (:) :: WorkVec, ReadingVector
    integer :: nReadSnp
    character(lengan) :: dummy
    
    open (unit = 3, file = trim(GenotypeFile), status = "old")

    nReadSnp = endSnp - startSnp + 1

    allocate(Genos(nAnisG, nReadSnp))
    Genos = MissingGenotypeCode
    
    !allocate(Phase(nAnisG, nReadSnp, 2))
    allocate(WorkVec(nSnp * 2))
    allocate(ReadingVector(nSnp))

    !allocate(HapLib(nAnisG * 2, nSnp))

    !Phase = 9

    Genos = MissingGenotypeCode
    do i = 1, nAnisG
      if (GenotypeFileFormat == 1) then
	read (3, *) dummy, ReadingVector(:)
	!do j = 1, nSnp
	do j = startSnp, endSnp
	  if ((ReadingVector(j) /= 0).and.(ReadingVector(j) /= 1).and.(ReadingVector(j) /= 2)) ReadingVector(j) = MissingGenotypeCode
	  Genos(i, j - startSnp + 1) = ReadingVector(j)
!	  if (Genos(i, j) == 0) Phase(i, j,:) = 0
!	  if (Genos(i, j) == 2) Phase(i, j,:) = 1
	end do
      end if
!      if (GenotypeFileFormat == 2) then
!	!read (3, *) GenotypeId(i), Phase(i,:, 1)
!	!read (3, *) GenotypeId(i), Phase(i,:, 2)
!	read (3, *) dummy, ReadingVector(:)
!	Phase(i,:,1) = ReadingVector(startSnp:endSnp)
!	read (3, *) dummy, ReadingVector(:)
!	Phase(i,:,2) = ReadingVector(startSnp:endSnp)
!      end if
      if (GenotypeFileFormat == 3) then
	read (3, *) dummy, WorkVec(:)
	k = 0
	!do j = 1, nSnp * 2
	do j = startSnp*2-1,endSnp*2
	  if (mod(j, 2) == 0) then
	    k = k + 1
	    if ((WorkVec(j - 1) == 1).and.(WorkVec(j) == 1)) Genos(i, k) = 0
	    if ((WorkVec(j - 1) == 1).and.(WorkVec(j) == 2)) Genos(i, k) = 1
	    if ((WorkVec(j - 1) == 2).and.(WorkVec(j) == 1)) Genos(i, k) = 1
	    if ((WorkVec(j - 1) == 2).and.(WorkVec(j) == 2)) Genos(i, k) = 2
	  endif
	end do
      endif
    enddo
    
    close(3)
  end function ParseGenotypeData
  
  function ParsePhaseData(startSnp, endSnp, nAnisG) result(Phase)
    use Parameters, only: GenotypeFile, nSnp, GenotypeFileFormat
    use Constants
    implicit none

    integer, intent(in) :: startSnp, endSnp
    integer, intent(in) :: nAnisG
    integer(kind=1), allocatable, dimension(:,:,:) :: Phase

    integer :: i, j, k
    integer, allocatable, dimension (:) :: WorkVec, ReadingVector
    integer :: nReadSnp
    character(lengan) :: dummy
    
    open (unit = 3, file = trim(GenotypeFile), status = "old")

    nReadSnp = endSnp - startSnp + 1

    allocate(Phase(nAnisG, nReadSnp, 2))
    allocate(ReadingVector(nSnp))

    Phase = 9

    do i = 1, nAnisG
      read (3, *) dummy, ReadingVector(:)
      Phase(i,:,1) = ReadingVector(startSnp:endSnp)
      read (3, *) dummy, ReadingVector(:)
      Phase(i,:,2) = ReadingVector(startSnp:endSnp)
    enddo
    
    close(3)
  end function ParsePhaseData

  subroutine ReadInParameterFile(filename)
    use Parameters
    implicit none
    
    character(*), intent(in) :: filename

    double precision :: PercSurrDisagree
    integer :: i, TempInt, Graphics, status, cl
    character (len = 300) :: dumC, FileFormat, OffsetVariable, hold

    !open (unit = 1, file = "AlphaPhaseSpec.txt", status = "old")
    open (unit = 1, file = filename, status = "old")

    read (1, *) dumC, PedigreeFile
    read (1, *) dumC, GenotypeFile, FileFormat
    if (trim(FileFormat) == 'GenotypeFormat') then
      GenotypeFileFormat = 1
    elseif (trim(FileFormat) == 'PhaseFormat') then
      GenotypeFileFormat = 2
    elseif (trim(FileFormat) == 'UnorderedFormat') then
      GenotypeFileFormat = 3
    else
      print*, "The genotype file format is not correctly specified"
      stop
    endif

    print *, " Parameter file read"
    print *, " "
    print *, " Using ", trim(PedigreeFile), " as the pedigree file"
    print *, " Using ", trim(GenotypeFile), " as the genotype file"
    print *, " "

    read (1, *) dumC, nSnp
    read (1, *) dumC, CoreAndTailLength
    if (CoreAndTailLength > nSnp) then
      print*, "GeneralCoreAndTailLength is too long"
      stop
    endif
    read (1, *) dumC, Jump, OffsetVariable
    if (Jump > nSnp) then
      print*, "GeneralCoreLength is too long"
      stop
    endif

    if (CoreAndTailLength < Jump) then
      print *, "GeneralCoreAndTailLength is shorted than GenralCoreLength"
      stop
    end if

    if (OffsetVariable == "Offset") then
      Offset = 1
    endif
    if (OffsetVariable == "NotOffset") then
      Offset = 0
    endif

    if ((OffsetVariable /= "Offset").and.(OffsetVariable /= "NotOffset")) then
      print*, "Offset variable is not properly parameterised"
      stop
    endif
    read (1, *) dumC, UseSurrsN
    read (1, *) dumC, PercSurrDisagree
    read (1, *) dumC, PercGenoHaploDisagree
    read (1, *) dumC, GenotypeMissingErrorPercentage
    read (1, *) dumC, NrmThresh
    read (1, *) dumC, FullFileOutput
    read (1, *) dumC, Graphics
    read (1, *) dumC, Simulation
    read (1, *) dumC, TruePhaseFile

    read (1, *, iostat=status) dumC, tempInt
    if (status == 0) then
      readCoreAtTime = (tempInt==1)
    else
      readCoreAtTime = .false.
    end if

    read (1, *, iostat=status) dumC, itterateType
    if (status == 0) then
      consistent = ((itterateType .eq. "Off") .and. (.not. readCoreAtTime))
    else
      consistent = .true.
      itterateType = "Off"
    end if

!    read (1, *, iostat=status) dumC, itterateNumber
!    if (status /= 0) then
!      itterateNumber = 200
!    end if
!    read (1, *, iostat=status) dumC, numIter
!    if (status /= 0) then
!      numIter = 1
!    end if
    read(1, *, iostat=status) dumC, hold
    if (status == 0) then
      if (hold(1:1) == "*") then
	read(hold,"(X,I2)") cl
	call get_command_argument(cl,hold)
	read(hold,*) itterateNumber
      else
	read(hold,*) itterateNumber
      end if
    else
      itterateNumber = 200
    end if
    
    read(1, *, iostat=status) dumC, hold
    if (status == 0) then
      if (hold(1:1) == "*") then
	read(hold,"(X,I2)") cl
	call get_command_argument(cl,hold)
	read(hold,*) numIter
      else
	read(hold,*) numIter
      end if
    else
      numIter = 1
    end if
    if (itterateType .eq. "Off") then
      numIter = 1
    end if

    read (1, *, iostat=status) dumC, startCoreChar, endCoreChar
    if (status /= 0) then
      startCoreChar = "1"
      endCoreChar = "Combine"
    end if
    
    read (1, *, iostat=status) dumC, hold
    if (status == 0) then
      if (hold(1:1) == "*") then
	read(hold,"(X,I2)") cl
	call get_command_argument(cl,hold)
	read(hold,*) minHapFreq
      else
	read(hold,*) minHapFreq
      end if
    else
      minHapFreq = 1
    end if

    PercSurrDisagree = PercSurrDisagree/100
    NumSurrDisagree = int(UseSurrsN * PercSurrDisagree)
    PercGenoHaploDisagree = PercGenoHaploDisagree/100
    GenotypeMissingErrorPercentage = GenotypeMissingErrorPercentage/100

    close (1)

    if (Graphics == 1) then
      print*, "Graphics option is not yet functional"
      stop
    end if

    !if (nSnp>32767) then
    !        print*, "Kind=2 is not sufficient for this number of SNP.... Contact John Hickey because there is a simple solution!"
    !        stop
    !end if

  end subroutine ReadInParameterFile

  subroutine HapCommonality(library, OutputPoint)
    use Parameters, only: FullFileOutput
    use Constants
    use HaplotypeLibraryDefinition
    implicit none

    type(HaplotypeLibrary), intent(in) :: library
    integer, intent(in) :: OutputPoint

    integer :: i, SizeCore, nHaps
    character(len = 300) :: filout

    integer, allocatable, dimension (:,:) :: HapRel

    SizeCore = library%getNumSnps()
    nHaps = library%getSize()

    if (FullFileOutput == 1) then
      if (WindowsLinux == 1) then
	write (filout, '(".\PhasingResults\HaplotypeLibrary\Extras\HapCommonality",i0,".txt")') OutputPoint
	open (unit = 27, FILE = filout, status = 'unknown')

      else
	write (filout, '("./PhasingResults/HaplotypeLibrary/Extras/HapCommonality",i0,".txt")') OutputPoint
	open (unit = 27, FILE = filout, status = 'unknown')

      endif
    endif

    HapRel = library%getHapRel()

    if (FullFileOutput == 1) then
      do i = 1, nHaps
	write (27, '(i10,20000F5.2,20000F5.2,20000F5.2,20000F5.2)') i, float(HapRel(i,:))/SizeCore
      enddo
    endif

    close(27)

  end subroutine HapCommonality

  !########################################################################################################################################################################

  subroutine WriteSurrogates(definition, threshold, OutputPoint, p)
    use SurrogateDefinition
    use PedigreeDefinition
    use Parameters, only : FullFileOutput
    use Constants

    implicit none

    character(len = 300) :: filout
    integer :: i, j, nSurrogates

    type(Surrogate), intent(in) :: definition
    integer, intent(in) :: threshold
    integer, intent(in) :: OutputPoint
    type(Pedigree), intent(in) :: p

    integer :: nAnisG

    nAnisG = size(definition%numOppose,1)

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
	  write (13, '(a20,20000i6,20000i6,20000i6,20000i6)') p%getID(i), definition%partition(i,:)
	else
	  write (13, *) p%getID(i), definition%partition(i,:)
	end if
	do j = 1, nAnisG
	  if (definition%numOppose(i, j) <= threshold) nSurrogates = nSurrogates + 1
	enddo
	write (19, '(a20,20000i6,20000i6,20000i6,20000i6)') &
	p%getID(i), count(definition%partition(i,:) == 1), count(definition%partition(i,:) == 2)&
	, count(definition%partition(i,:) == 3), nSurrogates, definition%method(i)
      enddo
      close(13)
      close(19)
    end if

  end subroutine WriteSurrogates

!#################################################################################################################################################################

  subroutine CountInData(nAnisRawPedigree, nAnisG)
    use Parameters, only: PedigreeFile, GenotypeFile, GenotypeFileFormat
    implicit none

    integer, intent(out) :: nAnisRawPedigree, nAnisG

    integer :: k
    character (len = 300) :: dumC

    nAnisRawPedigree = 0
    if (trim(PedigreeFile) /= "NoPedigree") then
      open (unit = 2, file = trim(PedigreeFile), status = "old")
      do
	read (2, *, iostat = k) dumC
	nAnisRawPedigree = nAnisRawPedigree + 1
	if (k /= 0) then
	  nAnisRawPedigree = nAnisRawPedigree - 1
	  exit
	endif
      enddo
      close(2)
      print*, " ", nAnisRawPedigree, " individuals in the pedigree file"
    endif

    nAnisG = 0
    open (unit = 3, file = trim(GenotypeFile), status = "old")
    do
      read (3, *, iostat = k) dumC
      nAnisG = nAnisG + 1
      if (k /= 0) then
	nAnisG = nAnisG - 1
	exit
      endif
    enddo
    if (GenotypeFileFormat == 2) then
      nAnisG = nAnisG /2
    end if
    close(3)
    print*, " ", nAnisG, " individuals in the genotype file"

    if (trim(PedigreeFile) == "NoPedigree") nAnisRawPedigree = nAnisG

  end subroutine CountInData

!########################################################################################################################################################################################################

  subroutine MakeDirectories
    use Parameters, only: FullFileOutput, Simulation
    use Constants
    implicit none

    print*, ""
    if (WindowsLinux == 1) then
      call system("rmdir /s /q Miscellaneous")
      call system("rmdir /s /q PhasingResults")
    else
      call system("rm -r Miscellaneous")
      call system("rm -r PhasingResults")
    endif

    call system("mkdir PhasingResults")
    if (WindowsLinux == 1) then
      call system("mkdir PhasingResults\HaplotypeLibrary")
      if (FullFileOutput == 1) call system("mkdir PhasingResults\HaplotypeLibrary\Extras")
      !open (unit = 29, file = ".\PhasingResults\PhasingYield.txt", status = "unknown")
    else
      call system("mkdir PhasingResults/HaplotypeLibrary")
      if (FullFileOutput == 1) call system("mkdir PhasingResults/HaplotypeLibrary/Extras")
      !open (unit = 29, file = "./PhasingResults/PhasingYield.txt", status = "unknown")
    endif


    if (WindowsLinux == 1) then
      if (FullFileOutput == 1) then
	call system("mkdir Miscellaneous")
      endif
    else
      if (FullFileOutput == 1) then
      call system("mkdir Miscellaneous")
      endif
    endif
    
    if (Simulation == 1) then
      if (WindowsLinux==1) then
	call system("rmdir /s /q Simulation")
      else
	call system("rm -r Simulation")
      endif

      if (WindowsLinux==1) then
	if (FullFileOutput==1) then
	  call system("mkdir Simulation")
	endif
      else
	if (FullFileOutput==1) then
	  call system("mkdir Simulation")
	endif
      endif
    end if

  end subroutine MakeDirectories  
  
  subroutine WriteHapLib(library, currentcore, c)
    use Parameters, only: fullfileoutput, ItterateType
    use HaplotypeLibraryDefinition
    use Constants
    use CoreDefinition
    implicit none

    type(HaplotypeLibrary), intent(in) :: library
    type(Core), intent(in) :: c
    integer, intent(in) :: currentcore

    integer :: i, j, k, counter, SizeCore, nHaps !, nAnisG
    character(len = 300) :: filout

    SizeCore = library%getNumSnps()

    nHaps = library%getSize()

    if (FullFileOutput == 1) then
      if (WindowsLinux == 1) then
	write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".txt")') currentcore
	open (unit = 24, FILE = filout, status = 'unknown')
      else
	write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".txt")') currentcore
	open (unit = 24, FILE = filout, status = 'unknown')
      endif
    endif
    if (WindowsLinux == 1) then
      write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".bin")') currentcore
      open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
    else
      write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') currentcore
      open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
    endif


    write (34) nHaps, SizeCore
    do i = 1, nHaps
      if (FullFileOutput == 1)&
      write (24, '(2i6,a2,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)') &
      i, library%getHapFreq(i), " ", library%getHap(i)
      write (34) library%getHap(i)
    end do

    if (FullFileOutput == 1) then
      close(24)
    end if
    close(34)

    if (ItterateType .eq. "Off") then
      print*, "   ", "Final iteration found ", nHaps, "haplotypes"
    
      print*, ""
      write (*, '(a4,a30,f5.2,a1)') "   ", "Final yield for this core was ", c%getTotalYield(), "%"
    end if
    
    if (WindowsLinux == 1) then
      open (unit = 29, file = ".\PhasingResults\PhasingYield.txt", status = "unknown", position = "append")
    else
      open (unit = 29, file = "./PhasingResults/PhasingYield.txt", status = "unknown", position = "append")
    endif

    write (29, '(i10,f7.2)') CurrentCore, c%getTotalYield()

    close(29)

  end subroutine WriteHapLib
  
  subroutine InsertionSort(array, pos)
    character(*), dimension(:), intent(inout) :: array
    integer, dimension(size(array)), intent(out) :: pos
    
    integer :: i, j, tempP
    character(len(array)) :: tempA
    
    do i = 1, size(array)
      pos(i) = i
    end do
    
    do i = 2, size(array)
      j = i
      do while (array(j-1) > array(j))
	tempA = array(j-1)
	array(j-1) = array(j)
	array(j) = tempA
	
	tempP = pos(j-1)
	pos(j-1) = pos(j)
	pos(j) = tempP
	
	j = j - 1
	if (j == 1) exit
      end do
    end do
  end subroutine InsertionSort
  
  function BinarySearch(array, val) result(pos)
    character(*), dimension(:), intent(in) :: array
    character(*), intent(in) :: val
    integer :: pos
    
    integer :: low, high, mid
    
    low = 1
    high = size(array)
    
    do while (low <= high)
      mid = (low + high) / 2
      if (array(mid) > val) then
	high = mid - 1
      end if
      if (array(mid) < val) then
	low = mid + 1
      end if
      if (array(mid) == val) then
	pos = mid
	return 
      end if
    end do
    
    pos = 0
  end function
  
end module InputOutput