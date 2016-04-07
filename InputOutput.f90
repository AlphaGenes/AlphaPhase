module InputOutput
  
contains

  subroutine WriteOutResults(phase, allHapAnis)
    use Global, only: CoreIndex, WindowsLinux, nCores, GenotypeID
    implicit none

    integer(kind=1), dimension(:,:,:), intent(in) :: phase
    integer, dimension(:,:,:), intent(in) :: allHapAnis

    integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp
    integer, allocatable, dimension(:) :: WorkOut
    double precision, allocatable, dimension(:) :: CoreCount

    nAnisG = size(phase,1)
    nSnp = size(phase,2)

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
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
      Phase(i,:, 1)
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
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
      write (33, '(i10,20000i5,20000i5,20000i5,20000i5,20000i5)') i, WorkOut(:)
    end do
    
    close(15)
    close(25)
    close(28)
    close(30)
    close(33)

  end subroutine WriteOutResults

  subroutine writeOutCore(phase, hapAnis, coreID)
    use Global, only: CoreIndex, WindowsLinux, nCores, GenotypeID
    implicit none
    
    integer(kind=1), dimension(:,:,:), intent(in) :: phase
    integer, dimension(:,:), intent(in) :: hapAnis
    integer, intent(in) :: coreID
  
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
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
      Phase(i,:, 1)
      write(15, '(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') GenotypeId(i), &
      Phase(i,:, 2)
    end do
  
    do i = 1, nSnp
      counter = 0
      do j = 1, nAnisG
        if ((Phase(j, i, 1) == 0).or.(Phase(j, i, 1) == 1)) counter = counter + 1
        if ((Phase(j, i, 2) == 0).or.(Phase(j, i, 2) == 1)) counter = counter + 1
      end do
      write (28, '(i10,f7.2)') i + CoreIndex(coreID,1) - 1, (100 * (float(counter)/(2 * nAnisG)))
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
  
  subroutine CombineResults(nAnisG)
    use Global, only: CoreIndex, WindowsLinux, nCores, GenotypeID
    implicit none    
        
    integer, intent(in) :: nAnisG
    
    integer, dimension(:), allocatable :: inUnits
    integer :: i, j, coreLength, inUnit
    character(:), allocatable :: coreIDtxt
    integer(kind=1), dimension(:), allocatable :: tempPhase
    integer, dimension(2) :: tempHapInd
    double precision, dimension(2) :: tempIndivPhase
    double precision :: tempSnpPhase
    character(len=20) :: id

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
	open (unit = 33, file = ".\PhasingResults\FinalPhase" // coreIDtxt // ".txt", status = "unknown")
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
	open (unit = 33, file = ".\PhasingResults\FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
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
	open (unit = 33, file = ".\PhasingResults\IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
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
	open (unit = 33, file = ".\PhasingResults\SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
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
  
  subroutine ParsePedigreeData
    use GlobalPedigree
    use Global
    implicit none

    integer :: i, j, k, SumPseudoNrmS, SumPseudoNrmD, truth, counter, CountMissingGenotype, SireGen, DamGen
    real(kind = 4) :: value, valueS, valueD, SumNrm, SumDiag
    integer, allocatable, dimension (:) :: GenoInPed, WorkVec, ReadingVector

    ! Removing Pedigree global variable as first step to moving to seperate subroutine
    character(lengan), allocatable :: ped(:,:)
    character(lengan), allocatable :: Id(:)

    interface PedigreeViewerRecode
      subroutine PedigreeViewerRecode(nobs, nAnisPedigree, ped, id)
	use GlobalPedigree
	implicit none

	integer :: nobs, nAnisPedigree
	character(len = lengan), dimension(:,:), intent(in) :: ped
	character(lengan), allocatable :: Id(:)
      end subroutine PedigreeViewerRecode
    end interface PedigreeViewerRecode

    allocate(GenotypeId(nAnisG))
    allocate(GenoInPed(nAnisG))
    allocate(RecodeGenotypeId(nAnisG))
    allocate(Ped(nAnisRawPedigree, 3))
    allocate(WorkVec(nSnp * 2))
    allocate(ReadingVector(nSnp))
    
    open (unit = 3, file = trim(GenotypeFile), status = "old")

    !allocate(HapLib(nAnisG * 2, nSnp))

    if (trim(PedigreeFile) /= "NoPedigree") then
      open (unit = 2, file = trim(PedigreeFile), status = "old")
      do i = 1, nAnisRawPedigree
	read(2, *) ped(i,:)
      enddo
      close(2)
    else
      do i = 1, nAnisRawPedigree
	ped(i, 2:3) = "0"
	read (3, *) ped(i, 1)
      enddo
      rewind (3)
    endif

    GenoInPed = 0

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
      do j = 1, nAnisRawPedigree
	if (GenotypeId(i) == ped(j, 1)) then
	  truth = 1
	  exit
	endif
      enddo
      if (truth == 0) GenoInPed(i) = 1
    enddo
    deallocate(Ped)

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
    call PedigreeViewerRecode(nAnisRawPedigree, nAnisP, ped, id)
    deallocate(Ped)
    deallocate(GenoInPed)

    do i = 1, nAnisG
      do j = 1, nAnisP
	if (Id(j) == GenotypeId(i)) then
	  RecodeGenotypeId(i) = j
	  exit
	end if
      end do
    end do
    rewind (3)

    allocate(RecSire(0:nAnisP))
    allocate(RecDam(0:nAnisP))
    RecSire(0) = 0
    RecDam(0) = 0
    do i = 1, nAnisP
      RecSire(i) = seqsire(i)
      RecDam(i) = seqdam(i)
    enddo

    allocate(SireGenotyped(nAnisG))
    allocate(DamGenotyped(nAnisG))

    do i = 1, nAnisG
      SireGenotyped(i) = seqsire(RecodeGenotypeId(i))
      DamGenotyped(i) = seqdam(RecodeGenotypeId(i))
      truth = 0
      do j = 1, nAnisG
	if (SireGenotyped(i) == RecodeGenotypeId(j)) then
	  truth = 1
	  SireGenotyped(i) = j
	  exit
	endif
      enddo
      if (truth == 0) SireGenotyped(i) = 0
      truth = 0
      do j = 1, nAnisG
	if (DamGenotyped(i) == RecodeGenotypeId(j)) then
	  truth = 1
	  DamGenotyped(i) = j
	  exit
	endif
      enddo
      if (truth == 0) DamGenotyped(i) = 0
    enddo

    deallocate(seqid)
  end subroutine ParsePedigreeData
  
  function ParseGenotypeData(startSnp, endSnp) result(Genos)
    use GlobalPedigree
    use Global
    implicit none

    integer, intent(in) :: startSnp, endSnp
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

    
end module InputOutput