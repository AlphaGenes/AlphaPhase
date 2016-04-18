!     Last change:  JM   21 Jun 2011   11:54 am
!v On the 6th of May 2015 Roberton Antolin and John Hickey agreed that the version Roberto has becomes the master version. The expectation is that this should include all changes made by John when visiting Holland and prior to Roberto starting to work here. John will keep his old version in case this turns out not to be true
!v Roberto will take the fastPHASE HMM and add it to this code and also do the shuffling of markers Speed up to this code. 
!v We have diffed the version of the first commit and John's vesion and they are exactly the same file (i.e., Roberto's master is 100% certain to include Johns original edits in holland etc)
!####################################################################################################################################################################

module Global
  implicit none

  !!!!! CONSTANTS !!!!!
  integer, parameter :: WindowsLinux = 0 !If 1 then compile for Windows / If 0 then compile for Linux
  integer, parameter :: lengan = 20, MissingGenotypeCode = 3
  integer, parameter :: nMAxRounds = 100
  
  !!!!! INPUT PARAMETERS !!!!!
  character(len=300) GenotypeFile
  integer :: GenotypeFileFormat
  integer :: nSnp   ! Possibly doesn't need to be a parameter  
  integer :: CoreAndTailLength
  integer :: Jump, Offset
  integer :: UseSurrsN
  integer :: NumSurrDisagree
  double precision :: PercGenoHaploDisagree
  double precision :: GenotypeMissingErrorPercentage
  double precision :: NrmThresh
  integer :: FullFileOutput
  integer :: Graphics
  integer :: Simulation
  character (len = 300) :: PedigreeFile ! Used in a really weird way that should probably be refactored
  
  logical :: readCoreAtTime
  character (len = 300) :: itterateType
  integer :: itterateNumber
  integer :: numIter
  character (len = 10) :: startCoreChar, endCoreChar
  
  logical :: consistent
  
  !!!!! INPUT DATA !!!!!
  integer(kind = 4), allocatable, dimension (:) :: SireGenotyped, DamGenotyped
  !integer(kind = 1), allocatable, dimension (:,:) :: Genos
  character(lengan), allocatable :: GenotypeId(:)

  
  !!!!! MOVE OUT? !!!!!
  !integer :: nAnisG, nAnisRawPedigree, nAnisP, nCores
  integer :: nAnisG, nAnisP 

  !integer(kind = 1), allocatable, dimension (:,:,:) :: Phase
  !integer, allocatable, dimension (:,:) :: CoreIndex, TailIndex
  !integer :: nGlobalHapsIter

end module Global

!####################################################################################################################################################################

!module GlobalClustering
!  implicit none

!  integer, parameter :: nCLusters = 2, nMaxRounds = 100
!  integer, parameter :: nMaxRounds = 100
!  integer :: Change, rounds, SurrCounter
!  integer :: Change
!  double precision :: CountCluster(nClusters), Dist(nClusters)
!  integer, allocatable, dimension (:,:) :: TempSurrArray
!  integer, allocatable, dimension (:) :: TempSurrVector, ClusterMember
!  integer, allocatable, dimension (:) :: ClusterMember
!  double precision, allocatable, dimension (:,:) :: Medoids
!  double precision, allocatable, dimension (:) :: MinClust

!end module GlobalClustering

!####################################################################################################################################################################

!module GlobalClusteringHaps
!  implicit none
!
!  integer, parameter :: nCLusters = 2, nMaxRounds = 100
!  integer :: Change, rounds, SurrCounter, nHapsCluster, SnpInCore
!  double precision :: CountCluster(nClusters), Dist(nClusters)
!  integer, allocatable, dimension (:,:) :: TempHapArray
!  integer, allocatable, dimension (:) :: TempHapVector, ClusterMember
!  double precision, allocatable, dimension (:,:) :: Medoids
!  double precision, allocatable, dimension (:) :: MinClust
!
!end module GlobalClusteringHaps

!####################################################################################################################################################################

module GlobalPedigree
  use Global
  implicit none

  !real(kind = 4), allocatable :: xnumrelmatHold(:)
  !integer :: NRMmem, shell, shellmax, shellWarning
  !integer :: NRMmem, shell
  !integer, allocatable :: seqid(:), seqsire(:), seqdam(:), RecodeGenotypeId(:), RecSire(:), RecDam(:)
  integer, allocatable :: seqsire(:), seqdam(:), RecodeGenotypeId(:)!, RecSire(:), RecDam(:)
  !character(lengan), allocatable :: ped(:,:), Id(:)
end module GlobalPedigree

!####################################################################################################################################################################

program Rlrplhi
  use Global
  use HaplotypeLibrary
  use SurrogateDefinition
  use CoreSubsetDefinition
  use Phasing
  use CoreDefinition
  use NRMcode
  use InputOutput
  use MemberManagerDefinition
  implicit none

  integer :: h, i, j, counter, SizeCore, nGlobalHapsOld, nCount, threshold
  double precision :: value, Yield
  
  type(HapLib) :: library
  type(SurrDef) :: surrogates
  type(Core) :: c
  type(CoreSubset) :: cs
  
  integer, allocatable, dimension (:,:,:) :: AllHapAnis
  integer(kind=1), allocatable, dimension(:,:,:) :: AllPhase
  integer(kind=1), allocatable, dimension(:,:) :: AllGenos, Genos
  integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp
  integer, allocatable, dimension (:,:) :: CoreIndex, TailIndex
  integer :: nCores
  integer :: nGlobalHapsIter
  
  integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM
  
  integer :: startCore, endCore
  logical :: combine
  
  type(MemberManager) :: manager
  
  interface calculateCores
    subroutine calculateCores(CoreIndex, TailIndex)
      use Global
      implicit none

      integer, dimension(:,:), allocatable :: CoreIndex, TailIndex
    end subroutine calculateCores
  end interface calculateCores
  
  ! Create a seed for RNG
  ! call system_clock(nCount)
  ! secs = mod(nCount,int(1e6))
  
  readCoreAtTime = .false.
  
  call Titles
  call ReadInParameterFile
  call MakeDirectories
  call ParsePedigreeData
  call CalculateCores(CoreIndex, TailIndex)
  nCores = size(CoreIndex,1)
  if (.not. readCoreAtTime) then
    allocate(AllGenos(nAnisG,nSnp))
    AllGenos = ParseGenotypeData(1,nSnp)
  end if
  
  if (consistent) then
    allocate(PseudoNRM(nAnisG,nAnisG))
    PseudoNRM = createNRM()
  end if
  
  if (.not. readCoreAtTime) then
    allocate(AllHapAnis(nAnisG, 2, nCores))
    allocate(AllPhase(nAnisG, nSnp, 2))
  end if
    
  threshold = int(GenotypeMissingErrorPercentage*CoreAndTailLength)
  
  if (startCoreChar .eq. "Combine") then
    startCore = nCores
    combine = .true.
  else
    read(startCoreChar, '(i10)') startCore
    combine = .false.
  end if
  
  if (endCoreChar .eq. "Combine") then
    endCore = nCores
    combine = .true.
  else
    read(endCoreChar, '(i10)') endCore
    combine = .false.
  end if
    

  !do h = 1, nCores
  do h = startCore, endCore
    !CurrentCore = h
    !nGlobalHaps = 0
    print*, " "
    print*, " "
    print*, " Starting Core", h
    !OutputPoint = h
    StartCoreSnp = CoreIndex(h, 1)
    EndCoreSnp = CoreIndex(h, 2)
    StartSurrSnp = TailIndex(h, 1)
    EndSurrSnp = TailIndex(h, 2)
    
    allocate(Genos(nAnisG, max(EndSurrSnp,EndCoreSnp)-startSurrSnp+1))
    if (readCoreAtTime) then
      Genos = ParseGenotypeData(StartSurrSnp,max(EndSurrSnp,EndCoreSnp))
    else
      Genos = AllGenos(:,StartSurrSnp:max(EndSurrSnp,EndCoreSnp))
    end if
    ! Fudge below
    call c%create(Genos, startCoreSnp-startSurrSnp+1, endCoreSnp-startSurrSnp+1, endSurrSnp-startSurrSnp+1)    
    !call manager%createInputOrder(c, 200)
    call manager%create(c)
    
    do while (manager%hasNext())
      call cs%create(c, SireGenotyped, DamGenotyped, manager%getNext())
    
      call surrogates%calculate(cs, threshold, consistent, pseudoNRM)
      call writeSurrogates(surrogates,threshold, h)
      call Erdos(surrogates, threshold, cs)
      call CheckCompatHapGeno(cs)      
      
      call library%initalise(EndCoreSnp-StartCoreSnp+1,500,500)      
      nGlobalHapsIter = 1
      call MakeHapLib(library, c)
      !nGlobalHapsOld = nGlobalHaps
      nGlobalHapsOld = library%getSize()
      print*, " "
      print*, "  ", "Haplotype library imputation step"
      do j = 1, 20
	call ImputeFromLib(library, c, nGlobalHapsIter)
	call MakeHapLib(library, c)
	!if (nGlobalHapsOld == nGlobalHaps) exit
	!nGlobalHapsOld = nGlobalHaps
	if (nGlobalHapsOld == library%getSize()) exit
	nGlobalHapsOld = library%getSize()
      end do
      
    end do  
      
    call WriteHapLib(library, h, c)
    
    call HapCommonality(library, h)
    
    if (readCoreAtTime) then
      call WriteOutCore(c%phase,c%hapAnis, h, CoreIndex(h,1))
    end if
    
    deallocate(Genos)
    
    AllPhase(:,startCoreSnp:endCoreSnp,:) = c%phase
    AllHapAnis(:,1,h) = c%hapAnis(:,1)
    AllHapAnis(:,2,h) = c%hapAnis(:,2)

    ! HIDDEN MARKOV MODEL SHOULD COME HERE

    !call RationaliseLibrary
    !if (Simulation == 1) then
    !  call Flipper
    !  call Checker
    !end if
  end do
  
  if (readCoreAtTime .and. combine) then
    call CombineResults(nAnisG,CoreIndex)
  else
    call WriteOutResults(AllPhase,AllHapAnis,CoreIndex)
  end if
  if (consistent) then
    deallocate(PseudoNRM)
  end if
  call PrintTimerTitles

end program Rlrplhi

!####################################################################################################################################################################

subroutine ReadInParameterFile
  use Global
  implicit none
  
  double precision :: PercSurrDisagree
  integer :: i, TempInt
  character (len = 300) :: dumC, TruePhaseFile, FileFormat, OffsetVariable
  
  open (unit = 1, file = "AlphaPhaseSpec.txt", status = "old")

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
  
  read (1, *) dumC, tempInt
  readCoreAtTime = (tempInt==1)
  
  read (1, *) dumC, itterateType
  consistent = ((itterateType .eq. "Off") .and. (.not. readCoreAtTime))
  
  read (1, *) dumC, itterateNumber
  read (1, *) dumC, numIter
  
  read (1, *) dumC, startCoreChar, endCoreChar
  
  PercSurrDisagree = PercSurrDisagree/100
  NumSurrDisagree = int(UseSurrsN * PercSurrDisagree)
  PercGenoHaploDisagree = PercGenoHaploDisagree/100
  GenotypeMissingErrorPercentage = GenotypeMissingErrorPercentage/100


  if (Simulation == 1) then
    open (unit = 16, file = trim(TruePhaseFile), status = "old")
  end if

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

subroutine calculateCores(CoreIndex, TailIndex)
  use Global
  implicit none
  
  integer, dimension(:,:), allocatable :: CoreIndex, TailIndex
  
  integer :: resid
  double precision :: corelength
  integer :: left, ltail, rtail, nCores
  integer :: i
  
  if (consistent) then
    if (Offset == 0) then
  !    StartCoreSnp = 1
  !    EndCoreSnp = CoreAndTailLength
  !    nSnpErrorThresh = int(GenotypeMissingErrorPercentage * CoreAndTailLength)
      !NumSurrDisagree = int(UseSurrsN * PercSurrDisagree)

      nCores = int(nSnp)/Jump
      allocate(CoreIndex(nCores, 2))
      allocate(TailIndex(nCores, 2))

      resid = int((CoreAndTailLength - Jump)/2)
      CoreIndex(1, 1) = 1
      CoreIndex(1, 2) = 1 + Jump - 1
      TailIndex(1, 1) = 1
      TailIndex(1, 2) = 1 + CoreAndTailLength - 1
      do i = 2, nCores
	CoreIndex(i, 1) = CoreIndex(i - 1, 1) + Jump
	CoreIndex(i, 2) = CoreIndex(i - 1, 2) + Jump
	TailIndex(i, 1) = CoreIndex(i, 1) - resid
	TailIndex(i, 2) = CoreIndex(i, 2) + resid
	if (TailIndex(i, 1) < 1) TailIndex(i, 1) = 1
	if (TailIndex(i, 2) > nSnp) TailIndex(i, 2) = nSnp
      end do
      CoreIndex(nCores, 2) = nSnp
    endif

    if (Offset == 1) then
  !    nSnpErrorThresh = int(GenotypeMissingErrorPercentage * CoreAndTailLength)
      !NumSurrDisagree = int(UseSurrsN * PercSurrDisagree)
      resid = int((CoreAndTailLength - Jump)/2)

      nCores = (int(nSnp)/Jump) + 1
      allocate(CoreIndex(nCores, 2))
      allocate(TailIndex(nCores, 2))

      CoreIndex(1, 1) = 1
      CoreIndex(1, 2) = int(Jump/2)
      TailIndex(1, 1) = 1
      TailIndex(1, 2) = nSnp
      do i = 2, nCores
	CoreIndex(i, 1) = CoreIndex(i - 1, 2) + 1
	CoreIndex(i, 2) = CoreIndex(i - 1, 2) + Jump
	TailIndex(i, 1) = CoreIndex(i, 1) - resid
	TailIndex(i, 2) = CoreIndex(i, 2) + resid
	if (TailIndex(i, 1) < 1) TailIndex(i, 1) = 1
	if (TailIndex(i, 2) > nSnp) TailIndex(i, 2) = nSnp
      end do
      CoreIndex(nCores, 2) = nSnp
      TailIndex(nCores, 1) = 1   
      TailIndex(nCores, 2) = nSnp
    endif
  else
    nCores = nSnp / Jump
    corelength = nSnp / nCores
    left = nSnp - nCores * corelength
    ltail = floor(dble(CoreAndTailLength - Jump) / 2.0)
    rtail = ceiling(dble(CoreAndTailLength - Jump) / 2.0)
    
    if (Offset == 0) then
      allocate(CoreIndex(nCores, 2))
      allocate(TailIndex(nCores, 2))
      CoreIndex(1, 1) = 1
      if (left /= 0) then
	CoreIndex(1, 2) = 1 + corelength
      else
	CoreIndex(1, 2) = corelength
      end if
    else
      nCores = nCores + 1
      allocate(CoreIndex(nCores, 2))
      allocate(TailIndex(nCores, 2))
      CoreIndex(1, 1) = 1
      if (left /= 0) then
	CoreIndex(1, 2) = 1 + floor(dble(corelength) / 2.0)
      else
	CoreIndex(1, 2) = floor(dble(corelength) / 2.0)
      end if
    end if
      
    do i = 2, nCores
      CoreIndex(i,1) = CoreIndex(i - 1, 2) + 1
      if (i < left) then
	CoreIndex(i, 2) = CoreIndex(i - 1, 2) + corelength + 1
      else
	CoreIndex(i, 2) = CoreIndex(i - 1, 2) + corelength
      end if
    end do
    
    if (Offset /= 0) then
      CoreIndex(nCores,2) = nSnp
    end if
    
    do i = 1, nCores
      TailIndex(i,1) = max(1,CoreIndex(i,1) - ltail)
      TailIndex(i,2) = min(nSnp,CoreIndex(i,2) + rtail)
    end do

  endif
end subroutine CalculateCores

!########################################################################################################################################################################################################

subroutine MakeDirectories
  use global
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
    open (unit = 29, file = ".\PhasingResults\PhasingYield.txt", status = "unknown")
  else
    call system("mkdir PhasingResults/HaplotypeLibrary")
    if (FullFileOutput == 1) call system("mkdir PhasingResults/HaplotypeLibrary/Extras")
    open (unit = 29, file = "./PhasingResults/PhasingYield.txt", status = "unknown")
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

end subroutine MakeDirectories

!########################################################################################################################################################################

!subroutine AllocateGlobalArrays
!  use Global
!  implicit none

  !allocate(Surrogates(nAnisG, nAnisG, 3))
!  allocate(FullyPhased(nAnisG, 2))
!  allocate(HapFreq(nAnisG * 2))
!  allocate(HapAnis(nAnisG, 2))
!  allocate(AllHapAnis(nAnisG, 2, nCores))
!  AllHapAnis = -99

!  allocate(nSnpErrorThreshAnims(nAnisG * (nAnisG + 1)/2))

!end subroutine AllocateGlobalArrays


!########################################################################################################################################################################
!pure function GetnSnpErrorThreshAnims(i, j)
!  use Global
!  implicit none
!
!  integer, intent(in) :: i, j
!  integer :: k
!  integer :: GetnSnpErrorThreshAnims
!
!  if (j > i) then
!    k = (j - 1) * j/2 + i
!  else
!    k = (i - 1) * i/2 + j
!  endif
!
!  GetnSnpErrorThreshAnims = nSnpErrorThreshAnims(k)
!
!end function GetnSnpErrorThreshAnims

!######################################################################################################################################################

!#################################################################################################################################################################


subroutine Flipper(phase,startCoreSnp,endCoreSnp)
  !use Global
  implicit none

  integer(kind=1), dimension(:,:,:), intent(inout) :: phase
  integer, intent(in) :: startCoreSnp, endCoreSnp
  
  integer :: i, j, SizeCore, CountAgreeStay1, CountAgreeStay2, CountAgreeSwitch1, CountAgreeSwitch2, truth, TempVal
  integer(kind = 1), allocatable, dimension(:,:,:) :: TruePhase
  integer(kind = 1), allocatable, dimension(:) :: W1, W2
  character(len = 300) :: dumC
  
  integer :: nAnisG, nSnp

  nAnisG = size(phase,1)
  nSnp = size(phase,2)
  
  TempVal = (EndCoreSnp - StartCoreSnp) + 1

  allocate(TruePhase(nAnisG, nSnp, 2))
  allocate(W1(TempVal))
  allocate(W2(TempVal))

  SizeCore = (EndCoreSnp - StartCoreSnp) + 1

  rewind (16)
  do i = 1, nAnisG
    read (16, *) dumC, TruePhase(i,:, 1)
    read (16, *) dumC, TruePhase(i,:, 2)
  end do
  rewind (16)

  do i = 1, nAnisG
    CountAgreeStay1 = 0
    CountAgreeStay2 = 0
    CountAgreeSwitch1 = 0
    CountAgreeSwitch2 = 0
    truth = 0
    do j = StartCoreSnp, EndCoreSnp
      if (TruePhase(i, j, 1) == Phase(i, j, 1)) CountAgreeStay1 = CountAgreeStay1 + 1
      if (TruePhase(i, j, 2) == Phase(i, j, 1)) CountAgreeSwitch1 = CountAgreeSwitch1 + 1
      if (TruePhase(i, j, 1) == Phase(i, j, 2)) CountAgreeSwitch2 = CountAgreeSwitch2 + 1
      if (TruePhase(i, j, 2) == Phase(i, j, 2)) CountAgreeStay2 = CountAgreeStay2 + 1
    end do
    if ((CountAgreeSwitch2 > CountAgreeStay2).and.(CountAgreeStay1 <= CountAgreeSwitch1)) truth = 1
    if ((CountAgreeSwitch1 > CountAgreeStay1).and.(CountAgreeStay2 <= CountAgreeSwitch2)) truth = 1
    if (truth == 1) then
      W1(:) = Phase(i, StartCoreSnp:EndCoreSnp, 1)
      W2(:) = Phase(i, StartCoreSnp:EndCoreSnp, 2)
      Phase(i, StartCoreSnp:EndCoreSnp, 1) = W2(:)
      Phase(i, StartCoreSnp:EndCoreSnp, 2) = W1(:)
    end if
  end do

end subroutine Flipper

!#################################################################################################################################################################


!#################################################################################################################################################################

subroutine HapCommonality(library, OutputPoint)
  use Global, only: FullFileOutput, WindowsLinux
  use HaplotypeLibrary
  implicit none
  
  type(HapLib), intent(in) :: library
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

!#################################################################################################################################################################

!subroutine RationaliseLibrary
!  use Global
!  implicit none
!
!  integer, parameter :: HapFreqThresh = 1
!  double precision, parameter :: HapAgreeThresh = 0.98
!
!  integer :: i, j, CountCompatHap, CompatHapId, SizeCore, Counter
!  SizeCore = (EndCoreSnp - StartCoreSnp) + 1
!
!  Counter = 0
!  do i = nHaps, 1, -1
!    CountCompatHap = 0
!    CompatHapId = 0
!    if (HapFreq(i) <= HapFreqThresh) then
!      CountCompatHap = 0
!      do j = 1, nHaps
!	if ((float(HapRel(i, j))/SizeCore) > HapAgreeThresh) then
!	  CountCompatHap = CountCompatHap + 1
!	  CompatHapId = j
!	end if
!      end do
!    end if
!    if (CountCompatHap > 0) then
!      print*, i, CompatHapId
!      Counter = Counter + 1
!
!    end if
!  end do
!
!  print*, Counter
!
!end subroutine RationaliseLibrary

!############################################################################################################################################################################################################################

!######################################################################################################################################################


!######################################################################################################################################################

!subroutine EvaluateMedoidsHaps
!  use GlobalClusteringHaps
!  implicit none
!
!  integer :: i, j, k, l
!  
!  ! Calculates the middle of each mediod (average allele?)
!
!  Medoids = 0
!  CountCluster = 0
!  do i = 1, nClusters
!      do k = 1, nHapsCluster
!	if (ClusterMember(k) == i) then
!	  do l = 1, SnpInCore
!	    Medoids(i, l) = Medoids(i, l) + TempHapArray(k, l)
!	  end do
!	  CountCluster(i) = CountCluster(i) + 1
!	end if
!      end do	
!    Medoids(i,:) = Medoids(i,:)/CountCluster(i)
!  end do
!
!end subroutine EvaluateMedoidsHaps

!######################################################################################################################################################

!subroutine RePartitionHaps
!  use GlobalClusteringHaps
!  implicit none
!
!  integer :: i, j, k, l
!
!  do i = 1, nHapsCluster
!    Dist = 0
!    do j = 1, SnpInCore
!      do k = 1, nClusters
!	Dist(k) = Dist(k) + abs(TempHapArray(i, j) - Medoids(k, j))
!      end do
!    end do
!    Dist = Dist/SnpInCore
!    ! Distance is average distance to each mediod for i
!    MinClust(i) = Dist(ClusterMember(i))
!    ! Minclust is the distance to it's current mediod for i
!    do k = 1, nClusters
!      ! if distance to another mediod is less or equal than the distance to the current mediod then swap hap to that mediod and stop
!      ! repartitioning
!      if ((Dist(k) <= MinClust(i)).and.(ClusterMember(i) /= k)) then
!	MinClust(i) = Dist(k)
!	ClusterMember(i) = k
!	!print *, "Swapped ", i, "to", k
!	Change = 1
!	goto 11
!      end if
!    end do
!  end do
!
!  11 rounds = rounds + 1
!
!end subroutine RePartitionHaps

!########################################################################################################################################################################


!########################################################################################################################################################################

subroutine CountInData(nAnisRawPedigree, nAnisG)
  use Global, only: PedigreeFile, GenotypeFile
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
  
  open (unit = 3, file = trim(GenotypeFile), status = "old")
  do
    read (3, *, iostat = k) dumC
    nAnisG = nAnisG + 1
    if (k /= 0) then
      nAnisG = nAnisG - 1
      exit
    endif
  enddo
  close(3)
  print*, " ", nAnisG, " individuals in the genotype file"

  if (trim(PedigreeFile) == "NoPedigree") nAnisRawPedigree = nAnisG

end subroutine CountInData

!####################################################################################################################################################################

!########################################################################################################################################################################################################

!subroutine MarkerNRMMaker
!  use Global
!  use GlobalPedigree
!  implicit none
!
!  integer :: i, j, k, nMissing
!  double precision :: Sumpq
!  double precision, allocatable, dimension (:,:) :: GenosR, tpose
!
!  allocate(GenosR(nAnisG, nSnp))
!  !allocate(tpose(nSnp,nAnisG))
!  allocate(AlleleFreq(nSnp))
!  !allocate(MarkerNRM(nAnisG, nAnisG))
!
!  Sumpq = 0.000001
!  allelefreq = 0.000001
!  do i = 1, nSnp
!    nMissing = 0
!    do j = 1, nAnisG
!      if (Genos(j, i) /= MissingGenotypeCode) then
!	AlleleFreq(i) = AlleleFreq(i) + Genos(j, i)
!      else
!	nMissing = nMissing + 1
!      end if
!    end do
!    AlleleFreq(i) = AlleleFreq(i)/(2 * (nAnisG - nMissing))
!    Sumpq = Sumpq + (AlleleFreq(i)*(1 - AlleleFreq(i)))
!    write (11, '(i10,f7.4)') i, AlleleFreq(i)
!  end do
!  print*, " "
!  print*, " Making marker derived NRM"
!
!  do i = 1, nAnisG
!    do j = 1, nSnp
!      !if (Genos(i,j)/=MissingGenotypeCode) then
!      !        GenosR(i,j)=float((Genos(i,j)-1))-(2.0*(AlleleFreq(j)-0.5))
!      !else
!      !        GenosR(i,j)=(AlleleFreq(j)-1.0)-(2.0*(AlleleFreq(j)-0.5))
!      !end if
!    end do
!  end do
!
!  !tpose=transpose(GenosR)
!  !MarkerNRM=matmul(GenosR,tpose)
!  deallocate(AlleleFreq)
!
!  do i = 1, nAnisG
!    do j = 1, nAnisG
!      !MarkerNRM(i,j)=MarkerNRM(i,j)/(2*Sumpq)
!    end do
!  end do
!  MarkerNRM = 0.0001
!
!  do i = 1, nAnisG
!    write (12, '(a20,20000f6.2,20000f6.2,20000f6.2,20000f6.2,20000f6.2)') GenotypeId(i), MarkerNRM(i,:)
!  end do
!
!end subroutine MarkerNRMMaker

!########################################################################################################################################################################################################



!########################################################################################################################################################################

subroutine PedigreeViewerRecode(nobs, nAnisPedigree, ped, id)
  use GlobalPedigree
  implicit none
  
  integer :: nobs, nAnisPedigree
  character(len = lengan), dimension(:,:), intent(in) :: ped
  character(lengan), allocatable :: Id(:)

  integer, allocatable :: passedorder(:)
  character (len = lengan), allocatable :: holdid(:), holdsireid(:), holddamid(:)
  character (len = lengan), allocatable :: Sortedid(:), Sortedsire(:), Sorteddam(:)
  character(lengan), allocatable :: sire(:), dam(:)
  character (len = lengan) :: IDhold
  integer, allocatable :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
  integer, allocatable :: OldN(:), NewN(:), holdsire(:), holddam(:)
  integer :: mode ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
  integer :: i, j, k, kk, newid, itth, itho, ihun, iten, iunit
  integer :: nsires, ndams, newsires, newdams, nbisexuals, flag
  integer :: ns, nd, iextra, oldnobs, kn, kb, oldkn, ks, kd
  integer :: Noffset, Limit, Switch, ihold, ipoint

  mode = 1
!  allocate(id(0:nobs), sire(nobs), dam(nobs), seqid(nobs), seqsire(nobs), seqdam(nobs))
  allocate(id(0:nobs), sire(nobs), dam(nobs), seqsire(nobs), seqdam(nobs))

  do i = 1, nobs
    id(i) = ped(i, 1)
    sire(i) = ped(i, 2)
    dam(i) = ped(i, 3)
  end do

  do j = 1, nobs
    if (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') then
      dam(j) = '0'
      seqdam(j) = 0
    endif
    if (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') then
      sire(j) = '0'
      seqsire(j) = 0
    endif
  enddo
  if (mode .eq. 1) then
    !print*,  ' Inserting dummy IDs into Pedigree ... '
    newid = 0
    do j = 1, nobs
      if (((sire(j) == '0').and.(dam(j) .ne. '0')) .or. ((sire(j) .ne. '0').and.(dam(j) == '0'))) then
	newid = newid + 1
	if (newid .gt. 99999) then
	  PRINT*, newid, ' ...'
	  stop 'too many dummy single parent IDs'
	endif
	itth = int(newid/10000)
	itho = int(newid/1000) - 10 * itth
	ihun = int(newid/100) - 10 * itho - 100 * itth
	iten = int(newid/10) - 10 * ihun - 100 * itho - 1000 * itth
	iunit = newid - 10 * iten - 100 * ihun - 1000 * itho - 10000 * itth
	if (sire(j) == '0') sire(j) = 'dum' // achar(48 + itth) // achar(48 + itho) // achar(48 + ihun) // achar(48 + iten) // achar(48 + iunit)
	if (dam(j) == '0') dam(j) = 'dum' // achar(48 + itth) // achar(48 + itho) // achar(48 + ihun) // achar(48 + iten) // achar(48 + iunit)
      endif
    enddo
  endif
  !print*,  ' Sorting Sires ... '
  allocate(Sortedid(nobs), SortedIdIndex(nobs))
  Sortedid(1:nobs) = sire(1:nobs)
  Noffset = int(nobs/2)
  do while (Noffset > 0)
    Limit = nobs - Noffset
    switch = 1
    do while (Switch .ne. 0)
      Switch = 0
      do i = 1, Limit
	if (Sortedid(i) .gt. Sortedid(i + Noffset)) then
	  IDhold = Sortedid(i)
	  Sortedid(i) = Sortedid(i + Noffset)
	  Sortedid(i + Noffset) = IDhold
	  Switch = i
	endif
      enddo
      Limit = Switch - Noffset
    enddo
    Noffset = int(Noffset/2)
  enddo
  nsires = 0
  if (Sortedid(1) /= '0') nsires = 1
  do i = 2, nobs
    if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') nsires = nsires + 1
  end do
  allocate(Sortedsire(0:nsires), SortedSireIndex(nsires))
  Sortedsire(0) = '0'
  nsires = 0
  if (Sortedid(1) /= '0') then
    nsires = 1
    Sortedsire(1) = Sortedid(1)
  endif
  do i = 2, nobs
    if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') then
      nsires = nsires + 1
      Sortedsire(nsires) = Sortedid(i)
    endif
  end do
  !print*,  ' Sorting Dams ... '
  Sortedid(1:nobs) = dam(1:nobs)
  Noffset = int(nobs/2)
  do while (Noffset > 0)
    Limit = nobs - Noffset
    switch = 1
    do while (Switch .ne. 0)
      Switch = 0
      do i = 1, Limit
	if (Sortedid(i) .gt. Sortedid(i + Noffset)) then
	  IDhold = Sortedid(i)
	  Sortedid(i) = Sortedid(i + Noffset)
	  Sortedid(i + Noffset) = IDhold
	  Switch = i
	endif
      enddo
      Limit = Switch - Noffset
    enddo
    Noffset = int(Noffset/2)
  enddo
  nDams = 0
  if (Sortedid(1) /= '0') nDams = 1
  do i = 2, nobs
    if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') nDams = nDams + 1
  end do
  allocate(Sorteddam(0:nDams), SortedDamIndex(ndams))
  Sorteddam(0) = '0'
  nDams = 0
  if (Sortedid(1) /= '0') then
    nDams = 1
    Sorteddam(1) = Sortedid(1)
  endif
  do i = 2, nobs
    if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') then
      nDams = nDams + 1
      Sorteddam(nDams) = Sortedid(i)
    endif
  end do
  !print*,  ' Sorting IDs ... '
  Sortedid(1:nobs) = id(1:nobs)
  do i = 1, nobs
    SortedIdIndex(i) = i
  end do
  Noffset = int(nobs/2)
  do while (Noffset > 0)
    Limit = nobs - Noffset
    switch = 1
    do while (Switch .ne. 0)
      Switch = 0
      do i = 1, Limit
	if (Sortedid(i) .gt. Sortedid(i + Noffset)) then
	  IDhold = Sortedid(i)
	  Sortedid(i) = Sortedid(i + Noffset)
	  Sortedid(i + Noffset) = IDhold
	  ihold = SortedIdIndex(i)
	  SortedIdIndex(i) = SortedIdIndex(i + Noffset)
	  SortedIdIndex(i + Noffset) = ihold
	  Switch = i
	endif
      enddo
      Limit = Switch - Noffset
    enddo
    Noffset = int(Noffset/2)
  enddo
  !print*,  ' Check for duplicate IDs ... '
  flag = -1
  do i = 2, nobs
    if (Sortedid(i) == Sortedid(i - 1)) then
      if (flag == -1) then
	open (4, file = 'ID_err.txt', status = 'unknown')
	write(4, *) 'Duplicated IDs ...'
	flag = 0
      end if
      write(4, *) Sortedid(i)
      flag = flag + 1
    end If
  enddo
  if (flag > -1) then
    close (4)
    print*, flag, ' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
    !pause
  end if
  !print*,  ' Males ... '
  !print*,  '  Find or set sire indices ... '
  newsires = 0
  do j = 1, nsires
    !check if already listed as an individual
    ipoint = int(nobs/2)
    Noffset = int(ipoint/2)
    do while (Noffset > 1)
      if (Sortedsire(j) .lt. Sortedid(ipoint)) then
	ipoint = ipoint - Noffset
	Noffset = int(Noffset/2)
      else
	ipoint = ipoint + Noffset
	Noffset = int(Noffset/2)
      endif
    enddo
    kn = 0
    if (Sortedsire(j) == Sortedid(ipoint)) kn = 1
    do while (ipoint < nobs .and. kn == 0 .and. Sortedsire(j) > Sortedid(ipoint))
      ipoint = ipoint + 1
    enddo
    if (Sortedsire(j) == Sortedid(ipoint)) kn = 1
    do while (ipoint > 1 .and. kn == 0 .and. Sortedsire(j) < Sortedid(ipoint))
      ipoint = ipoint - 1
    enddo
    if (Sortedsire(j) == Sortedid(ipoint)) kn = 1
    if (kn == 1) then
      SortedSireIndex(j) = SortedIdIndex(ipoint)
    else ! sire is unlisted base sire
      newsires = newsires + 1
      SortedSireIndex(j) = nobs + newsires ! for now
    endif
  end do
  allocate(holdsireid(newsires))
  kn = 0
  do j = 1, nsires
    if (SortedSireIndex(j) > nobs) then
      kn = kn + 1
      holdsireid(SortedSireIndex(j) - nobs) = Sortedsire(j)
    end if
  enddo
  if (kn /= newsires) stop'newsires error'
  !print*,  '  Find seqsire ... '
  do j = 1, nobs
    If (sire(j) == '0') Then
      seqsire(j) = 0
    else
      ipoint = int(nsires/2)
      Noffset = int(ipoint/2)
      do while (Noffset > 1)
	if (sire(j) .lt. Sortedsire(ipoint)) then
	  ipoint = ipoint - Noffset
	  Noffset = int(Noffset/2)
	else
	  ipoint = ipoint + Noffset
	  Noffset = int(Noffset/2)
	endif
      enddo
      kn = 0
      if (sire(j) == Sortedsire(ipoint)) kn = 1
      do while (ipoint < nsires .and. kn == 0 .and. sire(j) > Sortedsire(ipoint))
	ipoint = ipoint + 1
      enddo
      if (sire(j) == Sortedsire(ipoint)) kn = 1
      do while (ipoint > 1 .and. kn == 0 .and. sire(j) < Sortedsire(ipoint))
	ipoint = ipoint - 1
      enddo
      if (sire(j) == Sortedsire(ipoint)) kn = 1
      if (kn == 1) then
      seqsire(j) = SortedSireIndex(ipoint)
    else
      print*, ' Error: Sire missing: ', sire(j)
      stop
      endif
    endif
  enddo
  !print*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'
  !print*,  ' Females ... '
  !print*,  '  Find or set dam indices ... '
  newdams = 0
  nbisexuals = 0
  do j = 1, ndams
    !check if already listed as an individual
    ipoint = int(nobs/2)
    Noffset = int(ipoint/2)
    do while (Noffset > 1)
      if (Sorteddam(j) .lt. Sortedid(ipoint)) then
	ipoint = ipoint - Noffset
	Noffset = int(Noffset/2)
      else
	ipoint = ipoint + Noffset
	Noffset = int(Noffset/2)
      endif
    enddo
    kn = 0
    if (Sorteddam(j) == Sortedid(ipoint)) kn = ipoint ! store ipoint here as ipoint can change with bisexuals
    do while (ipoint < nobs .and. kn == 0 .and. Sorteddam(j) > Sortedid(ipoint))
      ipoint = ipoint + 1
    enddo
    if (Sorteddam(j) == Sortedid(ipoint)) kn = ipoint
    do while (ipoint > 1 .and. kn == 0 .and. Sorteddam(j) < Sortedid(ipoint))
      ipoint = ipoint - 1
    enddo
    if (Sorteddam(j) == Sortedid(ipoint)) kn = ipoint
    !check if already listed as a sire (and therefore bisexual)
    ipoint = int(nsires/2)
    Noffset = int(ipoint/2)
    do while (Noffset > 1)
      if (Sorteddam(j) .lt. Sortedsire(ipoint)) then
	ipoint = ipoint - Noffset
	Noffset = int(Noffset/2)
      else
	ipoint = ipoint + Noffset
	Noffset = int(Noffset/2)
      endif
    enddo
    kb = 0
    if (Sorteddam(j) == Sortedsire(ipoint)) kb = 1
    do while (ipoint < nsires .and. kb == 0 .and. Sorteddam(j) > Sortedsire(ipoint))
      ipoint = ipoint + 1
    enddo
    if (Sorteddam(j) == Sortedsire(ipoint)) kb = 1
    do while (ipoint > 1 .and. kb == 0 .and. Sorteddam(j) < Sortedsire(ipoint))
      ipoint = ipoint - 1
    enddo
    if (Sorteddam(j) == Sortedsire(ipoint)) kb = 1
    if (kb == 1) then
      nbisexuals = nbisexuals + 1
      open (5, file = 'bisex.txt', position = 'append')
      write(5, *) Sorteddam(j)
      close(5)
    endif
    if (kb == 1) then
      SorteddamIndex(j) = SortedSireIndex(ipoint)
    elseif (kn >= 1) then
      SorteddamIndex(j) = SortedIdIndex(kn)
    else ! dam is unlisted base dam
      newdams = newdams + 1
      SorteddamIndex(j) = nobs + newsires + newdams ! for now
    endif
  end do
  if (nbisexuals > 0) then
    print*, nbisexuals, ' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'
    !pause
  endif
  allocate(holddamid(newdams))
  kn = 0
  do j = 1, ndams
    if (SortedDamIndex(j) > nobs + newsires) then
      kn = kn + 1
      holddamid(SortedDamIndex(j) - nobs - newsires) = Sorteddam(j)
    end if
  enddo
  if (kn /= newdams) stop'newdams error'
  !print*,  '  Find seqdam ... '
  do j = 1, nobs
    if (dam(j) == '0') then
      seqdam(j) = 0
    else
      ipoint = int(ndams/2)
      Noffset = int(ipoint/2)
      do while (Noffset > 1)
	if (dam(j) .lt. Sorteddam(ipoint))then
	  ipoint = ipoint - Noffset
	  Noffset = int(Noffset/2)
	else
	  ipoint = ipoint + Noffset
	  Noffset = int(Noffset/2)
	endif
      enddo
      kn = 0
      if (dam(j) == Sorteddam(ipoint)) kn = 1
      do while (ipoint < ndams .and. kn == 0 .and. dam(j) > Sorteddam(ipoint))
	ipoint = ipoint + 1
      enddo
      if (dam(j) == Sorteddam(ipoint)) kn = 1
      do while (ipoint > 1 .and. kn == 0 .and. dam(j) < Sorteddam(ipoint))
	ipoint = ipoint - 1
      enddo
      if (dam(j) == Sorteddam(ipoint)) kn = 1
      if (kn == 1) then
      seqdam(j) = SorteddamIndex(ipoint)
    else
      print*, ' Error: dam missing: ', dam(j)
      stop
      endif

    endif

  enddo
  !print*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'
  !print*,  ' Arranging unlisted base parents ... '
  iextra = newsires + newdams
  if (iextra > 0) then
    !print*, ' ', iextra, ' unlisted base parents found.'
    ! SortedId and SortedIdIndex just used as a holder while redimensioning
    Sortedid(1:nobs) = id(1:nobs)
    deallocate (id)
    allocate(id(nobs + iextra))
    id(1 + iextra:nobs + iextra) = Sortedid(1:nobs)

    Sortedid(1:nobs) = sire(1:nobs)
    deallocate (sire)
    allocate(sire(nobs + iextra))
    sire(1 + iextra:nobs + iextra) = Sortedid(1:nobs)

    Sortedid(1:nobs) = dam(1:nobs)
    deallocate (dam)
    allocate(dam(nobs + iextra))
    dam(1 + iextra:nobs + iextra) = Sortedid(1:nobs)

    SortedIdIndex(1:nobs) = seqsire(1:nobs)
    deallocate (seqsire)
    allocate(seqsire(nobs + iextra))
    seqsire(1 + iextra:nobs + iextra) = SortedIdIndex(1:nobs)

    SortedIdIndex(1:nobs) = seqdam(1:nobs)
    deallocate (seqdam)
    allocate(seqdam(nobs + iextra))
    seqdam(1 + iextra:nobs + iextra) = SortedIdIndex(1:nobs)

!    deallocate (seqid)
!    allocate(seqid(nobs + iextra))

  endif

  !print*, ' Inserting unlisted base parents ...'

  oldnobs = nobs
  nobs = nobs + iextra
  !print*, ' Total number of animals = ',nobs

  allocate (passedorder(nobs))
  passedorder = 0

  do i = 1 + iextra, nobs
    passedorder(i) = i - iextra

    if (sire(i) == '0')then
      seqsire(i) = 0
    else
      seqsire(i) = iextra + seqsire(i)
      If (seqsire(i) > nobs) seqsire(i) = seqsire(i) - nobs ! for unlisted sires
    end If

    if (dam(i) == '0') then
      seqdam(i) = 0
    else
      seqdam(i) = iextra + seqdam(i)
      if (seqdam(i) > nobs) seqdam(i) = seqdam(i) - nobs
    end if
  enddo

  do i = 1, newsires
    id(i) = holdsireid(i)
    passedorder(i) = 0
    seqsire(i) = 0
    seqdam(i) = 0
  enddo !i
  do i = newsires + 1, newsires + newdams
    id(i) = holddamid(i - newsires)
    passedorder(i) = 0
    seqsire(i) = 0
    seqdam(i) = 0
  enddo !i
  deallocate(holdsireid, holddamid, SortedIdIndex, SortedId)
  flag = 0
  do i = 1, nobs
    if (i <= seqsire(i) .Or. i <= seqdam(i)) flag = 1
  enddo
  if (flag == 0) then
    goto 8012
  endif
  !print*, ' Re-Ordering pedigree ...'
  allocate ( OldN(0:nobs), NewN(0:nobs))
  allocate ( holdid(0:nobs), holdsire(nobs), holddam(nobs))
  OldN(0) = 0
  NewN = 0
  !seqsire(0) = 0 !not needed !
  !seqdam(0) = 0
  holdid(1:nobs) = id(1:nobs)
  holdsire = seqsire
  holddam = seqdam
  !Find base ancestors ...
  kn = 0
  do i = 1, nobs
    if (seqsire(i) == 0 .And. seqdam(i) == 0) then
      kn = kn + 1
      NewN(i) = kn
      OldN(kn) = i
    end if
  enddo !i
  !Re-order pedigree ...
  NewN(0) = nobs + 1
  flag = 0
  do while (kn < nobs)
    oldkn = kn
    do i = 1, nobs
      if (NewN(i) == 0) then !And id(i) <> 'UniqueNULL' Then
	Ks = seqsire(i)
	Kd = seqdam(i)
	if (NewN(Ks) > 0 .And. NewN(Kd) > 0) then
	  kn = kn + 1
	  NewN(i) = kn
	  OldN(kn) = i
	end if
      end if
    enddo !i
    ! to avoid hang on unexpected problem ...
    if (kn == oldkn) then
      flag = flag + 1
    else
      flag = 0
    endif

    if (flag > 10) then
      open(1, file = 'ped_err.txt', status = 'unknown')
      write(6, *) 'Pedigree errors found involving two or more of the following relationships ...'
      write(6, *)
      write(6, *) '       Index numbers are followed by names.'
      write(6, *) '       Index number 0 means unknown, whence name is blank.'
      write(6, *)
      do i = 1, nobs
	if (NewN(i) == 0) then
	  write(6, *) 'Individual:', i, ':  ', id(i)
	  write(6, *) '    Father:', seqsire(i), ':  ', id(seqsire(i))
	  write(6, *) '    Mother:', seqdam(i), ':  ', id(seqdam(i))
	  write(6, *)
	end if
      enddo !i
      close (6)
      print*, 'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
      stop
    endif
  enddo
  NewN(0) = 0
  do i = 1, nobs
    id(i) = holdid(OldN(i))
  enddo
  do i = 1, nobs
    seqsire(i) = NewN(holdsire(OldN(i)))
    seqdam(i) = NewN(holddam(OldN(i)))
    if (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
      print*, 'out of order'
      stop
    endif
  enddo

  deallocate ( OldN, NewN, holdid, holdsire, holddam)

  8012 do i = 1, nobs
!  seqid(i) = i
enddo
nAnisPedigree = nobs

deallocate(Sortedsire, Sorteddam, sire, dam, SortedSireIndex, SortedDamIndex)

end subroutine PedigreeViewerRecode

!####################################################################################################################################################################

subroutine Titles

  print*, ""
  print*, "                              **********************                         "
  print*, "                              *                    *                         "
  print*, "                              *   AlphaPhase 1.1   *                         "
  print*, "                              *                    *                         "
  print*, "                              **********************                         "
  print*, "                                                                              "
  print*, "                    Software For Phasing and Imputing Genotypes               "
  print*, ""
  print*, "                     Written by John Hickey and Brian Kinghorn                "
  print*, ""
  print*, ""
  print*, ""
  print*, ""

end subroutine Titles

!###################################################################################################################################################

subroutine PrintTimerTitles
  use Global

  implicit none
  real :: etime ! Declare the type of etime()
  real :: elapsed(2) ! For receiving user and system time
  real :: total, Minutes, Hours, Seconds

  print*, ""
  print*, ""
  print*, ""
  print*, "                              **********************                         "
  print*, "                              *                    *                         "
  print*, "                              *   AlphaPhase 1.1   *                         "
  print*, "                              *                    *                         "
  print*, "                              **********************                         "
  print*, "                                                                              "
  print*, "                    Software For Phasing and Imputing Genotypes               "
  print*, ""
  print*, "                     Written by John Hickey and Brian Kinghorn                "
  PRINT*, ""
  PRINT*, "                                  No Liability"
  PRINT*, "                          Bugs to John.Hickey@une.edu.au"
  PRINT*, ""
  PRINT*, "                Analysis Finished                         "

  total = etime(elapsed)
  Minutes = total/60
  Seconds = Total - (INT(Minutes) * 60)
  Hours = Minutes/60
  Minutes = INT(Minutes)-(INT(Hours) * 60)

  PRINT '(A107,A7,I3,A9,I3,A9,F6.2)', "Time Elapsed", "Hours", INT(Hours), "Minutes", INT(Minutes), "Seconds", Seconds

  if (WindowsLinux == 1) then
    open (unit = 32, file = ".\PhasingResults\Timer.txt", status = "unknown")
  else
    open (unit = 32, file = "./PhasingResults/Timer.txt", status = "unknown")
  endif

  write(32, '(A27,A7,I3,A9,I3,A9,F6.2)') "Time Elapsed", "Hours", INT(Hours), "Minutes", INT(Minutes), "Seconds", Seconds

end subroutine PrintTimerTitles

!############################################################################################################################################################################################################################

subroutine WriteSurrogates(definition, threshold, OutputPoint)
  use SurrogateDefinition
  use Global, only : FullFileOutput, WindowsLinux, GenotypeId
  
  implicit none
  
  character(len = 300) :: filout
  integer :: i, j, nSurrogates
  
  type(SurrDef), intent(in) :: definition
  integer, intent(in) :: threshold
  integer, intent(in) :: OutputPoint
  
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
	write (13, '(a20,20000i6,20000i6,20000i6,20000i6)') GenotypeId(i), definition%partition(i,:)
      else
	write (13, *) GenotypeId(i), definition%partition(i,:)
      end if
      do j = 1, nAnisG
	if (definition%numOppose(i, j) <= threshold) nSurrogates = nSurrogates + 1
      enddo
      write (19, '(a20,20000i6,20000i6,20000i6,20000i6)') &
      GenotypeId(i), count(definition%partition(i,:) == 1), count(definition%partition(i,:) == 2)&
      , count(definition%partition(i,:) == 3), nSurrogates, definition%method(i)
    enddo
  end if
    
end subroutine WriteSurrogates



subroutine RandomOrder(order, n, start, idum)
  !     Generate a random ordering of the integers 1 ... n.
  implicit none

  integer, INTENT(IN) :: n, start
  !integer, INTENT(OUT) :: order(n)
  integer, allocatable, INTENT(OUT) :: order(:)
  integer :: idum
  double precision ran1

  !     Local variables
  integer :: i, j, k
  double precision :: wk

  allocate(order(n))

  do i = 1, n
    order(i) = start - 1 + i
  end do

  !     Starting at the end, swap the current last indicator with one
  !     randomly chosen from those preceeding it.

  do i = n, 2, -1
    wk = ran1(idum)
    j = 1 + i * wk
    if (j < i) then
      k = order(i)
      order(i) = order(j)
      order(j) = k
    end if
  end do

  RETURN
end subroutine RandomOrder



FUNCTION ran1(idum)
  ! This Function returns a uniform random deviate between 0.0 and 1.0.
  ! Set IDUM to any negative value to initialize or reinitialize the sequence.
  !MODIFIED FOR REAL
  IMPLICIT NONE
  INTEGER idum, IA, IM, IQ, IR, NTAB, NDIV
  DOUBLE PRECISION ran1, AM, EPS, RNMX
  PARAMETER (IA = 16807, IM = 2147483647, AM = 1./IM, IQ = 127773, IR = 2836, NTAB = 32, NDIV = 1 + (IM - 1)/NTAB, EPS = 1.2e-7, RNMX = 1. - EPS)
  INTEGER j, k, iv(NTAB), iy
  SAVE iv, iy
  DATA iv /NTAB * 0/, iy /0/
  IF (idum .le. 0.or.iy .eq. 0) then
    idum = max(-idum, 1)
    DO 11 j = NTAB + 8, 1, -1
      k = idum/IQ
      idum = IA * (idum - k * IQ) - IR * k
      IF (idum .lt. 0) idum = idum + IM
      IF (j .le. NTAB) iv(j) = idum

      11 CONTINUE
      iy = iv(1)
    END IF
    k = idum/IQ
    idum = IA * (idum - k * IQ) - IR * k
    IF (idum .lt. 0) idum = idum + IM
    j = 1 + iy/NDIV
    iy = iv(j)
    iv(j) = idum
    ran1 = min(AM * iy, RNMX)
    RETURN
  END function ran1
