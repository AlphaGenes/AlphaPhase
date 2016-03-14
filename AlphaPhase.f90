!     Last change:  JM   21 Jun 2011   11:54 am
!v On the 6th of May 2015 Roberton Antolin and John Hickey agreed that the version Roberto has becomes the master version. The expectation is that this should include all changes made by John when visiting Holland and prior to Roberto starting to work here. John will keep his old version in case this turns out not to be true
!v Roberto will take the fastPHASE HMM and add it to this code and also do the shuffling of markers Speed up to this code. 
!v We have diffed the version of the first commit and John's vesion and they are exactly the same file (i.e., Roberto's master is 100% certain to include Johns original edits in holland etc)
!####################################################################################################################################################################

module Global
  implicit none

  integer, parameter :: WindowsLinux = 0 !If 1 then compile for Windows / If 0 then compile for Linux

  integer, parameter :: lengan = 20, MissingGenotypeCode = 3, NRMmemTemp = 20000

  integer :: nAnisG, nAnisRawPedigree, nAnisP, nCores, Offset

  integer :: nSnp
  integer :: CoreAndTailLength
  integer :: Jump
  integer :: UseSurrsN
  integer :: FullFileOutput
  integer :: Simulation
  integer :: Graphics
  integer :: GenotypeFileFormat
  integer :: nPruneIterations

  double precision :: GenotypeMissingErrorPercentage, PercSurrDisagree, PercGenoHaploDisagree
  double precision :: NrmThresh

  integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp, nSnpErrorThresh, OutputPoint, CurrentLoop, NumSurrDisagree, CurrentCore
  integer, allocatable, dimension(:) :: nSnpErrorThreshAnims

  integer :: ErdosNumber, HighestErdos
  integer(kind = 1), allocatable, dimension (:) :: Visited
  integer :: AlleleCount(2)
  integer, allocatable, dimension (:) :: SurrAveDiff

  integer(kind = 4), allocatable, dimension (:) :: SireGenotyped, DamGenotyped
  integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM, Genos
  integer(kind = 1), allocatable, dimension (:,:,:) :: Phase
  !integer(kind = 2), allocatable, dimension (:,:,:) :: Surrogates
  real(kind = 4), allocatable, dimension (:,:) :: NRM
  double precision, allocatable, dimension (:) :: AlleleFreq
  double precision, allocatable, dimension (:,:) :: MarkerNRM
  character(lengan), allocatable :: GenotypeId(:)

  !integer(kind = 1), allocatable, dimension (:,:) :: HapLib
  integer, allocatable, dimension (:,:,:) :: AllHapAnis
  integer, allocatable, dimension (:,:) :: FullyPhased, HapAnis, CoreIndex, TailIndex !, HapRel
  integer, allocatable, dimension (:) :: HapFreq
  integer :: nHaps, nGlobalHaps, nGlobalHapsIter
  character (len = 300) :: PedigreeFile

  integer :: secs

end module Global

!####################################################################################################################################################################

module GlobalClustering
  implicit none

  integer, parameter :: nCLusters = 2, nMaxRounds = 100
  integer :: Change, rounds, SurrCounter
  double precision :: CountCluster(nClusters), Dist(nClusters)
  integer, allocatable, dimension (:,:) :: TempSurrArray
  integer, allocatable, dimension (:) :: TempSurrVector, ClusterMember
  double precision, allocatable, dimension (:,:) :: Medoids
  double precision, allocatable, dimension (:) :: MinClust

end module GlobalClustering

!####################################################################################################################################################################

module GlobalClusteringHaps
  implicit none

  integer, parameter :: nCLusters = 2, nMaxRounds = 100
  integer :: Change, rounds, SurrCounter, nHapsCluster, SnpInCore
  double precision :: CountCluster(nClusters), Dist(nClusters)
  integer, allocatable, dimension (:,:) :: TempHapArray
  integer, allocatable, dimension (:) :: TempHapVector, ClusterMember
  double precision, allocatable, dimension (:,:) :: Medoids
  double precision, allocatable, dimension (:) :: MinClust

end module GlobalClusteringHaps

!####################################################################################################################################################################

module GlobalPedigree
  use Global
  implicit none

  real(kind = 4), allocatable :: xnumrelmatHold(:)
  integer :: NRMmem, shell, shellmax, shellWarning
  integer, allocatable :: seqid(:), seqsire(:), seqdam(:), RecodeGenotypeId(:), RecSire(:), RecDam(:)
  character(lengan), allocatable :: ped(:,:), Id(:)

end module GlobalPedigree

!####################################################################################################################################################################

program Rlrplhi
  use Global
  use HaplotypeLibrary
  use SurrogateDefinition
  implicit none

  integer :: h, i, j, counter, SizeCore, nGlobalHapsOld, nCount
  double precision :: value, Yield
  
  type(HapLib) :: library
  type(SurrDef) :: surrogates

  ! Create a seed for RNG
  ! call system_clock(nCount)
  ! secs = mod(nCount,int(1e6))

  call Titles
  call ReadInParameterFile
  call MakeDirectories
  call CountInData
  call ParseData
  call AllocateGlobalArrays
  Phase = 9
  nPruneIterations = 1

  do h = 1, nCores
    CurrentCore = h
    nGlobalHaps = 0
    nGlobalHapsIter = 1
    print*, " "
    print*, " "
    print*, " Starting Core", CurrentCore
    OutputPoint = CurrentCore
    StartCoreSnp = CoreIndex(h, 1)
    EndCoreSnp = CoreIndex(h, 2)
    StartSurrSnp = TailIndex(h, 1)
    EndSurrSnp = TailIndex(h, 2)
    
    !call MakeStartSurrogateArray
    !!!!!!NEED TO SORT OUT THE THRESHOLD!!!!!
    call surrogates%calculate(Genos(:,StartSurrSnp:EndSurrSnp), SireGenotyped, DamGenotyped, 0)
    call Erdos(surrogates, 0)
    call CheckCompatHapGeno
    call library%initalise(EndCoreSnp-StartCoreSnp+1,500,500)
    call MakeHapLib(library)
    nGlobalHapsOld = nGlobalHaps
    print*, " "
    print*, "  ", "Haplotype library imputation step"
    do j = 1, 20
      call ImputeFromLib(library)
      call MakeHapLib(library)
      if (nGlobalHapsOld == nGlobalHaps) exit
      nGlobalHapsOld = nGlobalHaps
    end do
    call WriteHapLib(library)
    call HapCommonality(library)

    ! HIDDEN MARKOV MODEL SHOULD COME HERE

    !call RationaliseLibrary
!    if (Simulation == 1) then
!      call Flipper
!      call Checker
!    end if
  end do

  call WriteOutResults
  call PrintTimerTitles

end program Rlrplhi

!####################################################################################################################################################################

subroutine ReadInParameterFile
  use Global
  implicit none

  integer :: i, resid
  character (len = 300) :: dumC, GenotypeFile, TruePhaseFile, FileFormat, OffsetVariable

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
  if (trim(PedigreeFile) /= "NoPedigree") open (unit = 2, file = trim(PedigreeFile), status = "old")
  open (unit = 3, file = trim(GenotypeFile), status = "old")

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

  PercSurrDisagree = PercSurrDisagree/100
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

  if (Offset == 0) then
    StartCoreSnp = 1
    EndCoreSnp = CoreAndTailLength
    nSnpErrorThresh = int(GenotypeMissingErrorPercentage * CoreAndTailLength)
    NumSurrDisagree = int(UseSurrsN * PercSurrDisagree)

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
    nSnpErrorThresh = int(GenotypeMissingErrorPercentage * CoreAndTailLength)
    NumSurrDisagree = int(UseSurrsN * PercSurrDisagree)
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

end subroutine ReadInParameterFile

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

subroutine AllocateGlobalArrays
  use Global
  implicit none

  !allocate(Surrogates(nAnisG, nAnisG, 3))
  allocate(FullyPhased(nAnisG, 2))
  allocate(HapFreq(nAnisG * 2))
  allocate(HapAnis(nAnisG, 2))
  allocate(AllHapAnis(nAnisG, 2, nCores))
  AllHapAnis = -99

  allocate(nSnpErrorThreshAnims(nAnisG * (nAnisG + 1)/2))

end subroutine AllocateGlobalArrays


!########################################################################################################################################################################
pure function GetnSnpErrorThreshAnims(i, j)
  use Global
  implicit none

  integer, intent(in) :: i, j
  integer :: k
  integer :: GetnSnpErrorThreshAnims

  if (j > i) then
    k = (j - 1) * j/2 + i
  else
    k = (i - 1) * i/2 + j
  endif

  GetnSnpErrorThreshAnims = nSnpErrorThreshAnims(k)

end function GetnSnpErrorThreshAnims

!######################################################################################################################################################

subroutine EvaluateMedoids
  use GlobalClustering
  implicit none

  integer :: i, j, k, l

  Medoids = 0
  CountCluster = 0
  do i = 1, nClusters
    do j = 1, SurrCounter
      if (ClusterMember(j) == i) then
	do k = 1, SurrCounter
	  if (ClusterMember(k) == i) then
	    do l = 1, SurrCounter
	      Medoids(i, l) = Medoids(i, l) + TempSurrArray(k, l)
	    end do
	  end if
	end do
	CountCluster(i) = CountCluster(i) + 1
      end if
    end do
    Medoids(i,:) = Medoids(i,:)/(CountCluster(i)**2)
  end do

end subroutine EvaluateMedoids

!######################################################################################################################################################

subroutine RePartition
  use GlobalClustering
  implicit none

  integer :: i, j, k, l

  do i = 1, SurrCounter
    Dist = 0
    do j = 1, SurrCounter
      do k = 1, nClusters
	Dist(k) = Dist(k) + abs(TempSurrArray(i, j) - Medoids(k, j))
      end do
    end do
    Dist = Dist/SurrCounter
    MinClust(i) = Dist(ClusterMember(i))
    do k = 1, nClusters
      if ((Dist(k) <= MinClust(i)).and.(ClusterMember(i) /= k)) then
	MinClust(i) = Dist(k)
	ClusterMember(i) = k
	Change = 1
	!goto 10    !here
      end if
    end do
  end do

  10 rounds = rounds + 1

end subroutine RePartition

!######################################################################################################################################################

subroutine Erdos(surrogates, threshold)
  use Global
  use SurrogateDefinition
  implicit none

  type(SurrDef), intent(in) :: surrogates
  integer :: threshold
  
  integer :: i, j
  integer :: counter, IterAllele, SizeCore
  double precision :: value

  integer :: GetnSnpErrorThreshAnims

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
	Phase(i, j, 1) = IterAllele(i, j, 1, surrogates, threshold)
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
	Phase(i, j, 2) = IterAllele(i, j, 2, surrogates, threshold)
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

function IterAllele(animal, snp, SideOn, surrogates, threshold) result (iAllele)
  use Global
  use SurrogateDefinition
  implicit none
  
  type(SurrDef), intent(in) :: surrogates
  integer :: threshold

  integer :: i, j, animal, snp, iAllele, SideOn
  integer(kind = 1), allocatable, dimension(:) :: ErdosNowVec, ErdosNextVec

  integer :: GetnSnpErrorThreshAnims


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

!#################################################################################################################################################################

subroutine CheckCompatHapGeno
  use Global
  implicit none

  integer :: i, j, CountError, SizeCore, ErrorAllow, Disagree, counter, counterMissing
  double precision :: value, Yield

  SizeCore = (EndCoreSnp - StartCoreSnp) + 1
  ErrorAllow = int(PercGenoHaploDisagree * SizeCore)

  do i = 1, nAnisG
    CountError = 0
    counterMissing = 0
    do j = StartCoreSnp, EndCoreSnp
      if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) /= 9)) then
	counterMissing = counterMissing + 1
	if ((Genos(i, j) /= MissingGenotypeCode).and.(sum(Phase(i, j,:)) /= Genos(i, j))) CountError = CountError + 1
      end if
    end do
    ErrorAllow = int(PercGenoHaploDisagree * counterMissing)
    if (CountError >= ErrorAllow) then
      do j = StartCoreSnp, EndCoreSnp
	if (Genos(i, j) /= MissingGenotypeCode) then
	  if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) /= 9).and.(sum(Phase(i, j,:)) /= Genos(i, j))) then
	    if (Genos(i, j) == 1) Phase(i, j, 2) = 9
	    if (Genos(i, j) == MissingGenotypeCode) Phase(i, j, 2) = 9
	    if (Genos(i, j) == 0) Phase(i, j,:) = 0
	    if (Genos(i, j) == 2) Phase(i, j,:) = 1
	  endif
	endif
      enddo
    endif
  end do

  if (CurrentLoop == nPruneIterations) then
    counter = count(Phase(:, StartCoreSnp:EndCoreSnp, 1) == 0)
    counter = count(Phase(:, StartCoreSnp:EndCoreSnp, 1) == 1) + counter
    Yield = (float(counter)/(nAnisG * SizeCore)) * 100
    print*, " "
    write (*, '(a3,f6.2,a45)') "  ", Yield, "% was the Paternal allele yield for this core"
    counter = count(Phase(:, StartCoreSnp:EndCoreSnp, 2) == 0)
    counter = count(Phase(:, StartCoreSnp:EndCoreSnp, 2) == 1) + counter
    Yield = (float(counter)/(nAnisG * SizeCore)) * 100
    write (*, '(a3,f6.2,a45)') "  ", Yield, "% was the Maternal allele yield for this core"
  end if

end subroutine CheckCompatHapGeno

!#################################################################################################################################################################

subroutine Flipper
  use Global
  implicit none

  integer :: i, j, SizeCore, CountAgreeStay1, CountAgreeStay2, CountAgreeSwitch1, CountAgreeSwitch2, truth, TempVal
  integer(kind = 1), allocatable, dimension(:,:,:) :: TruePhase
  integer(kind = 1), allocatable, dimension(:) :: W1, W2
  character(len = 300) :: dumC

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

subroutine WriteHapLib(library)
  use Global
  use HaplotypeLibrary
  implicit none
  
  type(HapLib), intent(in) :: library

  integer :: i, j, k, counter, SizeCore
  character(len = 300) :: filout
  SizeCore = (EndCoreSnp - StartCoreSnp) + 1

  if (FullFileOutput == 1) then
    if (WindowsLinux == 1) then
      write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".txt")') OutputPoint
      open (unit = 24, FILE = filout, status = 'unknown')
    else
      write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".txt")') OutputPoint
      open (unit = 24, FILE = filout, status = 'unknown')
    endif
  endif
  if (WindowsLinux == 1) then
    write (filout, '(".\PhasingResults\HaplotypeLibrary\HapLib",i0,".bin")') OutputPoint
    open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
  else
    write (filout, '("./PhasingResults/HaplotypeLibrary/HapLib",i0,".bin")') OutputPoint
    open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
  endif


  write (34) nHaps, SizeCore
  do i = 1, nHaps
    if (FullFileOutput == 1)&
    write (24, '(2i6,a2,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)') &
    i, HapFreq(i), " ", library%getHap(i)
    write (34) library%getHap(i)
  end do

  print*, "   ", "Final iteration found ", nHaps, "haplotypes"

  counter = 0
  do i = 1, nAnisG
    do j = StartCoreSnp, EndCoreSnp
      do k = 1, 2
	if (Phase(i, j, k) == 0) counter = counter + 1
	if (Phase(i, j, k) == 1) counter = counter + 1
      end do
    end do
  end do
  print*, ""
  write (*, '(a4,a30,f5.2,a1)') "   ", "Final yield for this core was ", 100 * (float(counter)/(2 * nAnisG * SizeCore)), "%"

  write (29, '(i10,f7.2)') CurrentCore, 100 * (float(counter)/(2 * nAnisG * SizeCore))

end subroutine WriteHapLib

!#################################################################################################################################################################

subroutine HapCommonality(library)
  use Global
  use HaplotypeLibrary
  implicit none
  
  type(HapLib), intent(in) :: library

  integer :: i, SizeCore
  character(len = 300) :: filout
  
  integer, allocatable, dimension (:,:) :: HapRel

  SizeCore = library%getNumSnps()
  
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

subroutine MakeHapLib(library)
  use Global
  use GlobalClusteringHaps
  use HaplotypeLibrary
  implicit none

  type(HapLib), intent(in) :: library
  
  integer :: i, j, k, l, m, truth, truth1
  integer :: nSNPcore, nCount
  integer, allocatable :: Shuffle(:)

  INTERFACE
    subroutine RandomOrder(order, n, start, idum)
      !     Generate a random ordering of the integers 1 ... n.

      integer, INTENT(IN) :: n, start
      integer, allocatable, INTENT(OUT) :: order(:)
    end subroutine RandomOrder
  END INTERFACE

  nHaps = 0
  HapFreq = 0
  FullyPhased = 0
  HapAnis = -99

  ! Create a seed for RNG
  call system_clock(nCount)
  secs = mod(nCount, int(1e6))

  ! Create random indexes
  nSNPcore = EndCoreSnp - StartCoreSnp + 1 ! Total number of markers in the core

  !THIS IS HORRIBLE!
  call library%initalise(nSNPcore,500,500)
  
  allocate(Shuffle(nSNPcore))
  call RandomOrder(Shuffle, nSNPcore, StartCoreSnp, -abs(secs))

  do i = 1, nAnisG
    !Paternal Haps
    truth = 0
    ! do j=StartCoreSnp,EndCoreSnp
    ! if (Phase(i,j,1)==9) then
    do j = 1, nSNPcore
      if (Phase(i, Shuffle(j), 1) == 9) then
	truth = 1
	exit
      endif
    enddo
    if (truth == 0) then
      if (nHaps == 0) then
	nHaps = 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	HapAnis(i, 1) = nHaps
	call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 1))
      else
	do k = 1, nHaps
	  Truth1 = 0
	  ! do j=StartCoreSnp,EndCoreSnp
	  !         if (HapLib(k,j)/=Phase(i,j,1)) then
	  do j = 1, nSNPcore
	    if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) /= Phase(i, Shuffle(j), 1)) then
	      Truth1 = 1
	      exit
	    end if
	  end do
	  if (Truth1 == 0) then
	    HapFreq(k) = HapFreq(k) + 1
	    HapAnis(i, 1) = k
	    exit
	  end if
	end do
	if (Truth1 == 1) then
	nHaps = nHaps + 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 1))
	HapAnis(i, 1) = nHaps
	end if
      end if
      FullyPhased(i, 1) = 1
    endif
    !Maternal Haps
    truth = 0
    ! do j=StartCoreSnp,EndCoreSnp
    !         if (Phase(i,j,2)==9) then
    do j = 1, nSNPcore
      if (Phase(i, Shuffle(j), 2) == 9) then
	truth = 1
	exit
      endif
    enddo
    if (truth == 0) then
      if (nHaps == 0) then
	nHaps = 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	HapAnis(i, 2) = nHaps
	call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 2))
      else
	do k = 1, nHaps
	  Truth1 = 0
	  ! do j=StartCoreSnp,EndCoreSnp
	  !         if (HapLib(k,j)/=Phase(i,j,2)) then
	  do j = 1, nSNPcore
	    if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) /= Phase(i, Shuffle(j), 2)) then
	      Truth1 = 1
	      exit
	    end if
	  end do
	  if (Truth1 == 0) then
	    HapFreq(k) = HapFreq(k) + 1
	    HapAnis(i, 2) = k
	    exit
	  end if
	end do
	if (Truth1 == 1) then
	nHaps = nHaps + 1
	HapFreq(nHaps) = HapFreq(nHaps) + 1
	call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 2))
	HapAnis(i, 2) = nHaps
	end if
	FullyPhased(i, 2) = 1
      endif
    endif
  enddo
  nGlobalHaps = nHaps

  AllHapAnis(:, 1, CurrentCore) = HapAnis(:, 1)
  AllHapAnis(:, 2, CurrentCore) = HapAnis(:, 2)

end subroutine MakeHapLib

!######################################################################################################################################################

subroutine ImputeFromLib(library)
  ! Impute the phase for gametes that are not completely phased by LRP 
  ! by matching their phased loci to haplotypes in the Haplotype Library,
  ! following strategies listed in the section Step 2e of Hickey et al 2011.

  use Global
  use GlobalClusteringHaps
  use HaplotypeLibrary
  implicit none
  
  type(HapLib), intent(in) :: library

  integer :: i, j, k, l, m, truth, truth1, HapLibIter, nHapsOld, Disagree, SizeCore, ErrorAllow, HapM, HapP, nCand, nCandPat, Miss, nHapsTmp
  integer :: CompatPairs, value, WorkScaler, CountZero, CountOne, Switch
  integer :: CountA, CountB, ErrorCountAB
  integer, allocatable, dimension(:) :: CandGenos, CandHaps, WorkVec!,ErrorAllow
  integer, allocatable, dimension(:,:) :: CandPairs
  character(len = 300) :: filout

  integer :: nSNPcore, nCount
  integer, allocatable :: Shuffle(:)


  INTERFACE
    subroutine RandomOrder(order, n, start, idum)
      !     Generate a random ordering of the integers 1 ... n.

      integer, INTENT(IN) :: n, start
      integer, allocatable, INTENT(OUT) :: order(:)
    end subroutine RandomOrder
  END INTERFACE


  allocate(CandGenos(nSnp))
  allocate(CandHaps(nAnisG * 2))
  allocate(WorkVec(nAnisG * 2))
  allocate(CandPairs(nAnisG * 2, 2))

  ! Create a seed for RNG
  call system_clock(nCount)
  secs = mod(nCount, int(1e6))

  ! Create random indexes
  nSNPcore = EndCoreSnp - StartCoreSnp + 1 ! Total number of markers in the core
  allocate(Shuffle(nSNPcore))
  call RandomOrder(Shuffle, nSNPcore, StartCoreSnp, -abs(secs))

  SizeCore = (EndCoreSnp - StartCoreSnp) + 1
  ErrorAllow = int(PercGenoHaploDisagree * SizeCore)
  SnpInCore = SizeCore
  ErrorCountAB = int(SizeCore * 0.09)

  HapLibIter = 1
  nGlobalHaps = nHaps
  if (nGlobalHapsIter == 1) then
    print*, "   ", "Iteration ", nGlobalHapsIter, "found ", nHaps, "haplotypes"
    nHapsOld = 0
  else
    nHapsOld = 0
  endif

  do while (nHapsOld /= nHaps)
    HapLibIter = HapLibIter + 1
    nGlobalHapsIter = nGlobalHapsIter + 1
    nHapsOld = nHaps
    do i = 1, nAnisG
      CandHaps = 0
      nCand = 0
      ErrorAllow = int(PercGenoHaploDisagree * count(Genos(i, StartCoreSnp:EndCoreSnp) /= MissingGenotypeCode))
      if (sum(FullyPhased(i,:)) /= 2) then

	! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): PATERNAL HAPLOTYPE
	if (FullyPhased(i, 1) == 1) then
	  CandHaps = 0
	  nCand = 0
	  truth1 = 0
	  HapM = 0
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if (Genos(i,j)/=MissingGenotypeCode) then
	      ! if (HapLib(k,j)+Phase(i,j,1)/=Genos(i,j)) then
	      if (Genos(i, Shuffle(j)) /= MissingGenotypeCode) then
		if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) + Phase(i, Shuffle(j), 1) /= Genos(i, Shuffle(j))) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    truth = 1
		    exit
		  endif
		end if
	      end if
	    end do
	    if (truth == 0) then
	      HapM = k
	      nCand = nCand + 1
	      CandHaps(nCand) = k
	    end if
	  end do

	  if (nCand > 1) then
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      Disagree = 0
	      do k = 1, nCand
		do l = k, nCand
		  ! if (HapLib(CandHaps(k),j)/=HapLib(CandHaps(l),j)) then
		  if (library%getPhase(CandHaps(k), Shuffle(j) - StartCoreSnp + 1) /= library%getPhase(CandHaps(l), Shuffle(j) - StartCoreSnp + 1)) then
		    Disagree = 1
		  end if
		end do
	      end do
	      if (Disagree == 0) then
		! Phase(i,j,2)=HapLib(CandHaps(1),j)
		Phase(i, Shuffle(j), 2) = library%getPhase(CandHaps(1), Shuffle(j) - StartCoreSnp + 1)
	      end if
	    end do
	  endif

	  if (nCand == 1) then
	    Phase(i, StartCoreSnp:EndCoreSnp, 2) = library%getHap(HapM)
	    FullyPhased(i, 2) = 1
	    HapFreq(HapM) = HapFreq(HapM) + 1
	    HapAnis(i, 2) = HapM
	  end if

	  if (nCand == 0) then
	    Miss = 0
	    do j = StartCoreSnp, EndCoreSnp
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 1)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 2) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 2) = 1
	      truth = 0
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		! if (HapLib(k,j)/=Phase(i,j,2)) then
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) /= Phase(i, Shuffle(j), 2)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 2) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 2))
		HapAnis(i, 2) = nHaps
	      end if
	    end if
	  endif
	end if

	! If one of the gametes is completely phased (Section Step 2e.i Hickey et al. 2011): MATERNAL HAPLOTYPE
	if (FullyPhased(i, 2) == 1) then
	  truth1 = 0
	  HapP = 0
	  CandHaps = 0
	  nCand = 0
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if (Genos(i,j)/=MissingGenotypeCode) then
	      !     if (HapLib(k,j)+Phase(i,j,2)/=Genos(i,j)) then
	      if (Genos(i, Shuffle(j)) /= MissingGenotypeCode) then
		if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) + Phase(i, Shuffle(j), 2) /= Genos(i, Shuffle(j))) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    truth = 1
		    exit
		  endif
		end if
	      end if
	    enddo
	    if (truth == 0) then
	      HapP = k
	      nCand = nCand + 1
	      CandHaps(nCand) = k
	    endif
	  enddo

	  if (nCand > 1) then
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      Disagree = 0
	      do k = 1, nCand
		do l = k, nCand
		  ! if (HapLib(CandHaps(k),j)/=HapLib(CandHaps(l),j)) then
		  if (library%getPhase(CandHaps(k), Shuffle(j) - StartCoreSnp + 1) /= library%getPhase(CandHaps(l), Shuffle(j) - StartCoreSnp + 1)) then
		    Disagree = 1
		  end if
		end do
	      end do
	      if (Disagree == 0) then
		! Phase(i,j,1)=HapLib(CandHaps(1),j)
		Phase(i, Shuffle(j), 1) = library%getPhase(CandHaps(1), Shuffle(j) - StartCoreSnp + 1)
	      end if
	    enddo
	  endif

	  if (nCand == 1) then
	    Phase(i, StartCoreSnp:EndCoreSnp, 1) = library%getHap(HapP)
	    FullyPhased(i, 1) = 1
	    HapFreq(HapP) = HapFreq(HapP) + 1
	    HapAnis(i, 1) = HapP
	  endif

	  if (nCand == 0) then
	    Miss = 0
	    do j = StartCoreSnp, EndCoreSnp
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 2)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 1) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 1) = 1
	      truth = 0
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		do j = 1, nSNPcore
		  ! if (HapLib(k,j)/=Phase(i,j,1)) then
		  if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) /= Phase(i, Shuffle(j), 1)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 1) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 1))
		HapAnis(i, 1) = nHaps
	      end if
	    endif
	  endif
	end if

	! If neither of the gametes is completely phased (Section Step 2e.ii Hickey et al. 2011)
	if (sum(FullyPhased(i,:)) == 0) then
	  HapP = 0
	  HapM = 0
	  CandHaps = 0
	  nCand = 0

	  ! Find candidates for paternal haplotype
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if ((Phase(i,j,1)/=9).and.(Phase(i,j,1)/=HapLib(k,j))) Disagree=Disagree+1
	      if ((Phase(i, Shuffle(j), 1) /= 9).and.(Phase(i, Shuffle(j), 1) /= library%getPhase(k, Shuffle(j) - StartCoreSnp + 1))) Disagree = Disagree + 1
	      if (Disagree > ErrorAllow) then
		truth = 1
		exit
	      end if
	    end do

	    ! If the there is no disagreement, we've found a new candidate haplotype
	    if (truth == 0) then
	      nCand = nCand + 1
	      CandHaps(nCand) = k
	    end if
	  enddo

	  ! Update the number of candidates for paternal haplotype
	  nCandPat = nCand

	  ! If only have one paternal candidate haplotype, then
	  ! the paternal haplotype is nCand
	  if (nCand == 1) HapP = CandHaps(nCand)

	  ! Find candidates for maternal haplotype
	  do k = 1, nHaps
	    truth = 0
	    Disagree = 0
	    ! do j=StartCoreSnp,EndCoreSnp
	    do j = 1, nSNPcore
	      ! if ((Phase(i,j,2)/=9).and.(Phase(i,j,2)/=HapLib(k,j))) Disagree=Disagree+1
	      if ((Phase(i, Shuffle(j), 2) /= 9).and.(Phase(i, Shuffle(j), 2) /= library%getPhase(k, Shuffle(j) - StartCoreSnp + 1))) Disagree = Disagree + 1
	      if (Disagree > ErrorAllow) then
		truth = 1
		exit
	      end if
	    end do

	    ! If the there is no disagreement, we've found a new candidate haplotype
	    if (truth == 0) then
	      truth1 = 0
	      do j = 1, nCand
		if (k == CandHaps(j)) then
		  truth1 = 1
		  exit
		end if
	      enddo
	      if (truth1 == 0) then
		nCand = nCand + 1
		CandHaps(nCand) = k
	      endif
	    endif
	  enddo

	  ! If only have one maternal candidate haplotype, then
	  ! the maternal haplotype is nCand
	  if ((nCand - nCandPat) == 1) HapM = CandHaps(nCand)

	  ! If only one maternal candidate haplotype and many paternal candidate haplotypes
	  if ((HapM > 0).AND.(HapP == 0)) then
	    truth1 = 0
	    do k = 1, nCandPat
	      Disagree = 0
	      truth = 1
	      ! do j=StartCoreSnp,EndCoreSnp
	      ! if ((Genos(i,j)/=MissingGenotypeCode).and.&
	      !         (Genos(i,j)/=(HapLib(HapM,j)+HapLib(CandHaps(k),j)))) then
	      do j = 1, nSNPcore
		if ((Genos(i, Shuffle(j)) /= MissingGenotypeCode).and.&
		  (Genos(i, Shuffle(j)) /= (library%getPhase(HapM, Shuffle(j) - StartCoreSnp + 1) + library%getPhase(CandHaps(k), Shuffle(j) - StartCoreSnp + 1)))) then
		  Disagree = Disagree + 1
		  if (Disagree > ErrorAllow) then
		    truth = 0
		    exit
		  end if
		endif
	      enddo
	      if (truth == 1) then
		truth1 = truth1 + 1
		HapP = CandHaps(k)
		if (truth1 > 1) then
		  HapP = 0
		end if
	      end if
	    end do
	  end if

	  ! If only one paternal candidate haplotype and many maternal candidate haplotypes 
	  if ((nCandPat == 1).and.(nCand - nCandPat > 0)) then
	    truth1 = 0
	    do k = nCandPat + 1, nCand
	      Disagree = 0
	      truth = 1
	      ! do j=StartCoreSnp,EndCoreSnp
	      ! if ((Genos(i,j)/=MissingGenotypeCode).and.(Genos(i,j)/=(HapLib(HapP,j)+HapLib(CandHaps(k),j)))) then                                                             
	      do j = 1, nSNPcore
		if ((Genos(i, Shuffle(j)) /= MissingGenotypeCode).and.(Genos(i, Shuffle(j)) /= (library%getPhase(HapP, Shuffle(j) - StartCoreSnp + 1) &
		  + library%getPhase(CandHaps(k),Shuffle(j) - StartCoreSnp + 1)))) then                                                             
		Disagree = Disagree + 1
		if (Disagree > ErrorAllow) then
		  truth = 0
		  exit
		end if
	      endif
	    enddo
	    if (truth == 0) HapM = 0
	    if (truth == 1) then
	      truth1 = truth1 + 1
	      HapM = CandHaps(k)
	      if (truth1 > 1) then
		HapM = 0
	      end if
	    end if
	  end do
	end if

	! If only have one paternal candidate haplotype
	if (HapP /= 0) then
	  Phase(i, StartCoreSnp:EndCoreSnp, 1) = library%getHap(HapP)
	  FullyPhased(i, 1) = 1

	  ! Update the Library
	  HapFreq(HapP) = HapFreq(HapP) + 1
	  HapAnis(i, 1) = HapP

	  ! If no haplotype has been found for the maternal gamete, or 
	  ! there are more than one maternal candidate haplotype
	  if (HapM == 0) then
	    Miss = 0
	    do j = StartCoreSnp, EndCoreSnp
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 1)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 2) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 2) = 1
	      truth = 0

	      ! Update (if necessary) Haplotype Library with the new maternal gamete found
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		! if (HapLib(k,j)/=Phase(i,j,2)) then
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) /= Phase(i, Shuffle(j), 2)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 2) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 2))
		HapAnis(i, 2) = nHaps
	      end if
	    end if
	  end if
	end if

	! If only have one maternal candidate haplotype 
	if (HapM /= 0) then
	  ! Phase the maternal haplotype with the only paternal candidate haplotype
	  ! NOTE: This is always necessary because
	  !   - If there is only one paternal candidate haplotype, 
	  !     we have already check they are compatible (line 3307-3316) (Step 2e.ii.A)
	  !   - If there is more than one paternal candidate haplotype, 
	  !     the paternal gamete is phased as the complementary of the maternal gamete
	  Phase(i, StartCoreSnp:EndCoreSnp, 2) = library%getHap(HapM)
	  FullyPhased(i, 2) = 1

	  ! Update the Library
	  HapFreq(HapM) = HapFreq(HapM) + 1
	  HapAnis(i, 2) = HapM

	  ! If no haplotype has been found for the paternal gamete, or 
	  ! there are more than one paternal candidate haplotype
	  if (HapP == 0) then
	    Miss = 0
	    do j = StartCoreSnp, EndCoreSnp
	      if (Genos(i, j) /= MissingGenotypeCode)then
		value = Genos(i, j) - Phase(i, j, 2)
		if ((value == 0).or.(value == 1)) then
		  Phase(i, j, 1) = value
		else
		  Miss = Miss + 1
		endif
	      else
		Miss = Miss + 1
	      endif
	    enddo
	    if (Miss == 0) then
	      FullyPhased(i, 1) = 1
	      truth = 0

	      ! Update (if necessary) Haplotype Library with the new maternal gamete found
	      do k = 1, nHaps
		Disagree = 0
		! do j=StartCoreSnp,EndCoreSnp
		! if (HapLib(k,j)/=Phase(i,j,1)) then
		do j = 1, nSNPcore
		  if (library%getPhase(k, Shuffle(j) - StartCoreSnp + 1) /= Phase(i, Shuffle(j), 1)) then
		    Disagree = 1
		    exit
		  endif
		enddo
		if (Disagree == 0) then
		  HapFreq(k) = HapFreq(k) + 1
		  HapAnis(i, 1) = k
		  truth = 1
		  exit
		end if
	      end do
	      if (truth == 0) then
		nHaps = nHaps + 1
		HapFreq(nHaps) = HapFreq(nHaps) + 1
		call library%addHap(Phase(i, StartCoreSnp:EndCoreSnp, 1))
		HapAnis(i, 1) = nHaps
	      end if
	    end if
	  end if
	end if

	! If the paternal and maternal gamete cannot be identify without ambiguity
	! (more than one or none at all)           
	if ((HapP == 0).and.(HapM == 0)) then
	  CandPairs = 0
	  CompatPairs = 0
	  do k = 1, nCand
	    do l = (k + 1), nCand
	      Disagree = 0
	      truth = 1

	      ! Check agreement between pairs
	      ! do j=StartCoreSnp,EndCoreSnp
	      ! if (Genos(i,j)/=MissingGenotypeCode)then
	      !     if ((HapLib(CandHaps(k),j)+HapLib(CandHaps(l),j))/=Genos(i,j)) then
	      do j = 1, nSNPcore
		if (Genos(i, Shuffle(j)) /= MissingGenotypeCode)then
		  if ((library%getPhase(CandHaps(k), Shuffle(j) - StartCoreSnp + 1) + library%getPhase(CandHaps(l), Shuffle(j) - StartCoreSnp + 1)) /= Genos(i, Shuffle(j))) then
		    Disagree = Disagree + 1
		    if (Disagree > ErrorAllow) then
		      truth = 0
		      exit
		    end if
		  end if
		endif
	      end do

	      ! If there is a pair that agrees, 
	      ! the haplotypes are consider to be the paternal and maternal gametes (arbitrarily)
	      if (truth == 1) then
		CompatPairs = CompatPairs + 1
		CandPairs(CompatPairs, 1) = CandHaps(k)
		CandPairs(CompatPairs, 2) = CandHaps(l)
		HapP = CandHaps(k)
		HapM = CandHaps(l)
		if ((CompatPairs * CompatPairs) > (nAnisG - 1)) exit
	      end if
	    end do
	    if ((CompatPairs * CompatPairs) > (nAnisG - 1)) exit
	  end do

	  ! If only one pair agrees...
	  if (CompatPairs == 1) then
	    ! Phase the paternal haplotype and update the library with the new frequency 
	    Phase(i, StartCoreSnp:EndCoreSnp, 1) = library%getHap(HapP)
	    FullyPhased(i, 1) = 1
	    HapFreq(HapP) = HapFreq(HapP) + 1
	    HapAnis(i, 1) = HapP
	    ! Phase the maternal haplotype and update the library with the new frequency 
	    Phase(i, StartCoreSnp:EndCoreSnp, 2) = library%getHap(HapM)
	    FullyPhased(i, 2) = 1
	    HapFreq(HapM) = HapFreq(HapM) + 1
	    HapAnis(i, 2) = HapM
	  end if

	  ! If more than one pair agrees...                
	  if ((CompatPairs > 1).and.((CompatPairs * CompatPairs) < nAnisG)) then !Note the 200 number is a fudge
	    truth = 1

	    ! Check how many paternal candidates haplotypes
	    value = CandPairs(1, 1)
	    do k = 2, CompatPairs
	      if (CandPairs(k, 1) /= value) truth = 0
	    end do
	    Switch = 0

	    ! If there is only one paternal haplotype in all the candidate pairs
	    if (truth == 1) then
	      ! Phase the paternal gamete with this haplotype
	      Phase(i, StartCoreSnp:EndCoreSnp, 1) = library%getHap(value)
	      FullyPhased(i, 1) = 1
	      HapFreq(value) = HapFreq(value) + 1
	      HapAnis(i, 1) = value

	      ! If only one haplotype is found for the paternal gamete 
	      ! and many for the maternal gamete, phase each loci only all pairs agree
	      ! (Step 2e.ii.B)
	      ! do j=StartCoreSnp,EndCoreSnp
	      do j = 1, nSNPcore
		! value=HapLib(CandPairs(1,2),j)
		value = library%getPhase(CandPairs(1, 2), Shuffle(j) - StartCoreSnp + 1)
		truth1 = 1
		do k = 2, CompatPairs
		  ! if (HapLib(CandPairs(k,2),j)/=value) then
		  if (library%getPhase(CandPairs(k, 2), Shuffle(j) - StartCoreSnp + 1) /= value) then
		    truth1 = 0
		    exit
		  end if
		end do
		if (truth1 == 1) then
		  ! Phase(i,j,2)=value
		  Phase(i, Shuffle(j), 2) = value
		end if
	      end do
	      Switch = 1
	    end if

	    ! Check how many maternal candidates haplotypes
	    truth = 1
	    value = CandPairs(1, 2)
	    do k = 2, CompatPairs
	      if (CandPairs(k, 2) /= value) truth = 0
	    end do

	    ! If there is only one maternal haplotype in all the candidate pairs
	    if (truth == 1) then
	      Phase(i, StartCoreSnp:EndCoreSnp, 2) = library%getHap(value)
	      FullyPhased(i, 2) = 1
	      HapFreq(value) = HapFreq(value) + 1
	      HapAnis(i, 2) = value

	      ! If only one haplotype is found for the paternal gamete 
	      ! and many for the maternal gamete, phase each loci only all pairs agree
	      ! (Step 2e.ii.C)
	      ! do j=StartCoreSnp,EndCoreSnp
	      do j = 1, nSNPcore
		! value=HapLib(CandPairs(1,1),j)
		value = library%getPhase(CandPairs(1, 1), Shuffle(j) - StartCoreSnp + 1)
		truth1 = 1
		do k = 2, CompatPairs
		  ! if (HapLib(CandPairs(k,1),j)/=value) then
		  if (library%getPhase(CandPairs(k, 1), Shuffle(j) - StartCoreSnp + 1) /= value) then
		    truth1 = 0
		    exit
		  end if
		end do
		if (truth1 == 1) then
		  ! Phase(i,j,1)=value
		  Phase(i, Shuffle(j), 1) = value
		endif
	      enddo
	      Switch = 1
	    endif

	    ! If proband is not completely phased and have more than one candidate 
	    ! for both paternal and maternal haplotype 
	    ! (Step 2e.iv)
	    if ((sum(FullyPhased(i,:)) /= 2).and.(Switch == 0)) then

	      ! Initialize procedure of k-medoids
	      WorkVec = 0
	      do k = 1, CompatPairs
		WorkVec(CandPairs(k, 1)) = 1
		WorkVec(CandPairs(k, 2)) = 1
	      enddo
	      WorkScaler = sum(WorkVec(:))
	      allocate(TempHapArray(WorkScaler, SizeCore))
	      allocate(TempHapVector(WorkScaler))
	      nHapsCluster = 0

	      ! Clusterize with k-medoids
	      do k = 1, nAnisG * 2
		if (WorkVec(k) == 1) then
		  nHapsCluster = nHapsCluster + 1
		  TempHapVector(nHapsCluster) = k
		  TempHapArray(nHapsCluster, 1:SnpInCore) = &
		  library%getHap(k)
		end if
	      end do
	      allocate(Medoids(nClusters, SnpInCore))
	      allocate(ClusterMember(nHapsCluster))
	      allocate(MinClust(nHapsCluster))
	      do j = 1, nHapsCluster
		if (mod(j, 2) == 0) then
		  ClusterMember(j) = 1
		else
		  ClusterMember(j) = 2
		endif
	      end do
	      call EvaluateMedoidsHaps
	      Change = 0
	      MinClust = 1
	      rounds = 1
	      call RePartitionHaps
	      do j = 1, nMaxRounds
		call EvaluateMedoidsHaps
		Change = 0
		call RePartitionHaps
		if (Change == 0) exit
	      enddo
	      if (rounds <= nMaxRounds) then
		if (count(ClusterMember(:) == 2) == 1) then
		  HapM = 0
		  do j = 1, nHapsCluster
		    if (ClusterMember(j) == 2) then
		      HapM = TempHapVector(j)
		      exit
		    endif
		  enddo
		  Phase(i, StartCoreSnp:EndCoreSnp, 2) = library%getHap(HapM)
		  FullyPhased(i, 2) = 1
		  HapFreq(HapM) = HapFreq(HapM) + 1
		  HapAnis(i, 2) = HapM
		end if
		if (count(ClusterMember(:) == 1) == 1) then
		  HapP = 0
		  do j = 1, nHapsCluster
		    if (ClusterMember(j) == 1) then
		      HapP = TempHapVector(j)
		      exit
		    endif
		  enddo
		  Phase(i, StartCoreSnp:EndCoreSnp, 1) = library%getHap(HapP)
		  FullyPhased(i, 1) = 1
		  HapFreq(HapP) = HapFreq(HapP) + 1
		  HapAnis(i, 1) = HapP
		end if
		if ((count(ClusterMember(:) == 2) > 1).and.(count(ClusterMember(:) == 2) > 1)) then
		  Phase(i, StartCoreSnp:EndCoreSnp,:) = 9
		  do j = StartCoreSnp, EndCoreSnp
		    if (Genos(i, j) == 0) Phase(i, j,:) = 0
		    if (Genos(i, j) == 2) Phase(i, j,:) = 1
		    CountZero = 0
		    CountOne = 0
		    do k = 1, nHapsCluster
		      if (ClusterMember(k) == 2) then
			if (library%getPhase(TempHapVector(k), j - StartCoreSnp + 1) == 0)&
			CountZero = CountZero + 1
			if (library%getPhase(TempHapVector(k), j - StartCoreSnp + 1) == 1)&
			CountOne = CountOne + 1
		      endif
		    end do
		    if ((CountZero == 0).and.(CountOne > 0)) Phase(i, j, 2) = 1
		    if ((CountZero > 0).and.(CountOne == 0)) Phase(i, j, 2) = 0
		  end do
		  do j = StartCoreSnp, EndCoreSnp
		    if (Genos(i, j) == 0) Phase(i, j,:) = 0
		    if (Genos(i, j) == 2) Phase(i, j,:) = 1
		    CountZero = 0
		    CountOne = 0
		    do k = 1, nHapsCluster
		      if (ClusterMember(k) == 1) then
			if (library%getPhase(TempHapVector(k), j - StartCoreSnp + 1) == 0) CountZero = CountZero + 1                                                  
			if (library%getPhase(TempHapVector(k), j - StartCoreSnp + 1) == 1) CountOne = CountOne + 1                                              
		      endif
		    end do
		    if ((CountZero == 0).and.(CountOne > 0)) Phase(i, j, 1) = 1
		    if ((CountZero > 0).and.(CountOne == 0)) Phase(i, j, 1) = 0
		  end do
		endif
	      end if
	      deallocate(ClusterMember)
	      deallocate(MinClust)
	      deallocate(Medoids)
	      deallocate(TempHapArray)
	      deallocate(TempHapVector)
	    endif
	  endif
	endif
      endif
    end if
  end do
  print*, "   ", "Iteration ", nGlobalHapsIter, "found ", nHaps, "haplotypes"
enddo

do i = 1, nAnisG
  CountA = 0
  CountB = 0
  do j = StartCoreSnp, EndCoreSnp
    if (Genos(i, j) /= MissingGenotypeCode) then
      if ((Phase(i, j, 1) /= 9).and.(Phase(i, j, 2) == 9)) then
	value = Genos(i, j) - Phase(i, j, 1)
	if ((value == 0).or.(value == 1)) then !here 7th april 2011
	  Phase(i, j, 2) = value
	else
	  CountA = CountA + 1
	endif
      endif
      if ((Phase(i, j, 2) /= 9).and.(Phase(i, j, 1) == 9)) then
	value = Genos(i, j) - Phase(i, j, 2)
	if ((value == 0).or.(value == 1)) then !here 7th april 2011
	  Phase(i, j, 1) = value
	else
	  CountB = CountB + 1
	endif
      endif
    end if
  end do
  if (CountA > ErrorCountAB) then
    Phase(i, StartCoreSnp:EndCoreSnp,:) = 9
    do j = StartCoreSnp, EndCoreSnp
      if (Genos(i, j) == 0) Phase(i, j,:) = 0
      if (Genos(i, j) == 2) Phase(i, j,:) = 1
    enddo
  endif
  if (CountB > ErrorCountAB) then
    Phase(i, StartCoreSnp:EndCoreSnp,:) = 9
    do j = StartCoreSnp, EndCoreSnp
      if (Genos(i, j) == 0) Phase(i, j,:) = 0
      if (Genos(i, j) == 2) Phase(i, j,:) = 1
    enddo
  endif
end do

do i = 1, nAnisG
  do j = StartCoreSnp, EndCoreSnp
    if (Genos(i, j) == 1) then
      if ((Phase(i, j, 1) == 9).and.(Phase(i, j, 2) /= 9)) Phase(i, j, 1) = Genos(i, j) - Phase(i, j, 2)
      if ((Phase(i, j, 2) == 9).and.(Phase(i, j, 1) /= 9)) Phase(i, j, 2) = Genos(i, j) - Phase(i, j, 1)
    endif
    if (Genos(i, j) == 0) then
      if ((Phase(i, j, 1) == 9).and.(Phase(i, j, 2) /= 9)) Phase(i, j, 1) = 0
      if ((Phase(i, j, 2) == 9).and.(Phase(i, j, 1) /= 9)) Phase(i, j, 2) = 0
    endif
    if (Genos(i, j) == 2) then
      if ((Phase(i, j, 1) == 9).and.(Phase(i, j, 2) /= 9)) Phase(i, j, 1) = 1
      if ((Phase(i, j, 2) == 9).and.(Phase(i, j, 1) /= 9)) Phase(i, j, 2) = 1
    endif
  enddo
enddo

HapFreq = 0
FullyPhased = 0
HapAnis = -99
!nGlobalHaps=nHaps

deallocate(CandGenos)
deallocate(CandHaps)
deallocate(WorkVec)
deallocate(CandPairs)

end subroutine ImputeFromLib

!######################################################################################################################################################

subroutine EvaluateMedoidsHaps
  use GlobalClusteringHaps
  implicit none

  integer :: i, j, k, l

  Medoids = 0
  CountCluster = 0
  do i = 1, nClusters
    do j = 1, nHapsCluster
      if (ClusterMember(j) == i) then
	do k = 1, nHapsCluster
	  if (ClusterMember(k) == i) then
	    do l = 1, SnpInCore
	      Medoids(i, l) = Medoids(i, l) + TempHapArray(k, l)
	    end do
	  end if
	end do
	CountCluster(i) = CountCluster(i) + 1
      end if
    end do
    Medoids(i,:) = Medoids(i,:)/(CountCluster(i)**2)
  end do

end subroutine EvaluateMedoidsHaps

!######################################################################################################################################################

subroutine RePartitionHaps
  use GlobalClusteringHaps
  implicit none

  integer :: i, j, k, l

  do i = 1, nHapsCluster
    Dist = 0
    do j = 1, SnpInCore
      do k = 1, nClusters
	Dist(k) = Dist(k) + abs(TempHapArray(i, j) - Medoids(k, j))
      end do
    end do
    Dist = Dist/SnpInCore
    MinClust(i) = Dist(ClusterMember(i))
    do k = 1, nClusters
      if ((Dist(k) <= MinClust(i)).and.(ClusterMember(i) /= k)) then
	MinClust(i) = Dist(k)
	ClusterMember(i) = k
	Change = 1
	goto 11
      end if
    end do
  end do

  11 rounds = rounds + 1

end subroutine RePartitionHaps

!########################################################################################################################################################################

subroutine WriteOutResults
  use Global
  implicit none

  integer :: i, j, k, l, counter, CounterM, CounterP
  integer, allocatable, dimension(:) :: WorkOut
  double precision, allocatable, dimension(:) :: CoreCount

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

end subroutine WriteOutResults

!########################################################################################################################################################################

subroutine CountInData
  use Global
  implicit none

  integer :: k
  character (len = 300) :: dumC


  if (trim(PedigreeFile) /= "NoPedigree") then
    do
      read (2, *, iostat = k) dumC
      nAnisRawPedigree = nAnisRawPedigree + 1
      if (k /= 0) then
	nAnisRawPedigree = nAnisRawPedigree - 1
	exit
      endif
    enddo
    rewind(2)
    print*, " ", nAnisRawPedigree, " individuals in the pedigree file"
  endif

  do
    read (3, *, iostat = k) dumC
    nAnisG = nAnisG + 1
    if (k /= 0) then
      nAnisG = nAnisG - 1
      exit
    endif
  enddo
  rewind(3)
  print*, " ", nAnisG, " individuals in the genotype file"

  if (trim(PedigreeFile) == "NoPedigree") nAnisRawPedigree = nAnisG

end subroutine CountInData

!####################################################################################################################################################################

subroutine ParseData
  use GlobalPedigree
  use Global
  implicit none

  integer :: i, j, k, SumPseudoNrmS, SumPseudoNrmD, truth, counter, CountMissingGenotype, SireGen, DamGen
  real, external :: xnumrelmat
  real(kind = 4) :: value, valueS, valueD, SumNrm, SumDiag
  integer, allocatable, dimension (:) :: GenoInPed, WorkVec, ReadingVector

  allocate(GenotypeId(nAnisG))
  allocate(GenoInPed(nAnisG))
  allocate(RecodeGenotypeId(nAnisG))
  allocate(PseudoNRM(nAnisG, nAnisG))
  allocate(Ped(nAnisRawPedigree, 3))
  allocate(Genos(nAnisG, nSnp))
  allocate(Phase(nAnisG, nSnp, 2))
  allocate(WorkVec(nSnp * 2))
  allocate(ReadingVector(nSnp))

  !allocate(HapLib(nAnisG * 2, nSnp))

  if (trim(PedigreeFile) /= "NoPedigree") then
    do i = 1, nAnisRawPedigree
      read(2, *) ped(i,:)
    enddo
  else
    do i = 1, nAnisRawPedigree
      ped(i, 2:3) = "0"
      read (3, *) ped(i, 1)
    enddo
    rewind (3)
  endif

  GenoInPed = 0

  Phase = 9

  Genos = MissingGenotypeCode
  do i = 1, nAnisG
    truth = 0
    if (GenotypeFileFormat == 1) then
      read (3, *) GenotypeId(i), ReadingVector(:)
      do j = 1, nSnp
	if ((ReadingVector(j) /= 0).and.(ReadingVector(j) /= 1).and.(ReadingVector(j) /= 2)) ReadingVector(j) = MissingGenotypeCode
	Genos(i, j) = ReadingVector(j)
	if (Genos(i, j) == 0) Phase(i, j,:) = 0
	if (Genos(i, j) == 2) Phase(i, j,:) = 1
      end do
    end if
    if (GenotypeFileFormat == 2) then
      read (3, *) GenotypeId(i), Phase(i,:, 1)
      read (3, *) GenotypeId(i), Phase(i,:, 2)
    end if
    if (GenotypeFileFormat == 3) then
      read (3, *) GenotypeId(i), WorkVec(:)
      k = 0
      do j = 1, nSnp * 2
	if (mod(j, 2) == 0) then
	  k = k + 1
	  if ((WorkVec(j - 1) == 1).and.(WorkVec(j) == 1)) Genos(i, k) = 0
	  if ((WorkVec(j - 1) == 1).and.(WorkVec(j) == 2)) Genos(i, k) = 1
	  if ((WorkVec(j - 1) == 2).and.(WorkVec(j) == 1)) Genos(i, k) = 1
	  if ((WorkVec(j - 1) == 2).and.(WorkVec(j) == 2)) Genos(i, k) = 2
	endif
      end do
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
    rewind(2)
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
  else
    do i = 1, nAnisG
      ped(i, 1) = GenotypeId(i)
      ped(i, 2) = "0"
      ped(i, 3) = "0"
    enddo

  endif
  call PedigreeViewerRecode(nAnisRawPedigree, nAnisP)
  deallocate(Ped)
  deallocate(GenoInPed)

  do i = 1, nAnisG
    do j = 1, nAnisP
      if (Id(j) == GenotypeId(i)) then
	RecodeGenotypeId(i) = seqid(j)
	exit
      end if
    end do
  end do
  close (2)
  close (3)

  allocate(RecSire(0:nAnisP))
  allocate(RecDam(0:nAnisP))
  RecSire(0) = 0
  RecDam(0) = 0
  do i = 1, nAnisP
    RecSire(i) = seqsire(i)
    RecDam(i) = seqdam(i)
  enddo
  NRMmem = NRMmemTemp
  if (NRMmem > nAnisP) NRMmem = nAnisP
  shellmax = 50000
  allocate(xnumrelmatHold(-1 * NRMmem: NRMmem * (NRMmem + 1)/2))
  xnumrelmatHold = -9.
  xnumrelmatHold(-1 * NRMmem: 0) = 0.

  if (FullFileOutput == 1) then
    allocate(NRM(nAnisG, nAnisG))
    if (WindowsLinux == 1) then
      open (unit = 8, file = ".\Miscellaneous\GenotypedNRM.txt", status = "unknown")
      open (unit = 9, file = ".\Miscellaneous\GenotypedPseudoNRM.txt", status = "unknown")
      open (unit = 10, file = ".\Miscellaneous\SummaryOfGenotypedNRM.txt", status = "unknown")
      open (unit = 11, file = ".\Miscellaneous\AlleleFrequency.txt", status = "unknown")
      open (unit = 12, file = ".\Miscellaneous\GenotypedMarkerNRM.txt", status = "unknown")
    else
      open (unit = 8, file = "./Miscellaneous/GenotypedNRM.txt", status = "unknown")
      open (unit = 9, file = "./Miscellaneous/GenotypedPseudoNRM.txt", status = "unknown")
      open (unit = 10, file = "./Miscellaneous/SummaryOfGenotypedNRM.txt", status = "unknown")
      open (unit = 11, file = "./Miscellaneous/AlleleFrequency.txt", status = "unknown")
      open (unit = 12, file = "./Miscellaneous/GenotypedMarkerNRM.txt", status = "unknown")
    endif
  end if

  print*, " "
  print*, " Making NRM"
  PseudoNRM = 0
  do i = 1, nAnisG
    if (mod(i, 400) == 0) print*, "   NRM done for genotyped individual --- ", i
    if (FullFileOutput == 1) then
      shellwarning = 0
      do j = i, nAnisG
	shell = 0
	value = xnumrelmat(RecodeGenotypeId(i), RecodeGenotypeId(j))
	NRM(i, j) = value
	NRM(j, i) = value
      enddo
    endif
    shellwarning = 0
    do j = i, nAnisG
      shell = 0
      valueS = xnumrelmat(seqsire(RecodeGenotypeId(i)), RecodeGenotypeId(j))
      shell = 0
      valueD = xnumrelmat(seqdam(RecodeGenotypeId(i)), RecodeGenotypeId(j))
      if ((valueS > NrmThresh).and.(valueD <= NrmThresh)) PseudoNRM(i, j) = 1
      if ((valueS <= NrmThresh).and.(valueD > NrmThresh)) PseudoNRM(i, j) = 2
      PseudoNRM(j, i) = PseudoNRM(i, j)
    enddo
    if (FullFileOutput == 1) then
      if (nAnisG < 20000) then
	write (8, "(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)") GenotypeId(i), NRM(i,:)
      else
	write (8, *) GenotypeId(i), NRM(i,:)
      end if
      write (9, *) GenotypeId(i), PseudoNRM(i,:)
    endif
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

  deallocate(xnumrelmatHold)
  deallocate(seqid)
  deallocate(seqsire)
  deallocate(seqdam)

  if (FullFileOutput == 1) then
    call MarkerNRMMaker
    print*, "   Finished making marker derived NRM"
    write (10, *) "Average of NRM elements amongst Genotyped Individuals"
    SumDiag = 0
    do i = 1, nAnisG
      SumDiag = SumDiag + NRM(i, i)
    enddo
    write (10, '(f7.4)') (sum(NRM(:,:)) - SumDiag)/(nAnisG * nAnisG)
    write (10, *) " "
    write (10, *) "Average of NRM elements amongst each genotyped individual and all other genotyped Individuals "
    write (10, *) "Diagonal element of NRM for each genotyped individual"
    write (10, *) "Number of genotyped indiviudals related to sire of each genotyped indiviudal above the NRMThresh"
    write (10, *) "Number of genotyped indiviudals related to dam of each genotyped indiviudal above the NRMThresh"
    write (10, *) "Missing genotypes for each indiviudal"
    write (10, *) "Sire genotyped"
    write (10, *) "Dam genotyped"
    do i = 1, nAnisG
      SumNrm = 0
      SumPseudoNrmS = 0
      SumPseudoNrmD = 0
      do j = 1, nAnisG
	if (PseudoNRM(i, j) == 1) SumPseudoNrmS = SumPseudoNrmS + 1
	if (PseudoNRM(i, j) == 2) SumPseudoNrmD = SumPseudoNrmD + 1
	SumNrm = SumNrm + NRM(i, j)
      end do
      CountMissingGenotype = 0
      do j = 1, nSnp
	if (Genos(i, j) == MissingGenotypeCode) CountMissingGenotype = CountMissingGenotype + 1
      end do
      SireGen = 0
      DamGen = 0
      if (SireGenotyped(i) /= 0) SireGen = 1
      if (DamGenotyped(i) /= 0) DamGen = 1
      write (10, '(a20,3f10.4,3i10,2i4)') GenotypeId(i), (SumNrm - NRM(i, i))/(nAnisG - 1), NRM(i, i), MarkerNRM(i, i), &
      SumPseudoNrmS, SumPseudoNrmD, CountMissingGenotype, SireGen, DamGen
    end do
    deallocate(MarkerNRM)
  endif

end subroutine ParseData

!########################################################################################################################################################################################################

subroutine MarkerNRMMaker
  use Global
  use GlobalPedigree
  implicit none

  integer :: i, j, k, nMissing
  double precision :: Sumpq
  double precision, allocatable, dimension (:,:) :: GenosR, tpose

  allocate(GenosR(nAnisG, nSnp))
  !allocate(tpose(nSnp,nAnisG))
  allocate(AlleleFreq(nSnp))
  allocate(MarkerNRM(nAnisG, nAnisG))

  Sumpq = 0.000001
  allelefreq = 0.000001
  do i = 1, nSnp
    nMissing = 0
    do j = 1, nAnisG
      if (Genos(j, i) /= MissingGenotypeCode) then
	AlleleFreq(i) = AlleleFreq(i) + Genos(j, i)
      else
	nMissing = nMissing + 1
      end if
    end do
    AlleleFreq(i) = AlleleFreq(i)/(2 * (nAnisG - nMissing))
    Sumpq = Sumpq + (AlleleFreq(i)*(1 - AlleleFreq(i)))
    write (11, '(i10,f7.4)') i, AlleleFreq(i)
  end do
  print*, " "
  print*, " Making marker derived NRM"

  do i = 1, nAnisG
    do j = 1, nSnp
      !if (Genos(i,j)/=MissingGenotypeCode) then
      !        GenosR(i,j)=float((Genos(i,j)-1))-(2.0*(AlleleFreq(j)-0.5))
      !else
      !        GenosR(i,j)=(AlleleFreq(j)-1.0)-(2.0*(AlleleFreq(j)-0.5))
      !end if
    end do
  end do

  !tpose=transpose(GenosR)
  !MarkerNRM=matmul(GenosR,tpose)
  deallocate(AlleleFreq)

  do i = 1, nAnisG
    do j = 1, nAnisG
      !MarkerNRM(i,j)=MarkerNRM(i,j)/(2*Sumpq)
    end do
  end do
  MarkerNRM = 0.0001

  do i = 1, nAnisG
    write (12, '(a20,20000f6.2,20000f6.2,20000f6.2,20000f6.2,20000f6.2)') GenotypeId(i), MarkerNRM(i,:)
  end do

end subroutine MarkerNRMMaker

!########################################################################################################################################################################################################

recursive function xnumrelmat(i, j) result (xA)
  use GlobalPedigree
  use Global
  implicit none

  real, external :: xnumrelmat_mem
  integer :: i, j
  real :: xA

  shell = shell + 1
  if (i .eq. 0 .or. j .eq. 0 .or. shell > shellmax) then
    xA = 0.0
    if (shell > shellmax) shellWarning = shellWarning + 1
    return
  endif

  If (i .le. NRMmem .And. j .le. NRMmem) Then
    xA = xnumrelmat_mem(i, j)
  else
    IF (i .eq. j)then
    xA = 1 + .5 * xnumrelmat(RecSire(i), RecDam(i))
  elseIF (i .lt. j)then
    xA = .5 * (xnumrelmat(i, RecSire(j)) + xnumrelmat(i, RecDam(j)))
  elseIF (j .lt. i)then
    xA = .5 * (xnumrelmat(j, RecSire(i)) + xnumrelmat(j, RecDam(i)))
    endif
  endif

end function xnumrelmat

!########################################################################################################################################################################################################

RECURSIVE function xnumrelmat_mem(i, j) RESULT (xA)
  use GlobalPedigree
  use Global


  INTEGER :: i, j, k
  REAL :: xA

  ! Addressing vector xnumrelmatHold:  k=(i-1)*(float(NRMmem)-float(i)/2)+j
  shell = shell + 1
  !IF (shell>4000 .AND. MOD(shell,1000).eq.0 )   PRINT*, shell
  if (i .eq. 0 .or. j .eq. 0 .or. shell > shellmax) then
    xA = 0.0
    if (shell > shellmax) shellWarning = shellWarning + 1
    return
  endif

  if (i .le. j) then
    k = (i - 1) * (float(NRMmem) - float(i) / 2) + j
  else
    k = (j - 1) * (float(NRMmem) - float(j) / 2) + i
  endif

  if (xnumrelmatHold(k) > -1) then
    xA = xnumrelmatHold(k)
    return
  endif


  IF (i .eq. j)then
    xA = 1 + .5 * xnumrelmat_mem(RecSire(i), RecDam(i))
  elseIF (i .lt. j)then
    xA = .5 * (xnumrelmat_mem(i, RecSire(j)) + xnumrelmat_mem(i, RecDam(j)))
  elseIF (j .lt. i)then
    xA = .5 * (xnumrelmat_mem(j, RecSire(i)) + xnumrelmat_mem(j, RecDam(i)))
  endif

  xnumrelmatHold(k) = xA


end function xnumrelmat_mem

!########################################################################################################################################################################

subroutine PedigreeViewerRecode(nobs, nAnisPedigree)
  use GlobalPedigree
  implicit none

  integer, allocatable :: passedorder(:)
  integer :: nobs, nAnisPedigree
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
  allocate(id(0:nobs), sire(nobs), dam(nobs), seqid(nobs), seqsire(nobs), seqdam(nobs))

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

    deallocate (seqid)
    allocate(seqid(nobs + iextra))

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
  seqid(i) = i
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

