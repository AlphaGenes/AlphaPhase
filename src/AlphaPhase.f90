program Rlrplhi
  use HaplotypeLibraryDefinition 
  use SurrogateDefinition
  use CoreSubsetDefinition
  use PedigreeDefinition
  use CoreDefinition
  use MemberManagerDefinition
  use ParametersDefinition
  use TestResultDefinition
  use CoreUtils
  use InputOutput
  
  use LongRangePhasing
  use HaplotypeLibraryPhasing
  
  use InputOutput
  
  use HaplotypeModule
  
  implicit none

  integer :: h, i, j, nGlobalHapsOld, threshold
  
  type(HaplotypeLibrary) :: library
  type(Surrogate) :: surrogates
  type(Core) :: c
  type(CoreSubset) :: cs
  type(Pedigree) :: p
  type(Parameters) :: params
  type(TestResults) :: results
  
  integer, allocatable, dimension (:,:,:) :: AllHapAnis
  integer(kind=1), allocatable, dimension(:,:,:) :: AllPhase, Phase, AllTruePhase, TruePhase
  type(Genotype), pointer, dimension(:) :: Genos
  integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp
  integer, dimension (:,:), pointer :: CoreIndex, TailIndex
  integer :: nCores
  integer :: nGlobalHapsIter
  integer :: nAnisG
  integer :: subsetCount
  type(HaplotypeLibrary) :: globalLibrary
  
  type(Core), allocatable, dimension(:) :: AllCores
  
  type(Haplotype), pointer :: hap
  
  !Linux max path length is 4096 which is more than windows or mac (all according to google)
  character(len=4096) :: specfile
  character(len=4096) :: cmd
  
 integer :: startCore, endCore
  logical :: combine, singleRun, printOldProgress, writeSwappable, outputSurrogates
    
  type(MemberManager) :: manager
  
  real :: total, elapsed(2), etime
  
  if (Command_Argument_Count() > 0) then
    call get_command_argument(1,cmd)
    if (cmd(1:2) .eq. "-v") then
      call PrintVersion
      call exit(0)
    end if
  end if

  call Titles
  
  if (Command_Argument_Count() > 0) then
    call Get_Command_Argument(1,specfile)
  else
    specfile="AlphaPhaseSpec.txt"
  end if
  params = ReadInParameterFile(specfile)

  call MakeDirectories(params)
  p = ParsePedigreeData(params)
  nAnisG = p%getNAnis()

  if (params%library .eq. "None") then
    CoreIndex => CalculateCores(params%nSnp, params%Jump, params%offset)    
  else
    CoreIndex => getCoresFromHapLib(params%library)
  end if
  TailIndex => calculateTails(CoreIndex, params%nSnp, params%Jump, params%CoreAndTailLength)
  nCores = size(CoreIndex,1)  
    
  allocate(AllHapAnis(nAnisG, 2, nCores))
  allocate(AllPhase(nAnisG, params%nSnp, 2))
  allocate(AllCores(nCores))
  if (params%GenotypeFileFormat /= 2) then
    Genos => ParseGenotypeData(nAnisG,params)
  else
    AllPhase = ParsePhaseData(params%GenotypeFile,1,params%nSnp,nAnisG,params%nSnp)
  end if

  if (params%Simulation) then
    allocate(AllTruePhase(nAnisG, params%nSnp, 2))
    AllTruePhase = ParsePhaseData(params%TruePhaseFile, 1, params%nSnp, nAnisG, params%nSnp)
  end if  
    
  threshold = int(params%GenotypeMissingErrorPercentage*params%CoreAndTailLength)
  
  singleRun = .true.
  if (params%startCoreChar .eq. "Combine") then
    startCore = nCores + 1
    combine = .true.
    singleRun = .false.
  else
    read(params%startCoreChar, '(i10)') startCore
    combine = .false.
    if (startCore /= 1) then
      singleRun = .false.
    end if
  end if
  
  if (params%endCoreChar .eq. "Combine") then
    endCore = nCores
    combine = .true.
  else
    read(params%endCoreChar, '(i10)') endCore
    combine = .false.
    singleRun = .false.
  end if
  
  printOldProgress = (params%ItterateType .eq. "Off")
  outputSurrogates = (params%ItterateType .eq. "Off") .and. (params%numIter == 1)
  writeSwappable = (params%GenotypeFileFormat /= 2)
  
 
  do h = startCore, endCore
    StartCoreSnp = CoreIndex(h, 1)
    EndCoreSnp = CoreIndex(h, 2)
    StartSurrSnp = TailIndex(h, 1)
    EndSurrSnp = TailIndex(h, 2)
    
    print*, " "
    print*, " "
    print*, " Starting Core", h, "/", nCores
    
    total = etime(elapsed)
    print *, total

    if (params%GenotypeFileFormat /= 2) then
      c = Core(Genos, startSurrSnp, startCoreSnp, endCoreSnp, endSurrSnp)
      if (params%library .eq. "None") then
	library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
      else
	library = getHaplotypeLibrary(params%library, h)
      end if

      do i = 1, params%NumIter	
	manager = MemberManager(c, params%itterateType, params%itterateNumber)

	subsetCount = 0
	do while (manager%hasNext())
	  cs = CoreSubSet(c, p, manager%getNext())

	  surrogates = Surrogate(cs, threshold, printOldProgress)
	  if (outputSurrogates) then
	    call writeSurrogates(surrogates, h, p, params)
	  end if
	  total = etime(elapsed)
	  print *, total
	  call Erdos(surrogates, cs, params%numsurrdisagree, params%useSurrsN, printOldProgress)
	  total = etime(elapsed)
	  print *, total
	  call CheckCompatHapGeno(cs, params%percgenohaplodisagree, printOldProgress)
	  total = etime(elapsed)
	  print *, total

	  subsetCount = subsetCount + 1
	  if (.not. printOldProgress) then
	    print '(8x, i5, a20, f6.2, a16, f6.2, a16)', subsetCount, " Subsets completed, ", c%getYield(1), "% Paternal yield, ", &
	      c%getYield(2), "% Maternal Yield"
	  end if
	end do

	globalLibrary = HaplotypeLibrary(c%getNCoreSnp(),500,500)
	nGlobalHapsIter = 1
	call UpdateHapLib(c, library)
	nGlobalHapsOld = library%getSize()
	if (params%ItterateType .eq. "Off") then
	  print*, " "
	  print*, "  ", "Haplotype library imputation step"
	end if
	do j = 1, 20
	  call ImputeFromLib(library, c, nGlobalHapsIter, params%PercGenoHaploDisagree, params%minHapFreq)
	  call UpdateHapLib(c,library)
	  if (nGlobalHapsOld == library%getSize()) exit
	  nGlobalHapsOld = library%getSize()
	end do

	if (.not. printOldProgress) then
	  print '(4x, a9, 20x, i6, a19, f6.2, a25)', "After HLI", library%getSize(), " Haplotypes found, ", &
	    c%getPercentFullyPhased(), "% Haplotypes fully phased"
	  print '(33x, f6.2, a16, f6.2, a16)', c%getYield(1), "% Paternal yield, ", c%getYield(2), "% Maternal Yield"
	end if

	total = etime(elapsed)
	print *, total
      end do
    else
      allocate(Phase(nAnisG,EndCoreSnp-StartCoreSnp+1,2))
      Phase = AllPhase(:,StartCoreSnp:EndCoreSnp,:)
      c = Core(Phase)
      do i = 1, nAnisG
	call c%setHaplotype(i,1,Haplotype(Phase(i,:,1)))
	call c%setHaplotype(i,2,Haplotype(Phase(i,:,2)))
      end do
      library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
      call UpdateHapLib(c,library)
      deallocate(Phase)
    end if
   
    call WriteHapLib(library, h, c, params)
    
    if (params%outputHapCommonality) then
      call HapCommonality(library, h, params)
    end if
    
    if (.not. singleRun) then
      call WriteOutCore(c, h, CoreIndex(h,1), p, writeSwappable, params)
    else
      do i = 1, size(AllPhase,1)
	hap => c%phase(i,1)
	AllPhase(i,startCoreSnp:endCoreSnp,1) = hap%toIntegerArray()
	hap => c%phase(i,2)
	AllPhase(i,startCoreSnp:endCoreSnp,2) = hap%toIntegerArray()
      end do
      AllHapAnis(:,1,h) = c%hapAnis(:,1)
      AllHapAnis(:,2,h) = c%hapAnis(:,2)

      AllCores(h) = c
    end if
    
    if (params%Simulation) then
      allocate(TruePhase(nAnisG,EndCoreSnp-StartCoreSnp+1,2))
      TruePhase = AllTruePhase(:,StartCoreSnp:EndCoreSnp,:)
      
      call Flipper(c,TruePhase)
      results = TestResults(c,TruePhase)
      call WriteTestResults(results,c,surrogates,p,h,outputSurrogates,params)
      call WriteMistakes(c,TruePhase,p,h,params)
      
      deallocate(TruePhase)
    end if
  end do
  
  if ((.not. SingleRun) .and. combine) then
    call CombineResults(nAnisG,CoreIndex,writeSwappable,params)
  else
    call WriteOutResults(AllCores,CoreIndex,p,writeSwappable, params)
  end if
  deallocate(AllHapAnis)
  deallocate(AllPhase)
  if (params%Simulation) then
    deallocate(AllTruePhase)
  end if
  if (params%Simulation) then
    call CombineTestResults(nCores,params)
  end if
  
  deallocate(CoreIndex,TailIndex)
  
  call PrintTimerTitles(params)

end program Rlrplhi