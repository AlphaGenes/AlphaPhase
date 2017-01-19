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
  use AlphaPhaseResultsDefinition
  
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
  type(TestResults) :: test
  
  type(Genotype), pointer, dimension(:) :: Genos
  type(Haplotype), pointer, dimension(:,:) :: Phase, TruePhase, CoreTruePhase
  integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp
  integer, dimension (:,:), pointer :: CoreIndex, TailIndex
  integer :: nCores
  integer :: nGlobalHapsIter
  integer :: nAnisG
  integer :: subsetCount
  type(HaplotypeLibrary) :: globalLibrary
  type(Haplotype), pointer :: hap
  type(Haplotype), target :: stupid, tempHap
  type(Core), allocatable, dimension(:) :: AllCores
  type(AlphaPhaseResults) :: results
  integer :: id
  
  !Linux max path length is 4096 which is more than windows or mac (all according to google)
  character(len=4096) :: specfile
  character(len=4096) :: cmd
  
  integer :: startCore, endCore
  logical :: combine, singleRun, printOldProgress, writeSwappable, singleSurrogates
    
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
    
  allocate(AllCores(nCores))
  if (params%GenotypeFileFormat /= 2) then
    Genos => ParseGenotypeData(nAnisG,params)
  else
    Phase => ParsePhaseData(params%GenotypeFile,nAnisG,params%nSnp)
  end if

  if (params%Simulation) then
    TruePhase => ParsePhaseData(params%TruePhaseFile,nAnisG,params%nSnp)
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
  singleSurrogates = (params%ItterateType .eq. "Off") .and. (params%numIter == 1)
  writeSwappable = (params%GenotypeFileFormat /= 2)
  
 
  results = AlphaPhaseResults(endCore-startCore+1,singleSurrogates,params%Simulation)
  
  do h = startCore, endCore
    StartCoreSnp = CoreIndex(h, 1)
    EndCoreSnp = CoreIndex(h, 2)
    StartSurrSnp = TailIndex(h, 1)
    EndSurrSnp = TailIndex(h, 2)
    
    print*, " "
    print*, " "
    print*, " Starting Core", h, "/", nCores

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
	  if (singleSurrogates) then
	    results%surrogates(h-startCore+1) = surrogates
	  end if
	  call Erdos(surrogates, cs, params%numsurrdisagree, params%useSurrsN, printOldProgress)
	  call CheckCompatHapGeno(cs, params%percgenohaplodisagree, printOldProgress)

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
      end do
    else
      c = Core(Phase,StartCoreSnp,EndCoreSnp)
      library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
      call UpdateHapLib(c,library)
    end if    
    
    if (params%ItterateType .eq. "Off") then
      print*, "   ", "Final iteration found ", library%getSize(), "haplotypes"
    
      print*, ""
      write (*, '(a4,a30,f5.2,a1)') "   ", "Final yield for this core was ", c%getTotalYield(), "%"
    end if
    
    if (params%Simulation) then
      allocate(CoreTruePhase(nAnisG,2))
      do i = 1, nAnisG
	hap => TruePhase(i,1)
	CoreTruePhase(i,1) = hap%subset(startCoreSnp,endCoreSnp)
	hap => TruePhase(i,2)
	CoreTruePhase(i,2) = hap%subset(startCoreSnp,endCoreSnp)
      end do
      call c%flipHaplotypes(CoreTruePhase)
      results%testResults(h-startCore+1) = TestResults(c,CoreTruePhase)
      deallocate(CoreTruePhase)
    end if
   
    results%libraries(h-startCore+1) = library
    results%cores(h-startCore+1) = c
    results%ids(h-startCore+1) = h
    results%startIndexes(h-startCore+1) = CoreIndex(h,1)
    results%endIndexes(h-startCore+1) = CoreIndex(h,2)
  end do
  
  if ((.not. SingleRun) .and. combine) then
    call CombineResults(nAnisG,CoreIndex,writeSwappable,params)
  else
    call WriteOutResults(results%cores,CoreIndex,p,writeSwappable, params)
  end if
  if (params%Simulation) then
    call CombineTestResults(nCores,params)
  end if
  
  do i = 1, endCore-startCore+1
    id = results%ids(i)
    call WriteHapLib(results%libraries(i), id, params)
    if (params%outputHapCommonality) then
      call HapCommonality(results%libraries(i), id, params)
    end if
    if (.not. singleRun) then
      call WriteOutCore(results%cores(i), id, results%startIndexes(i), p, writeSwappable, params)
    end if
    if (singleSurrogates) then
      call writeSurrogates(results%surrogates(i), id, p, params)
    end if
    if (params%Simulation) then
      call WriteTestResults(results%testResults(i),results%cores(i),p,id,params)
      allocate(CoreTruePhase(nAnisG,2))
      do j = 1, nAnisG
	hap => TruePhase(j,1)
	CoreTruePhase(j,1) = hap%subset(startCoreSnp,endCoreSnp)
	hap => TruePhase(j,2)
	CoreTruePhase(j,2) = hap%subset(startCoreSnp,endCoreSnp)
      end do
      call WriteMistakes(results%cores(i),CoreTruePhase,p,id,params)  
    end if
  end do
  
  deallocate(CoreIndex,TailIndex)
  
  call PrintTimerTitles(params)

end program Rlrplhi