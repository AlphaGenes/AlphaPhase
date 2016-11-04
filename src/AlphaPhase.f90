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
  
  use NRMcode
  
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
  integer(kind=1), allocatable, dimension(:,:) :: AllGenos, Genos
  integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp
  integer, dimension (:,:), pointer :: CoreIndex, TailIndex
  integer :: nCores
  integer :: nGlobalHapsIter
  integer :: nAnisG
  integer :: subsetCount
  
  type(Core), allocatable, dimension(:) :: AllCores
  
  !Linux max path length is 4096 which is more than windows or mac (all according to google)
  character(len=4096) :: specfile
  character(len=4096) :: cmd
  
  integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM
  
  integer :: startCore, endCore
  logical :: combine, printOldProgress, writeSwappable, outputSurrogates
  
  type(MemberManager) :: manager
  
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
    CoreIndex => CalculateCores(params%nSnp, params%Jump, params%offset, params%consistent)    
  else
    CoreIndex => getCoresFromHapLib(params%library)
  end if
  TailIndex => calculateTails(CoreIndex, params%nSnp, params%Jump, params%CoreAndTailLength, params%offset, params%consistent)
  nCores = size(CoreIndex,1)  
    
  if (.not. params%readCoreAtTime) then
    allocate(AllHapAnis(nAnisG, 2, nCores))
    allocate(AllPhase(nAnisG, params%nSnp, 2))
    allocate(AllCores(nCores))
    if (params%GenotypeFileFormat /= 2) then
      allocate(AllGenos(nAnisG,params%nSnp))
      AllGenos = ParseGenotypeData(1,params%nSnp,nAnisG,params)
    else
      AllPhase = ParsePhaseData(params%GenotypeFile,1,params%nSnp,nAnisG,params%nSnp)
    end if
    
    if (params%Simulation) then
      allocate(AllTruePhase(nAnisG, params%nSnp, 2))
      AllTruePhase = ParsePhaseData(params%TruePhaseFile, 1, params%nSnp, nAnisG, params%nSnp)
    end if
  end if
  
  if ((params%consistent) .and. (params%GenotypeFileFormat /= 2)) then
    allocate(PseudoNRM(nAnisG,nAnisG))
    PseudoNRM = createNRM(p, params)
  end if
    
  threshold = int(params%GenotypeMissingErrorPercentage*params%CoreAndTailLength)
  
  if (params%startCoreChar .eq. "Combine") then
    startCore = nCores
    combine = .true.
  else
    read(params%startCoreChar, '(i10)') startCore
    combine = .false.
  end if
  
  if (params%endCoreChar .eq. "Combine") then
    endCore = nCores
    combine = .true.
  else
    read(params%endCoreChar, '(i10)') endCore
    combine = .false.
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
    
    if (params%GenotypeFileFormat /= 2) then
      allocate(Genos(nAnisG, max(EndSurrSnp,EndCoreSnp)-startSurrSnp+1))
      if (params%readCoreAtTime) then
	Genos = ParseGenotypeData(StartSurrSnp,max(EndSurrSnp,EndCoreSnp),nAnisG,params)
      else
	Genos = AllGenos(:,StartSurrSnp:max(EndSurrSnp,EndCoreSnp))
      end if
      
      ! Fudge below
      c = Core(Genos, startCoreSnp-startSurrSnp+1, endCoreSnp-startSurrSnp+1, endSurrSnp-startSurrSnp+1)
      
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

	  surrogates = Surrogate(cs, threshold, params%consistent, pseudoNRM, printOldProgress)
	  if (outputSurrogates) then
	    call writeSurrogates(surrogates,threshold, h, p, params)
	  end if
	  call Erdos(surrogates, cs, threshold, params%numsurrdisagree, params%useSurrsN, params%consistent, printOldProgress)
	  call CheckCompatHapGeno(cs, params%percgenohaplodisagree, printOldProgress)     

	  subsetCount = subsetCount + 1
	  if (.not. printOldProgress) then
	    print '(8x, i5, a20, f6.2, a16, f6.2, a16)', subsetCount, " Subsets completed, ", c%getYield(1), "% Paternal yield, ", &
	      c%getYield(2), "% Maternal Yield"
	  end if
	end do

	nGlobalHapsIter = 1
	if (params%consistent) then
	  library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
	end if
	call MakeHapLib(c, library, params%consistent)
	nGlobalHapsOld = library%getSize()
	if (params%ItterateType .eq. "Off") then
	  print*, " "
	  print*, "  ", "Haplotype library imputation step"
	end if
	do j = 1, 20
	  call ImputeFromLib(library, c, nGlobalHapsIter, params%PercGenoHaploDisagree, params%minHapFreq, params%consistent)
	  if (params%consistent) then
	    library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
	  end if
	  call MakeHapLib(c,library,params%consistent)
	  if (nGlobalHapsOld == library%getSize()) exit
	  nGlobalHapsOld = library%getSize()
	end do

	if (.not. printOldProgress) then
	  print '(4x, a9, 20x, i6, a19, f6.2, a25)', "After HLI", library%getSize(), " Haplotypes found, ", &
	    c%getPercentFullyPhased(), "% Haplotypes fully phased"
	  print '(33x, f6.2, a16, f6.2, a16)', c%getYield(1), "% Paternal yield, ", c%getYield(2), "% Maternal Yield"
	end if

      end do

      deallocate(Genos)
    else
      allocate(Phase(nAnisG,EndCoreSnp-StartCoreSnp+1,2))
      if (params%readCoreAtTime) then
	Phase = ParsePhaseData(params%GenotypeFile,StartCoreSnp,EndCoreSnp,nAnisG,params%nSnp)
      else
	Phase = AllPhase(:,StartCoreSnp:EndCoreSnp,:)
      end if
      c = Core(Phase)
      do i = 1, nAnisG
	call c%setHaplotype(i,1,Phase(i,:,1))
	call c%setHaplotype(i,2,Phase(i,:,2))
      end do
      library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
      call MakeHapLib(c,library,params%consistent)
      deallocate(Phase)
    end if
   
    call WriteHapLib(library, h, c, params)
    
    if (params%consistent) then
      call HapCommonality(library, h, params)
    end if
    
    if (params%readCoreAtTime .or. .not. combine) then
      call WriteOutCore(c, h, CoreIndex(h,1), p, writeSwappable, params)
    else
      AllPhase(:,startCoreSnp:endCoreSnp,:) = c%getAllPhase()
      AllHapAnis(:,1,h) = c%hapAnis(:,1)
      AllHapAnis(:,2,h) = c%hapAnis(:,2)

      AllCores(h) = c
    end if
    
    if (params%Simulation) then
      allocate(TruePhase(nAnisG,EndCoreSnp-StartCoreSnp+1,2))
      if (params%readCoreAtTime) then
	TruePhase = ParsePhaseData(params%TruePhaseFile,StartCoreSnp,EndCoreSnp,nAnisG,params%nSnp)
      else
	TruePhase = AllTruePhase(:,StartCoreSnp:EndCoreSnp,:)
      end if
      
      call Flipper(c,TruePhase)
      results = TestResults(c,TruePhase)
      call WriteTestResults(results,c,surrogates,p,TruePhase,h,outputSurrogates,params)
      call WriteMistakes(c,TruePhase,p,h,params)
      
      deallocate(TruePhase)
    end if
  end do
  
  if (params%readCoreAtTime .and. combine) then
    call CombineResults(nAnisG,CoreIndex,p,writeSwappable,params)
  else
    call WriteOutResults(AllCores,CoreIndex,p,writeSwappable, params)
  end if
  if ((params%consistent) .and. (params%GenotypeFileFormat /= 2)) then
    deallocate(PseudoNRM)
  end if
  if (.not. params%readCoreAtTime) then
    deallocate(AllHapAnis)
    deallocate(AllPhase)
    if (params%GenotypeFileFormat /= 2) then
      deallocate(AllGenos)
    end if
    if (params%Simulation) then
      deallocate(AllTruePhase)
    end if
  end if
  if (params%Simulation) then
    call CombineTestResults(nCores,params)
  end if
  
  deallocate(CoreIndex,TailIndex)
  
  call PrintTimerTitles(params)

end program Rlrplhi