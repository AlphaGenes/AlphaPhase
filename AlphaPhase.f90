program Rlrplhi
  use HaplotypeLibraryDefinition
  use SurrogateDefinition
  use CoreSubsetDefinition
  use PedigreeDefinition
  use CoreDefinition
  use MemberManagerDefinition
  use ParametersDefinition
  use TestResultDefinition
  
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
  integer, allocatable, dimension (:,:) :: CoreIndex, TailIndex
  integer :: nCores
  integer :: nGlobalHapsIter
  integer :: nAnisG
  integer :: subsetCount
  
  type(Core), allocatable, dimension(:) :: AllCores
  
  !Linux max path length is 4096 which is more than windows or mac (all according to google)
  character(len=4096) specfile
  character(len=4096) :: cmd
  
  integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM
  
  integer :: startCore, endCore
  logical :: combine, printOldProgress, writeSwappable, outputSurrogates
  
  type(MemberManager) :: manager
  
  interface calculateCores
    subroutine calculateCores(CoreIndex, TailIndex, nSnp, Jump, CoreAndTailLength, offset, consistent)
      implicit none

      integer, dimension(:,:), allocatable :: CoreIndex, TailIndex
      integer, intent(in) :: nSnp, Jump, CoreAndTailLength
      logical, intent(in) :: offset, consistent
    end subroutine calculateCores
  end interface calculateCores
  
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
  call CalculateCores(CoreIndex, TailIndex, params%nSnp, params%Jump, params%CoreAndTailLength, params%offset, params%consistent)
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

      do i = 1, params%NumIter
	manager = MemberManager(c, params%itterateType, params%itterateNumber)

	subsetCount = 0
	do while (manager%hasNext())
	  cs = CoreSubSet(c, p, manager%getNext())

	  surrogates = Surrogate(cs, threshold, params%consistent, pseudoNRM, printOldProgress)
	  if (outputSurrogates) then
	    call writeSurrogates(surrogates,threshold, h, p)
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
	library = MakeHapLib(c, params%consistent)
	nGlobalHapsOld = library%getSize()
	if (params%ItterateType .eq. "Off") then
	  print*, " "
	  print*, "  ", "Haplotype library imputation step"
	end if
	do j = 1, 20
	  call ImputeFromLib(library, c, nGlobalHapsIter, params%PercGenoHaploDisagree, params%minHapFreq, params%consistent)
	  library = MakeHapLib(c,params%consistent)
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
      library = MakeHapLib(c,params%consistent)
      deallocate(Phase)
    end if
   
    call WriteHapLib(library, h, c, params)
    
    if (params%consistent) then
      call HapCommonality(library, h, params)
    end if
    
    if (params%readCoreAtTime .or. .not. combine) then
      call WriteOutCore(c, h, CoreIndex(h,1), p, writeSwappable)
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
      if (params%FullFileOutput) then
	call WriteTestResults(results,c,surrogates,p,TruePhase,h,outputSurrogates)
	call WriteMistakes(c,TruePhase,p,h)
      end if
      
      deallocate(TruePhase)
    end if
  end do
  
  if (params%readCoreAtTime .and. combine) then
    call CombineResults(nAnisG,CoreIndex,p,writeSwappable)
  else
    !call WriteOutResults(AllPhase,AllHapAnis,CoreIndex,p)
    call WriteOutResults(AllCores,CoreIndex,p,writeSwappable)
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
  if ((params%Simulation) .and. combine .and. (params%FullFileOutput)) then
    call CombineTestResults(nCores)
  end if
  
  deallocate(CoreIndex,TailIndex)
  
  call PrintTimerTitles

end program Rlrplhi

!####################################################################################################################################################################

subroutine calculateCores(CoreIndex, TailIndex, nSnp, Jump, CoreAndTailLength, offset, consistent)
  implicit none
  
  integer, dimension(:,:), allocatable :: CoreIndex, TailIndex
  integer, intent(in) :: nSnp, Jump, CoreAndTailLength
  logical, intent(in) :: offset, consistent
  
  integer :: resid
  double precision :: corelength
  integer :: left, ltail, rtail, nCores
  integer :: i
  
  if (consistent) then
    if (.not. Offset) then
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

    if (Offset) then
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
    
    if (.not. Offset) then
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

!####################################################################################################################################################################
subroutine Header
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
end subroutine Header

!####################################################################################################################################################################

subroutine Titles

  call Header
  print*, ""
  print*, ""
  print*, ""
  print*, ""

end subroutine Titles

!###################################################################################################################################################

subroutine PrintTimerTitles
  use Constants
  use InputOutput

  implicit none
  real :: etime ! Declare the type of etime()
  real :: elapsed(2) ! For receiving user and system time
  real :: total, Minutes, Hours, Seconds

  print*, ""
  print*, ""
  call Header
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
  call writeTimer(INT(Hours),INT(Minutes),Seconds)
end subroutine PrintTimerTitles

!###################################################################################################################################################

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

subroutine PrintVersion
  call Header
  print *
  print *, "                              Commit:   "//TOSTRING(COMMIT)
  print *, "                              Compiled: "//__DATE__//", "//__TIME__  
end subroutine PrintVersion
