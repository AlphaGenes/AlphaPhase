
program AlphaPhase

  implicit none

  ! Main is a seperate subroutine and called here due to compiler issues (Roberto / 16.0.3)
  call main

  contains

  subroutine main
    use HaplotypeLibraryDefinition
    use ProgramParametersDefinition
    use InputOutput
    use AlphaPhaseResultsDefinition
    use AlphaPhaseFunctions
    use InputOutput
    use PedigreeModule
    use HaplotypeModule

    implicit none
    type(ProgramParameters) :: params
    type(PedigreeHolder) :: p
    type(Haplotype), pointer, dimension(:,:) :: Phase, TruePhase
    type(HaplotypeLibrary), dimension(:), pointer :: existingLibraries
    integer :: nAnisG
    type(AlphaPhaseResults) :: results
    logical :: notPrephased
    integer, dimension(:,:), pointer :: CoreIndex
    integer :: specloc

    !Linux max path length is 4096 which is more than windows or mac (all according to google)
    character(len=4096) :: specfile
    character(len=4096) :: cmd
    character(len=4096) :: limit

    if (Command_Argument_Count() > 0) then
      call get_command_argument(1,cmd)
      if (cmd(1:2) .eq. "-v") then
        call PrintVersion
        call exit(0)
      end if
      if ((cmd(1:2) .eq. "-s") .or. (cmd(1:2) .eq. "-c")) then
	call Titles
	if (Command_Argument_Count() > 1) then
	  call Get_Command_Argument(2,specfile)
	else
	  specfile="AlphaPhaseSpec.txt"
	end if
	params = ReadInParameterFile(specfile)
	if (cmd(1:2) .eq. "-s") then
	  CoreIndex => calculateCores(params%nSnp,params%params%Jump,params%params%offset)
	  p = ParsePedigreeAndGenotypeData(params)
	  
	  print *
	  call printCoreInfo(CoreIndex)
	  
	  call makeDirectories(params%outputParams)
	  
	  call writeCoreIndex(params%outputParams, size(CoreIndex,1), p%nGenotyped, params%nSnp, CoreIndex(:,1), CoreIndex(:,2))
	  print *
	  print *, "Setup complete"
	  call exit(0)
	end if
	if (cmd(1:2) .eq. "-c") then
	  params%outputParams%outputPerCore = .false.
	  params%outputParams%outputCombined = .true.
	  
	  ! The following are already done on a per core basis so will already have been created and aren't created right anyway.
	  ! Bit hacky and should probably change in future but works for now.
	  params%outputParams%outputHaplotypeLibraryText = .false.
	  params%outputParams%outputHaplotypeLibraryBinary = .false.
	  params%outputParams%outputSurrogatesSummary = .false.
	  params%outputParams%outputSurrogates = .false.
	  params%outputParams%outputHapCommonality = .false.
	  
	  p = ParsePedigreeAndGenotypeData(params)
	  call readInPerCoreResults(results, params%outputParams, p)
	  call writeAlphaPhaseResults(results, p, params%outputParams)
	  print *, "Results combined"
	  call exit(0)
	end if
      end if
    end if

    call Titles

    if (Command_Argument_Count() > 0) then
      if (cmd(1:2) .eq. "-r") then 
	specloc = 3
      else
	specloc = 1
      end if
      if (specloc <= Command_Argument_Count()) then
	call Get_Command_Argument(specloc,specfile)
      else
	specfile="AlphaPhaseSpec.txt"
      end if
    else
      specfile="AlphaPhaseSpec.txt"
    end if
    params = ReadInParameterFile(specfile)

    p = ParsePedigreeAndGenotypeData(params)
    
    if (Command_Argument_Count() > 0) then
      if (cmd(1:2) .eq. "-r") then
	call get_command_argument(2,limit)
	if (index(limit,"-") > 0) then
	  params%params%startCoreChar = limit(1:index(limit,"-")-1)
	  params%params%endCoreChar = trim(limit(index(limit,"-")+1:4096))
	else
	  params%params%startCoreChar = limit
	  params%params%endCoreChar = limit
	end if
	params%outputParams%outputCombined = .false.
	params%outputParams%outputPerCore = .true.
	params%outputParams%outputTimer = .false.
      end if
    end if
    
    nAnisG = p%nHd

    notPrephased = (params%GenotypeFileFormat /= 2)
    
    if (params%CoreFile .ne. "None") then
      CoreIndex => readInCores(params%CoreFile)
    end if

    if (notPrephased) then
      if (params%Simulation) then
        TruePhase => ParsePhaseData(params%TruePhaseFile,nAnisG,params%nSnp)
        if (params%Library .ne. "None") then
	  existingLibraries => getHaplotypeLibraries(params%library)
	  results = phaseAndCreateLibraries(p, params%params, existingLibraries, TruePhase, quiet = .false.)
        else
	  if (params%CoreFile .ne. "None") then
	    results = phaseAndCreateLibraries(p, params%params, truePhase = TruePhase, userCoreIndex = coreIndex, quiet = .false.)
	  else
	    results = phaseAndCreateLibraries(p, params%params, truePhase = TruePhase, quiet = .false.)
	  end if
        end if
      else
        if (params%Library .ne. "None") then
	  existingLibraries => getHaplotypeLibraries(params%library)
	  results = phaseAndCreateLibraries(p, params%params, existingLibraries, quiet = .false.)
        else
	  if (params%CoreFile .ne. "None") then
	    results = phaseAndCreateLibraries(p, params%params, userCoreIndex = coreIndex, quiet = .false.)
	  else
	    results = phaseAndCreateLibraries(p, params%params, quiet = .false.)
	  end if
        end if
      end if
    else
      Phase => ParsePhaseData(params%GenotypeFile,nAnisG,params%nSnp)
      if (params%Library .ne. "None") then
	existingLibraries => getHaplotypeLibraries(params%library)
	results = createLibraries(Phase, params%params, existingLibraries)
      else
	if (params%CoreFile .ne. "None") then
	  results = createLibraries(Phase, params%params, userCoreIndex = coreIndex)
	else
	  results = createLibraries(Phase, params%params)
	end if
      end if
    end if

    call writeAlphaPhaseResults(results, p, params%outputParams)    
    
    call PrintTimerTitles(params%outputParams)
  end subroutine main

end program AlphaPhase