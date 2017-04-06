
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

    !Linux max path length is 4096 which is more than windows or mac (all according to google)
    character(len=4096) :: specfile
    character(len=4096) :: cmd

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

    p = ParsePedigreeAndGenotypeData(params)
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