
program AlphaPhase

  implicit none

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
    integer :: id
    integer :: nAnisG
    type(AlphaPhaseResults) :: results
    logical :: singleSurrogates, singleRun, combine, notPrephased
    integer :: i

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

    call MakeDirectories(params)
    p = ParsePedigreeAndGenotypeData(params)
    nAnisG = p%nHd

    notPrephased = (params%GenotypeFileFormat /= 2)

    if (notPrephased) then
      if (params%Simulation) then
        TruePhase => ParsePhaseData(params%TruePhaseFile,nAnisG,params%nSnp)
        if (params%Library .ne. "None") then
    existingLibraries => getHaplotypeLibraries(params%library)
    results = phaseAndCreateLibraries(p, params%params, existingLibraries, TruePhase, quiet = .false.)
        else
    results = phaseAndCreateLibraries(p, params%params, truePhase = TruePhase, quiet = .false.)
        end if
      else
        if (params%Library .ne. "None") then
    existingLibraries => getHaplotypeLibraries(params%library)
    results = phaseAndCreateLibraries(p, params%params, existingLibraries, quiet = .false.)
        else
    results = phaseAndCreateLibraries( p, params%params, quiet = .false.)
        end if
      end if
    else
      Phase => ParsePhaseData(params%GenotypeFile,nAnisG,params%nSnp)
      results = createLibraries(Phase, params%params)
    end if

    singleSurrogates = (params%params%ItterateType .eq. "Off") .and. (params%params%numIter == 1)
    singleRun = (params%params%StartCoreChar .eq. "1") .and. (params%params%EndCoreChar .eq. "Combine")
    combine = (params%params%EndCoreChar .eq. "Combine")

    do i = 1, results%nCores
      id = results%ids(i)
      call WriteHapLib(results%libraries(i), id, params)
      if (params%outputParams%outputHapCommonality) then
        call HapCommonality(results%libraries(i), id, params)
      end if
      if (.not. singleRun) then
        call WriteOutCore(results%cores(i), id, results%startIndexes(i), p, notPrephased, params)
      end if
      if (singleSurrogates) then
        call writeSurrogates(results%surrogates(i), id, p, params)
      end if
      if (params%Simulation) then
        call WriteTestResults(results%testResults(i),results%cores(i),p,id,params)
      end if
    end do

    if ((.not. SingleRun) .and. combine) then
      call CombineResults(nAnisG,notPrephased,params)
    else
      call WriteOutResults(results%cores,results%startIndexes,results%endIndexes,p,notPrephased,params)
    end if
    if (params%Simulation) then
      call CombineTestResults(results%nCores,params)
    end if

    call PrintTimerTitles(params)
  end subroutine main

end program AlphaPhase