
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaPhaseFunctions.f90
!
! DESCRIPTION:
!> @brief    Module allowing for functionality of alphaphase to be run independently
!
!
!> @date     January 4, 2017
!
!> @version  1.0.0
!
!
!-------------------------------------------------------------------------------

module AlphaPhaseFunctions
  use HaplotypeLibraryDefinition
  use SurrogateDefinition
  use CoreSubsetDefinition
  use PedigreeModule
  use CoreDefinition
  use MemberManagerDefinition
  use AlphaPhaseParametersDefinition
  use TestResultDefinition
  use CoreUtils
  use InputOutput
  use AlphaPhaseResultsDefinition

  use LongRangePhasing
  use HaplotypeLibraryPhasing

  use InputOutput

  use HaplotypeModule
  implicit none

contains
  function phaseAndCreateLibraries(p, params, existingLibraries, truePhase, quiet) result(results)
    ! Following use statements needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeLibraryDefinition
    use PedigreeModule
    use HaplotypeModule

    type(PedigreeHolder), intent(inout) :: p
    type(AlphaPhaseParameters) :: params
    type(HaplotypeLibrary), pointer, dimension(:), intent(in), optional :: existingLibraries
    type(Haplotype), pointer, dimension(:,:), intent(in), optional :: truePhase
    logical, optional :: quiet

    integer :: h, i, threshold

    type(HaplotypeLibrary) :: library
    type(Surrogate) :: surrogates
    type(Core) :: c
    type(CoreSubset) :: cs
    type(Haplotype), pointer, dimension(:,:) :: CoreTruePhase
    integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp
    integer, dimension (:,:), pointer :: CoreIndex, TailIndex
    integer :: nCores
    integer :: nAnisG, nSnp
    integer :: subsetCount
    type(Haplotype), pointer :: hap
    type(AlphaPhaseResults) :: results

    integer :: startCore, endCore
    logical :: printOldProgress, singleSurrogates, quietInternal

    type(MemberManager) :: manager

    if (.not. present(quiet)) then
      quietInternal = .true.
    else
      quietInternal = quiet
    end if

    if (p%nHd == 0) then
    ! TODO check if this is wanted behaviour
      p%nHd = p%nGenotyped
      p%hdMap = p%genotypeMap
      p%hdDictionary = p%genotypeDictionary
    endif
    nAnisG = p%nHd
    nSnp = p%pedigree(p%hdMap(1))%individualGenotype%getLength()

    if (.not. present(existingLibraries)) then
      CoreIndex => CalculateCores(nSnp, params%Jump, params%offset)
    else
      CoreIndex => getCoresFromLibraries(existingLibraries)
    end if
    TailIndex => calculateTails(CoreIndex, nSnp, params%Jump, params%CoreAndTailLength)
    nCores = size(CoreIndex,1)

    threshold = int(params%GenotypeMissingErrorPercentage*params%CoreAndTailLength)

    if (params%startCoreChar .eq. "0") then
      startCore = 1
    else
      read(params%startCoreChar, '(i10)') startCore
    end if

    if (params%endCoreChar .eq. "0") then
      endCore = nCores
    else
      read(params%endCoreChar, '(i10)') endCore
    end if

    printOldProgress = (params%iterateType .eq. "Off") .and. (.not. quietInternal)
    singleSurrogates = (params%iterateType .eq. "Off") .and. (params%numIter == 1)

    results = AlphaPhaseResults(endCore-startCore+1,singleSurrogates,present(TruePhase))

    do h = startCore, endCore
      StartCoreSnp = CoreIndex(h, 1)
      EndCoreSnp = CoreIndex(h, 2)
      StartSurrSnp = TailIndex(h, 1)
      EndSurrSnp = TailIndex(h, 2)

      if (.not. quietInternal) then
	print*, " "
	print*, "Starting Core", h, "/", nCores
      end if

      c = Core(p, startSurrSnp, startCoreSnp, endCoreSnp, endSurrSnp)
      if (.not. present(existingLibraries)) then
	library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
      else
	library = existingLibraries(h)
      end if

      do i = 1, params%NumIter
	manager = MemberManager(c, params%iterateType, params%iterateNumber)

	subsetCount = 0
	do while (manager%hasNext())
	  cs = CoreSubSet(c, p, manager%getNext())

	  surrogates = Surrogate(cs, threshold, params%minOverlap)
	  if (singleSurrogates) then
	    results%surrogates(h-startCore+1) = surrogates
	  end if
	  call Erdos(surrogates, cs, params%numsurrdisagree, params%useSurrsN)
	  call CheckCompatHapGeno(cs, params%percgenohaplodisagree)

	  subsetCount = subsetCount + 1
	  if ((.not. quietInternal) .and. (params%iterateType .ne. "Off")) then
	    print '(6x, i5, a21, f6.2, a7)', subsetCount, " Subsets completed ", c%getTotalYield(), "% Yield"
	  end if
	end do
	
	call UpdateHapLib(c, library, params%percminpresent, params%minoverlap, params%PercGenoHaploDisagree)
	if (.not. quietInternal) then
	  if (params%iterateType .eq. "Off") then
	    print '(6x, a15, 11x, f6.2, a8, i7, a11)', " LRP completed ", c%getTotalYield(), "% Yield ", &
	      library%numberPercentPhased(params%percMinToKeep), " Haplotypes"
	  end if
	  print*, "   Haplotype Library Imputation step"
	end if
	call imputeFromLib(library, c, params%PercGenoHaploDisagree, params%percMinPresent, params%minoverlap, params%minHapFreq, &
	  params%percMinToKeep, quietInternal)
      end do
      
      library = library%rationalise(params%percMinToKeep,c)

      if (.not. quietInternal) then
	print *, "   Core complete"
      end if

      if (present(TruePhase)) then
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
  end function phaseAndCreateLibraries

  function createLibraries(phase, params, existingLibraries) result (results)
    ! use HaplotypeModule needed here due to compiler issues (Roberto / 16.0.3)
    use HaplotypeLibraryDefinition
    use HaplotypeModule

    type(Haplotype), pointer, dimension(:,:), intent(in) :: phase
    type(AlphaPhaseParameters) :: params
    type(HaplotypeLibrary), pointer, dimension(:), intent(in), optional :: existingLibraries

    type(AlphaPhaseResults) :: results

    integer :: h

    type(HaplotypeLibrary) :: library
    type(Core) :: c
    integer :: StartSurrSnp, EndSurrSnp, StartCoreSnp, EndCoreSnp
    integer, dimension (:,:), pointer :: CoreIndex, TailIndex
    integer :: startCore, endCore, nCores, nSnp

    nSnp = phase(1,1)%getLength()

    if (.not. present(existingLibraries)) then
      CoreIndex => CalculateCores(nSnp, params%Jump, params%offset)
    else
      CoreIndex => getCoresFromLibraries(existingLibraries)
    end if
    TailIndex => calculateTails(CoreIndex, nSnp, params%Jump, params%CoreAndTailLength)
    nCores = size(CoreIndex,1)

    if (params%startCoreChar .eq. "Combine") then
      startCore = nCores + 1
    else
      read(params%startCoreChar, '(i10)') startCore
    end if

    if (params%endCoreChar .eq. "Combine") then
      endCore = nCores
    else
      read(params%endCoreChar, '(i10)') endCore
    end if

    results = AlphaPhaseResults(endCore-startCore+1)

    do h = startCore, endCore
      StartCoreSnp = CoreIndex(h, 1)
      EndCoreSnp = CoreIndex(h, 2)
      StartSurrSnp = TailIndex(h, 1)
      EndSurrSnp = TailIndex(h, 2)

      c = Core(Phase,StartCoreSnp,EndCoreSnp)
      if (.not. present(existingLibraries)) then
	library = HaplotypeLibrary(c%getNCoreSnp(),500,500)
      else
	library = existingLibraries(h)
      end if
      call UpdateHapLib(c,library,params%percminpresent,params%minoverlap,params%PercGenoHaploDisagree) 

      results%libraries(h-startCore+1) = library
      results%cores(h-startCore+1) = c
      results%ids(h-startCore+1) = h
      results%startIndexes(h-startCore+1) = CoreIndex(h,1)
      results%endIndexes(h-startCore+1) = CoreIndex(h,2)
    end do
  end function createLibraries

end module AlphaPhaseFunctions

