#ifndef _WIN32

#DEFINE DASH "/"
#DEFINE MD "mkdir -p "

#else

#DEFINE DASH "\"
#DEFINE MD "md "
#endif

module InputOutput
    use ConstantModule
    implicit none

    integer, parameter, private :: lengan = 20

contains

    subroutine WriteOutResults(allCores, startIndex, endIndex, p, params)
        use PedigreeModule
        use CoreModule
        use OutputParametersModule
        use HaplotypeModule

        type(CoreType), dimension(:), intent(in) :: allCores
        integer, dimension(:), intent(in) :: startIndex, endIndex
        type(PedigreeHolder), intent(in) :: p
        type(OutputParameters), intent(in) :: params

        integer(kind=1), dimension(:), allocatable :: tempPhase

        integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp, nCores
        integer, allocatable, dimension(:) :: WorkOut
        double precision, allocatable, dimension(:) :: CoreCount
        integer(kind=1), allocatable, dimension(:) :: TempSwap
        character(len=100) :: fmt
        integer :: finalPhaseUnit,snpPhaseRateUnit,swappatmatunit,PhasingYieldunit,finaHapIndCarryUnit,indivPhaseRateUnit
        type(Haplotype), pointer :: hap1, hap2

        nAnisG = allCores(1)%getNAnisG()
        nCores = size(allCores)
        nSnp = 0
        do i = 1, nCores
            nSnp = nSnp + allCores(i)%getNCoreSnp()
        end do

        if (params%outputFinalPhase) then
            open (newUnit = finalPhaseUnit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalPhase.txt", status = "unknown")
            allocate(tempPhase(nSnp))
            write(fmt, '(a,i10,a)') '(a20,', nSnp, 'i2)'
            do i = 1, nAnisG
                do j = 1, nCores
                    hap1 => allCores(j)%phase(i,1)
                    TempPhase(startIndex(j):endIndex(j)) = hap1%toIntegerArray()
                end do
                write(finalPhaseUnit, fmt) p%pedigree(p%hdMap(i))%originalId, &
                TempPhase
                do j = 1, nCores
                    hap2 => allCores(j)%phase(i,2)
                    TempPhase(startIndex(j):endIndex(j)) = hap2%toIntegerArray()
                end do
                write(finalPhaseUnit, fmt) p%pedigree(p%hdMap(i))%originalId, &
                TempPhase
            end do
            deallocate(tempPhase)
            close(finalPhaseUnit)
        end if

        if (params%outputCoreIndex) then      
            call writeCoreIndex(params, nCores, nAnisG, nSnp, startIndex, endIndex)
        end if

        if (params%outputSnpPhaseRate) then
            open (newUnit = snpPhaseRateUnit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"SnpPhaseRate.txt", status = "unknown")
            do i = 1, nCores
                do j = 1, allCores(i)%getNCoreSnp()
                    counter = 0
                    do k = 1, nAnisG
                        hap1 => allCores(i)%phase(k,1)
                        hap2 => allCores(i)%phase(k,2)
                        if ((hap1%getPhase(j) == 0).or.(hap1%getPhase(j) == 1)) counter = counter + 1
                        if ((hap2%getPhase(j) == 0).or.(hap2%getPhase(j) == 1)) counter = counter + 1
                    end do
                    write (snpPhaseRateUnit, '(i10,f7.2)') startIndex(i) + j - 1, (100 * (float(counter)/(2 * nAnisG)))
                end do
            end do
            close(snpPhaseRateUnit)
        end if

        if (params%outputIndivPhaseRate) then
            open (newunit = indivPhaseRateunit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"IndivPhaseRate.txt", status = "unknown")
            allocate(CoreCount(nCores * 2))
            write(fmt, '(a,i10,a)') '(a20,', nCores*2, 'f7.2)'
            do i = 1, nAnisG
                l = 0
                do j = 1, nCores
                    hap1 => allCores(j)%phase(i,1)
                    CounterP = allCores(j)%getNCoreSnp() - hap1%numberMissing()
                    l = l + 1
                    CoreCount(l) = (float(counterP)/allCores(j)%getNCoreSnp()) * 100
                    hap2 => allCores(j)%phase(i,1)
                    CounterM = allCores(j)%getNCoreSnp() - hap2%numberMissing()
                    l = l + 1
                    CoreCount(l) = (float(counterM)/allCores(j)%getNCoreSnp()) * 100
                end do
                write (indivPhaseRateunit, fmt) p%pedigree(p%hdMap(i))%originalId, CoreCount(:)
            end do
            deallocate(CoreCount)
            close(indivPhaseRateunit)
        end if

        if (params%outputHapIndex) then
            open (newunit = finaHapIndCarryUnit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalHapIndCarry.txt", status = "unknown")
            allocate(WorkOut(nCores * 2))
            write(fmt, '(a,i10,a)') '(a20,', nCores*2, 'i8)'
            do i = 1, nAnisG
                k = 0
                do j = 1, nCores
                    k = k + 2
                    WorkOut(k - 1) = AllCores(j)%hapAnis(i,1)
                    WorkOut(k) = AllCores(j)%hapAnis(i,2)
                end do
                write (finaHapIndCarryUnit, fmt) p%pedigree(p%hdMap(i))%originalId, WorkOut(:)
            end do
            deallocate(WorkOut)
            close(finaHapIndCarryUnit)
        end if

        if (params%outputSwappable) then
            open (newunit = swappatmatunit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"SwapPatMat.txt", status = "unknown")
            allocate(TempSwap(nCores))
            write(fmt, '(a,i10,a)') '(a20,', nCores, 'i2)'
            do i = 1, nAnisG
                do j = 1, nCores
                    TempSwap(j) = AllCores(j)%swappable(i)
                end do
                write (swappatmatunit, fmt) p%pedigree(p%hdMap(i))%originalId, TempSwap
            end do
            deallocate(TempSwap)
            close(swappatmatunit)
        end if

        if (params%outputPhasingYield) then
            open (newunit = PhasingYieldunit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"PhasingYield.txt", status = "unknown")
            do j = 1,nCores
                write (PhasingYieldunit, '(i10,f7.2)') j, AllCores(j)%getTotalYield()
            end do
            close(PhasingYieldunit)
        end if

    end subroutine WriteOutResults

    subroutine WriteCoreIndex(params, nCores, nAnisG, nSnp, startIndex, endIndex)
        use OutputParametersModule

        type(OutputParameters) :: params
        integer, intent(in) :: nCores, nAnisG, nSnp
        integer, dimension(:), intent(in) :: startIndex, endIndex
        integer :: unit
        integer :: i
        open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"information.txt", status = "unknown")
        write (unit, *) "numCores",nCores
        write (unit, *) "numAnimals", nAnisG
        write (unit, *) "numSnp", nSnp
        close(unit)
        open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"CoreIndex.txt", status = "unknown")
        do i = 1, nCores
            write (unit, *) i, startIndex(i), endIndex(i)
        end do
        close(unit)
    end subroutine WriteCoreIndex

    subroutine writeOutCore(c, coreID, coreStart, p, params)
        use PedigreeModule
        use CoreModule
        use OutputParametersModule
        use HaplotypeModule

        type(CoreType), intent(in) :: c
        integer, intent(in) :: coreID
        integer, intent(in) :: coreStart
        type(PedigreeHolder), intent(in) :: p
        class(OutputParameters) :: params

        integer :: i, j, counter, CounterM, CounterP, nAnisG, nSnp
        integer, allocatable, dimension(:) :: WorkOut
        double precision, allocatable, dimension(:) :: CoreCount
        character(len=100) :: fmt

        integer :: unit
        character(:), allocatable :: coreIDtxt

        type(Haplotype), pointer :: hap1, hap2

        nAnisG = c%getNAnisG()
        nSnp = c%getNCoreSnp()

        allocate(WorkOut(2))
        allocate(CoreCount(2))

        coreIDtxt = itoa(coreID)

        if (params%outputFinalPhase) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalPhase" // coreIDtxt // ".txt", status = "unknown")
            write(fmt, '(a,i10,a)') '(a20,', c%getNSnp(), 'i2)'
            do i = 1, nAnisG
                hap1 => c%phase(i, 1)
                hap2 => c%phase(i, 2)
                write(unit, fmt) p%pedigree(p%hdMap(i))%originalId, &
                hap1%toIntegerArray()
                write(unit, fmt) p%pedigree(p%hdMap(i))%originalId, &
                hap2%toIntegerArray()
            end do
            close(unit)
        end if

        if (params%outputSnpPhaseRate) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
            do i = 1, nSnp
                counter = 0
                do j = 1, nAnisG
                    hap1 => c%phase(j, 1)
                    hap2 => c%phase(j, 2)
                    if ((hap1%getPhase(i) == 0).or.(hap1%getPhase(i) == 1)) counter = counter + 1
                    if ((hap2%getPhase(i) == 0).or.(hap2%getPhase(i) == 1)) counter = counter + 1
                end do
                write (unit, '(i10,f7.2)') i + CoreStart - 1, (100 * (float(counter)/(2 * nAnisG)))
            end do
            close(unit)
        end if

        if (params%outputIndivPhaseRate) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
            do i = 1, nAnisG
                hap1 => c%phase(i, 1)
                hap2 => c%phase(i, 2)
                CounterP = c%getNCoreSnp() - hap1%numberMissing()
                CounterM = c%getNCoreSnp() - hap2%numberMissing()
                CoreCount(1) = (float(counterP)/(nSnp) * 100)
                CoreCount(2) = (float(counterM)/(nSnp) * 100)
                write(unit, '(a20,2f7.2)') p%pedigree(p%hdMap(i))%originalId, CoreCount(:)
            end do
            close(unit)
        end if

        if (params%outputHapIndex) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
            do i = 1, nAnisG
                WorkOut(1) = c%hapAnis(i, 1)
                WorkOut(2) = c%hapAnis(i, 2)
                write (unit, '(a20,2i8)') p%pedigree(p%hdMap(i))%originalId, WorkOut(:)
            end do
            close(unit)
        end if

        if (params%outputSwappable) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"SwapPatMat" // coreIDtxt // ".txt", status = "unknown")
            do i = 1, nAnisG
                write (unit, '(a20,i2)') p%pedigree(p%hdMap(i))%originalId, c%swappable(i)
            end do
            close(unit)
        end if
    end subroutine writeOutCore

    function itoa(i) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: i
        character(range(i)+2) :: tmp
        write(tmp,'(i0)') i
        res = trim(tmp)
    end function

    function ParsePedigreeAndGenotypeData(params) result(p)
        use ProgramParametersModule
        use PedigreeModule
        use SortingModule
        use PedigreeModule

        type(ProgramParameters), intent(inout) :: params
        type(PedigreeHolder) :: p

        integer :: nAnisG
        integer :: nAnisRawPedigree

        call CountInData(nAnisRawPedigree, nAnisG, params)

        if (trim(params%PedigreeFile) /= "NoPedigree") then
            call initPedigreeGenotypeFiles(p,params%GenotypeFile,nAnisG,params%nSnp,params%GenotypeFileFormat , params%PedigreeFile )
        else
            call initPedigreeGenotypeFiles(p,params%GenotypeFile,nAnisG,params%nSnp, params%GenotypeFileFormat )
        endif
    end function ParsePedigreeAndGenotypeData

    function ParsePhaseData(PhaseFile, nAnisG, nSnp) result(Phase)
        use HaplotypeModule

        character(len=300) :: PhaseFile
        integer, intent(in) :: nAnisG, nSnp
        type(Haplotype), pointer, dimension(:,:) :: Phase

        integer :: i,unit
        integer(kind=1), dimension (nSnp) :: ReadingVector
        character(lengan) :: dummy

        open (newunit = unit, file = trim(PhaseFile), status = "old")


        allocate(Phase(nAnisG, 2))

        do i = 1, nAnisG
            read (unit, *) dummy, ReadingVector(:)
            call Phase(i,1)%newHaplotypeInt(ReadingVector)
            read (unit, *) dummy, ReadingVector(:)
            call Phase(i,2)%newHaplotypeInt(ReadingVector)
        enddo

        close(unit)
    end function ParsePhaseData

    function ReadInParameterFile(filename) result (params)
        use ProgramParametersModule
        use AlphaHouseMod, only: parseToFirstWhitespace,splitLineIntoTwoParts,toLower
        character(*), intent(in) :: filename
        type(ProgramParameters) :: params

        double precision :: PercSurrDisagree
        integer :: TempInt, status
        integer :: unit
        character (len = 300) :: OffsetVariable, outputoption
        character(len=300) :: first, line
        character(len=:), allocatable::tag
        character(len=300),dimension(:),allocatable :: second

        logical :: singleRun, singleSurrogates, minOverlapSet

        params = ProgramParameters()

        open(newunit=unit, file=filename, action="read", status="old")

        status = 0
        READFILE: do while (status==0)
            read(unit,"(A)", IOStat=status)  line
            if (len_trim(line)==0) then
                CYCLE
            end if

            call splitLineIntoTwoParts(trim(line), first, second)
            tag = parseToFirstWhitespace(first)
            if (first(1:1)=="=" .or. len(trim(line))==0) then
                cycle
            else
                select case(trim(tag))

                case("pedigreefile")
                    if (.not. allocated(second)) then
                        write(*, "(A,A)") "No pedigree file specified. Using default filename: ", params%PedigreeFile
                    else
                        write(params%PedigreeFile, "(A)") second(1)
                    end if

                case("genotypefile")
                    if (.not. allocated(second)) then
                        write(*, "(A,A)") "No genotype file specified. Using default filename: ", params%Genotypefile
                    else
                        write(params%Genotypefile, "(A)") second(1)
                        if (size(second) >1) then
                            select case(trim(toLower(second(2))))
                            case('genotypeformat')
                                params%GenotypeFileFormat = 1
                            case('phaseformat')
                                params%GenotypeFileFormat = 2
                            case('unorderedformat')
                                params%GenotypeFileFormat = 3
                            end select
                        else
                            params%GenotypeFileFormat = 1
                        endif

                    endif

                case("numberofsnp")
                    read(second(1), *) params%nSnp

                case("generalcoreandtaillength")
                    read(second(1), *) params%params%CoreAndTailLength
                    write(*,*) "GeneralCoreAndTailLength is likely to be deprecated in a future release."
                    write(*,*) "Please consider using TailLength instead."

                case("taillength")
                    read(second(1), *) params%params%TailLength

                case("generalcorelength")

                    if (size(second) <2) then
                        write(*,*) "GeneralCoreLength incorrectly specified"
                    endif
                    read (second(1), *) params%params%Jump
                    read (second(2), *)OffsetVariable

                case("usethisnumberofsurrogates")
                    read(second(1), *) params%params%UseSurrsN

                case("percentagesurrdisagree")
                    read(second(1), *) PercSurrDisagree
                    PercSurrDisagree = PercSurrDisagree / 100


                case("percentagegenohaplodisagree")
                    read(second(1),*) params%params%PercGenoHaploDisagree
                    params%params%PercGenoHaploDisagree = params%params%PercGenoHaploDisagree / 100

                case("genotypemissingerrorpercentage")
                    read(second(1), *) params%params%GenotypeMissingErrorPercentage
                    params%params%GenotypeMissingErrorPercentage = params%params%GenotypeMissingErrorPercentage / 100


                case("fulloutput")
                    read(second(1),*) outputoption


                case("simulation")
                    read(second(1), *) tempInt
                    params%Simulation = (TempInt == 1)


                case("truephasefile")
                    write(params%TruePhaseFile, "(A)") second(1)


                case("iteratemethod")
                    if (allocated(second)) then
                        read(second(1), *) params%params%iterateType
                    endif

                case("iteratesubsetsize")
                    if(allocated(second)) then
                        read(second(1), *) params%params%iterateNumber
                    end if

                case("iterateiterations")
                    if (allocated(second)) then
                        read(second(1), *) params%params%numIter
                    end if
                    if (params%params%iterateType .eq. "Off") then
                        params%params%numIter = 1
                    end if


                case("cores")
                    if (size(second) == 2) then
                        read(second(1), *) params%params%startCoreChar
                        read(second(2), *) params%params%endCoreChar
                    endif

                case("minhapfreq")
                    if(allocated(second)) then
                        read (second(1), *) params%params%minHapFreq
                    end if

                case("library")
                    if (allocated(second)) then
                        read(second(1), *) params%library
                    endif

                case("minoverlap")
                    if(allocated(second)) then
                        read(second(1), *) params%params%minoverlap
                        minOverLapSet = .true.
                    end if

                case("percminpresent")
                    if(allocated(second)) then
                        read(second(1),*) params%params%percMinPresent
                        params%params%PercMinPresent = params%params%PercMinPresent / 100
                    end if

                case("percmintokeep")
                    if (allocated(second)) then
                        read(second(1),*) params%params%percMinToKeep
                        params%params%PercMinToKeep = params%params%PercMinToKeep / 100
                    end if  

                case("outputdirectory")
                    if (allocated(second)) then
                        read(second(1), *) params%outputParams%outputDirectory
                    end if

                case("corefile")
                    if (allocated(second)) then
                        read(second(1), *) params%CoreFile
                    end if

                case("presphasefile")
                    if (allocated(second)) then
                        read(second(1), *) params%PrePhaseFile
                    end if

                case("graphics")
                    write(*,"(A)") "graphics is no longer a valid option for the AlphaPhase Spec File."

                case("coreattime")
                    write(*,"(A)") "coreattime is no longer a valid option for the AlphaPhase Spec File."

                case("nrmthresh")
                    write(*,"(A)") "nrmthresh is no longer a valid option for the AlphaPhase Spec File."

                case default
                    write(*,"(A,A)") trim(tag), " is not valid for the AlphaPhase Spec File."
                    cycle
                end select
            end if
        end do READFILE


        close (unit)

        print *
        print *, " Parameter file read"
        print *, " "
        print *, " Using ", trim(params%PedigreeFile), " as the pedigree file"
        print *, " Using ", trim(params%GenotypeFile), " as the genotype file"
        print *, " "

        if (params%params%CoreAndTailLength > params%nSnp) then
            print*, "GeneralCoreAndTailLength is too long"
            stop
        endif
        if (params%params%Jump > params%nSnp) then
            print*, "GeneralCoreLength is too long"
            stop
        endif
        if ((params%params%CoreAndTailLength /= -1) .and. (params%params%CoreAndTailLength < params%params%Jump)) then
            print *, "GeneralCoreAndTailLength is shorted than GenralCoreLength"
            stop
        end if

        if (OffsetVariable == "Offset") then
            params%params%Offset = .true.
        endif
        if (OffsetVariable == "NotOffset") then
            params%params%Offset = .false.
        endif

        if ((OffsetVariable /= "Offset").and.(OffsetVariable /= "NotOffset")) then
            print*, "Offset variable is not properly parameterised"
            stop
        endif

        if ((params%params%percminpresent /= 1) .and. (params%params%percgenohaplodisagree /= 0)) then
            print *, "Not recommended to run MultiHD options (percminpresent not equal to 100%)"
            print *, "with PercentageGenoHaploDisagree set to non-zero.  AlphaPhase may be"
            print *, "SLOW."
        end if

        if ((params%params%percminpresent /= 1) .and. (params%params%minHapFreq > 1)) then
            print *, "Running MultiHD options (percminpresent not equal to 100%) with"
            print *, "with MinHapFreq greater than one is not tested or supported."
            print *, "Use at own risk!"
        end if

        if (outputoption .eq. "Impute") then
            params%outputParams%outputFinalPhase = .true.
            params%outputParams%outputCoreIndex = .true.
            params%outputParams%outputSnpPhaseRate = .false.
            params%outputParams%outputIndivPhaseRate = .false.
            params%outputParams%outputHapIndex = .false.
            params%outputParams%outputSwappable = .false.
            params%outputParams%outputHapCommonality = .false.
            params%outputParams%outputSurrogates = .false.
            params%outputParams%outputSurrogatesSummary = .false.
            params%outputParams%outputHaplotypeLibraryText = .false.
            params%outputParams%outputHaplotypeLibraryBinary = .true.
            params%outputParams%outputPhasingYield = .false.
            params%outputParams%outputTimer = .false.
            params%outputParams%outputIndivMistakes = .false.
            params%outputParams%outputIndivMistakesPercent = .false.
            params%outputParams%outputCoreMistakesPercent = .false.
            params%outputParams%outputGlobalCoreMistakesPercent = .false.
        end if
        if (outputoption .eq. "SeqOpt") then
            params%outputParams%outputFinalPhase = .false. 
            params%outputParams%outputCoreIndex = .false. 
            params%outputParams%outputSnpPhaseRate = .false. 
            params%outputParams%outputIndivPhaseRate = .false. 
            params%outputParams%outputHapIndex = .true. 
            params%outputParams%outputSwappable = .false. 
            params%outputParams%outputHapCommonality = .false. 
            params%outputParams%outputSurrogates = .false. 
            params%outputParams%outputSurrogatesSummary = .false. 
            params%outputParams%outputHaplotypeLibraryText = .true. 
            params%outputParams%outputHaplotypeLibraryBinary = .false. 
            params%outputParams%outputPhasingYield = .false. 
            params%outputParams%outputTimer = .false. 
            params%outputParams%outputIndivMistakes = .false. 
            params%outputParams%outputIndivMistakesPercent = .false. 
            params%outputParams%outputCoreMistakesPercent = .false. 
            params%outputParams%outputGlobalCoreMistakesPercent = .false.
        end if
        if ((outputoption .eq. "Full") .or. (outputoption .eq. "1")) then
            params%outputParams%outputFinalPhase = .true.
            params%outputParams%outputCoreIndex = .true.
            params%outputParams%outputSnpPhaseRate = .true.
            params%outputParams%outputIndivPhaseRate = .true.
            params%outputParams%outputHapIndex = .true.
            params%outputParams%outputSwappable = .true.
            params%outputParams%outputHapCommonality = .true.
            params%outputParams%outputSurrogates = .true.
            params%outputParams%outputSurrogatesSummary = .true.
            params%outputParams%outputHaplotypeLibraryText = .true.
            params%outputParams%outputHaplotypeLibraryBinary = .true.
            params%outputParams%outputPhasingYield = .true.
            params%outputParams%outputTimer = .true.
            params%outputParams%outputIndivMistakes = .true.
            params%outputParams%outputIndivMistakesPercent = .true.
            params%outputParams%outputCoreMistakesPercent = .true.
            params%outputParams%outputGlobalCoreMistakesPercent = .true.
        end if
        if ((outputoption .eq. "Standard") .or. (outputoption .eq. "0")) then
            params%outputParams%outputFinalPhase = .true.
            params%outputParams%outputCoreIndex = .true.
            params%outputParams%outputSnpPhaseRate = .true.
            params%outputParams%outputIndivPhaseRate = .true.
            params%outputParams%outputHapIndex = .true.
            params%outputParams%outputSwappable = .true.
            params%outputParams%outputHapCommonality = .false.
            params%outputParams%outputSurrogates = .false.
            params%outputParams%outputSurrogatesSummary = .false.
            params%outputParams%outputHaplotypeLibraryText = .false.
            params%outputParams%outputHaplotypeLibraryBinary = .true.
            params%outputParams%outputPhasingYield = .true.
            params%outputParams%outputTimer = .true.
            params%outputParams%outputIndivMistakes = .false.
            params%outputParams%outputIndivMistakesPercent = .false.
            params%outputParams%outputCoreMistakesPercent = .false.
            params%outputParams%outputGlobalCoreMistakesPercent = .false.
        end if

        singleRun = (params%params%StartCoreChar .eq. "0")
        ! Should probably have an input option to always output per core even when a single run - hence two lines here - but not
        ! currently implemented
        params%outputParams%outputPerCore = .not. singleRun
        params%outputParams%outputCombined = singleRun

        ! No purpose is served in writing out swappable info if data is prephased
        params%outputParams%outputSwappable = params%outputParams%outputSwappable .and. (.not. (params%GenotypeFileFormat /= 2))

        ! Only write out surrogate information if we only do one lot of surrogate calculations per core
        singleSurrogates = (params%params%iterateType .eq. "Off") .and. (params%params%numIter == 1)
        params%outputParams%outputSurrogates = params%outputParams%outputSurrogates .and. singleSurrogates
        params%outputParams%outputSurrogatesSummary = params%outputParams%outputSurrogatesSummary .and. singleSurrogates

        !Only write "simulation" results if we run in that mode
        params%outputParams%outputIndivMistakes = params%outputParams%outputIndivMistakes .and. params%Simulation
        params%outputParams%outputIndivMistakesPercent = params%outputParams%outputIndivMistakesPercent .and. params%Simulation
        params%outputParams%outputCoreMistakesPercent = params%outputParams%outputCoreMistakesPercent .and. params%Simulation
        params%outputParams%outputGlobalCoreMistakesPercent = params%outputParams%outputGlobalCoreMistakesPercent .and. params%Simulation

        params%params%NumSurrDisagree = int(params%params%UseSurrsN * PercSurrDisagree)

        !  if (.not. minOverlapSet) then 
        !    params%params%minOverlap = int(params%params%Jump * 0.4)
        !  end if

    end function ReadInParameterFile


    subroutine HapCommonality(library, OutputPoint, params)
        use HaplotypeLibraryModule
        use OutputParametersModule

        type(HaplotypeLibrary), intent(in) :: library
        integer, intent(in) :: OutputPoint
        class(OutputParameters), intent(in) :: params

        integer :: i, SizeCore, nHaps
        character(len = 300) :: filout
        character(len = 100) :: fmt

        integer, allocatable, dimension (:,:) :: HapRel
        integer :: HapCommonalityUnit
        SizeCore = library%getNumSnps()
        nHaps = library%getSize()

        write (filout, '(a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"Extras",a1,"HapCommonality",i0,".txt")') DASH, DASH, DASH, DASH, OutputPoint
        open (newunit = HapCommonalityUnit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')

        HapRel = library%getHapRel()
        write (fmt, '(a,i10,a)') '(i10,', size(HapRel,2), 'f5.2)'

        do i = 1, nHaps
            write (HapCommonalityUnit, fmt) i, float(HapRel(i,:))/SizeCore
        enddo

        close(HapCommonalityUnit)

    end subroutine HapCommonality

    subroutine WriteSurrogates(definition, OutputPoint, p, params)
        use SurrogateModule
        use PedigreeModule
        use OutputParametersModule

        character(len = 300) :: filout
        integer :: i, j, nSurrogates

        type(Surrogate), intent(in) :: definition
        integer, intent(in) :: OutputPoint
        type(PedigreeHolder), intent(in) :: p
        type(OutputParameters), intent(in) :: params

        integer :: nAnisG
        integer :: unit
        character(len=100) :: fmt

        nAnisG = size(definition%numOppose,1)

        if (params%outputSurrogates) then
            write (filout, '(a1,"Miscellaneous",a1,"Surrogates",i0,".txt")') DASH, DASH, OutputPoint
            open (newunit = unit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')
            write(fmt, '(a,i10,a)') '(a20,', size(definition%partition,2), 'i6)'
            do i = 1, nAnisG
                write (unit, fmt) p%pedigree(p%hdMap(i))%originalId, definition%partition(i,:)
            end do
            close(unit)
        end if


        if (params%outputSurrogatesSummary) then
            write (filout, '(a1,"Miscellaneous",a1,"SurrogatesSummary",i0,".txt")') DASH, DASH, OutputPoint
            open (newunit = unit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')
            do i = 1, nAnisG
                nSurrogates = 0
                do j = 1, nAnisG
                    if ((definition%numOppose(i, j) <= definition%threshold) .and. &
                        definition%enoughIncommon(i,j)) then
                        nSurrogates = nSurrogates + 1
                    end if
                enddo
                write (unit, '(a20,5i8)') &
                p%pedigree(p%hdMap(i))%originalId, count(definition%partition(i,:) == 1), count(definition%partition(i,:) == 2)&
                , count(definition%partition(i,:) == 3), nSurrogates, definition%method(i)
            enddo
            close(unit)
        end if

    end subroutine WriteSurrogates

    subroutine CountInData(nAnisRawPedigree, nAnisG, params)
        use ProgramParametersModule

        type(ProgramParameters), intent(in) :: params
        integer, intent(out) :: nAnisRawPedigree, nAnisG

        integer :: k,unit
        character (len = 300) :: dumC

        nAnisRawPedigree = 0
        if (trim(params%PedigreeFile) /= "NoPedigree") then
            open (newunit = unit, file = trim(params%PedigreeFile), status = "old")
            do
                read (unit, *, iostat = k) dumC
                nAnisRawPedigree = nAnisRawPedigree + 1
                if (k /= 0) then
                    nAnisRawPedigree = nAnisRawPedigree - 1
                    exit
                endif
            enddo
            close(unit)
            print*, " ", nAnisRawPedigree, " individuals in the pedigree file"
        endif

        nAnisG = 0
        open (newunit = unit, file = trim(params%GenotypeFile), status = "old")
        do
            read (unit, *, iostat = k) dumC
            nAnisG = nAnisG + 1
            if (k /= 0) then
                nAnisG = nAnisG - 1
                exit
            endif
        enddo
        if (params%GenotypeFileFormat == 2) then
            nAnisG = nAnisG /2
        end if
        close(unit)
        print*, " ", nAnisG, " individuals in the genotype file"

        if (trim(params%PedigreeFile) == "NoPedigree") nAnisRawPedigree = nAnisG

    end subroutine CountInData

    subroutine MakeDirectories(params)
        use OutputParametersModule

        type(OutputParameters), intent(in) :: params

        print*, ""

        call system(MD // trim(params%outputDirectory) // DASH // "PhasingResults")
        call system(MD // trim(params%outputDirectory) // DASH //"PhasingResults"//DASH//"HaplotypeLibrary")
        if (params%outputHapCommonality) call system(MD // trim(params%outputDirectory) // DASH //"PhasingResults"//DASH//"HaplotypeLibrary"//DASH//"Extras")
        call system(MD // trim(params%outputDirectory) // DASH //"Miscellaneous")

        if (params%outputIndivMistakes .or. params%outputIndivMistakesPercent .or. params%outputCoreMistakesPercent) then
            call system(MD // trim(params%outputDirectory) // DASH //"Simulation")
        endif

    end subroutine MakeDirectories  



    subroutine writeTimer(hours, minutes, seconds, params)
        use OutputParametersModule
        integer, intent(in) :: hours, minutes
        real, intent(in) :: seconds
        type(OutputParameters), intent(in) :: params
        integer :: unit

        if (params%outputTimer) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"Timer.txt", status = "unknown")
            write(unit, '(A27,A7,I3,A9,I3,A9,F6.2)') "Time Elapsed", "Hours", INT(Hours), "Minutes", INT(Minutes), "Seconds", Seconds
            close(unit)
        end if
    end subroutine writeTimer

    subroutine WriteTestResults(results, c, p, OutputPoint, params)
        use SurrogateModule
        use PedigreeModule
        use CoreModule
        use TestResultModule
        use OutputParametersModule

        type(TestResults), intent(in) :: results
        type(CoreType), intent(in) :: c
        type(PedigreeHolder), intent(in) :: p
        type(OutputParameters), intent(in) :: params
        integer, intent(in) :: OutputPoint

        integer :: i, unit
        integer :: nAnisG, nSnp
        character(len = 300) :: filout

        nAnisG = c % getNAnisG()
        nSNp = c % getNCoreSnp()

        if (params%outputIndivMistakes) then
            write (filout, '(a1,"Simulation",a1,"IndivMistakes",i0,".txt")') DASH, DASH, OutputPoint
            open (newunit = unit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')
            do i = 1, nAnisG
                write(unit, '(a20,6i6,a6,6i6,a6,6i6,a6,6i6)') p%pedigree(p%hdMap(i))%originalId, &
                results%counts(i,1,ALL_,CORRECT_), results%counts(i,2,ALL_,CORRECT_), &
                results%counts(i,1,ALL_,NOTPHASED_), results%counts(i,2,ALL_,NOTPHASED_), &
                results%counts(i,1,ALL_,INCORRECT_), results%counts(i,2,ALL_,INCORRECT_), "|", &
                results%counts(i,1,HET_,CORRECT_), results%counts(i,2,HET_,CORRECT_), &
                results%counts(i,1,HET_,NOTPHASED_), results%counts(i,2,HET_,NOTPHASED_), &
                results%counts(i,1,HET_,INCORRECT_), results%counts(i,2,HET_,INCORRECT_), "|", &
                results%counts(i,1,MISS_,CORRECT_), results%counts(i,2,MISS_,CORRECT_), &
                results%counts(i,1,MISS_,NOTPHASED_), results%counts(i,2,MISS_,NOTPHASED_), &
                results%counts(i,1,MISS_,INCORRECT_), results%counts(i,2,MISS_,INCORRECT_), "|", &
                results%counts(i,1,ERROR_,CORRECT_), results%counts(i,2,ERROR_,CORRECT_), &
                results%counts(i,1,ERROR_,NOTPHASED_), results%counts(i,2,ERROR_,NOTPHASED_), &
                results%counts(i,1,ERROR_,INCORRECT_), results%counts(i,2,ERROR_,INCORRECT_)
            end do    
            close(unit)
        end if

        if (params%outputIndivMistakesPercent) then
            write (filout, '(a1,"Simulation",a1,"IndivMistakesPercent",i0,".txt")') DASH, DASH, OutputPoint
            open (newunit = unit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')

            do i = 1, nAnisG
                write (unit, '(a20,6f7.1,a6,6f7.1,a6,6f7.1,a6,6f7.1)') p%pedigree(p%hdMap(i))%originalId, &
                results%percent(i,1,ALL_,CORRECT_), results%percent(i,2,ALL_,CORRECT_), &
                results%percent(i,1,ALL_,NOTPHASED_), results%percent(i,2,ALL_,NOTPHASED_), &
                results%percent(i,1,ALL_,INCORRECT_), results%percent(i,2,ALL_,INCORRECT_), "|", &
                results%percent(i,1,HET_,CORRECT_), results%percent(i,2,HET_,CORRECT_), &
                results%percent(i,1,HET_,NOTPHASED_), results%percent(i,2,HET_,NOTPHASED_), &
                results%percent(i,1,HET_,INCORRECT_), results%percent(i,2,HET_,INCORRECT_), "|", &
                results%percent(i,1,MISS_,CORRECT_), results%percent(i,2,MISS_,CORRECT_), &
                results%percent(i,1,MISS_,NOTPHASED_), results%percent(i,2,MISS_,NOTPHASED_), &
                results%percent(i,1,MISS_,INCORRECT_), results%percent(i,2,MISS_,INCORRECT_), "|", &
                results%percent(i,1,ERROR_,CORRECT_), results%percent(i,2,ERROR_,CORRECT_), &
                results%percent(i,1,ERROR_,NOTPHASED_), results%percent(i,2,ERROR_,NOTPHASED_), &
                results%percent(i,1,ERROR_,INCORRECT_), results%percent(i,2,ERROR_,INCORRECT_)
            end do
            close(unit)
        end if

        if (params%outputCoreMistakesPercent) then
            write (filout, '(a1,"Simulation",a1,"CoreMistakesPercent",i0,".txt")') DASH, DASH, OutputPoint
            open (newunit = unit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')
            write (unit, '(6f9.4)') &
            (results%percentAll(1,ALL_,CORRECT_) + results%percentAll(2,ALL_,CORRECT_)) / 2, &
            (results%percentAll(1,HET_,CORRECT_) + results%percentAll(2,HET_,CORRECT_)) / 2, &
            (results%percentAll(1,ALL_,NOTPHASED_) + results%percentAll(2,ALL_,NOTPHASED_)) / 2, &
            (results%percentAll(1,HET_,NOTPHASED_) + results%percentAll(2,HET_,NOTPHASED_)) / 2, &
            (results%percentAll(1,ALL_,INCORRECT_) + results%percentAll(2,ALL_,INCORRECT_)) / 2, &
            (results%percentAll(1,HET_,INCORRECT_) + results%percentAll(2,HET_,INCORRECT_)) / 2
            close(unit)
        end if
    end subroutine WriteTestResults

    function getHaplotypeLibraries(directory) result (libraries)
        use HaplotypeLibraryModule

        character(*), intent(in) :: directory
        type(HaplotypeLibrary), dimension(:), pointer :: libraries

        character(4096) :: filename
        integer :: nLibs, i
        logical :: ex

        i = 0
        nLibs = 0
        ex = .true.
        do while (ex)
            i = i + 1
            write (filename, '(a, a, "HapLib", i0, ".bin")') trim(directory), DASH, i
            inquire(FILE=filename,EXIST=ex)
            if (ex) then
                nLibs = nLibs + 1
            end if
        end do

        allocate(libraries(nLibs))
        do i = 1, nLibs
            write (filename, '(a, a, "HapLib", i0, ".bin")') trim(directory), DASH, i
            libraries(i) = HaplotypeLibrary(filename,500)
        end do
    end function getHaplotypeLibraries

    function getCoresFromHapLib(directory) result (CoreIndex)
        character(*), intent(in) :: directory
        integer, dimension(:,:), pointer :: CoreIndex

        integer :: i, numLibraries, start, nHaps, nSnps
        character(4096) :: filename
        logical :: ex

        i = 1
        write (filename, '(a, a, "HapLib", i0, ".bin")') trim(directory), DASH, i 
        inquire(FILE=filename, EXIST=ex)
        do while (ex)
            i = i + 1
            write (filename, '(a, a, "HapLib", i0, ".bin")') trim(directory), DASH, i
            inquire(FILE=filename, EXIST=ex)
        end do
        numLibraries = i - 1

        allocate(CoreIndex(numLibraries,2))
        start = 1

        do i = 1, numLibraries
            write (filename, '(a, a, "HapLib", i0, ".bin")') trim(directory), DASH, i
            open (unit=2001,file=trim(filename),status="old",form="unformatted")
            read(2001) nHaps,nSnps
            close(2001)

            CoreIndex(i,1) = start
            CoreIndex(i,2) = start + nSnps - 1

            start = start + nSnps
        end do
    end function getCoresFromHapLib

    subroutine Header
        print*, ""
        print*, "                              ************************                         "
        print*, "                              *                      *                         "
        print*, "                              *   AlphaPhase 1.3.6   *                         "
        print*, "                              *                      *                         "
        print*, "                              ************************                         "
        print*, "                                                                              "
        print*, "                    Software For Phasing and Imputing Genotypes               "
    end subroutine Header

    subroutine Titles

        call Header
        print*, ""
        print*, ""

    end subroutine Titles

    subroutine PrintTimerTitles(params)
        use OutputParametersModule

        implicit none

        class(OutputParameters), intent(in) :: params

        real :: etime ! Declare the type of etime()
        real :: elapsed(2) ! For receiving user and system time
        real :: total, Minutes, Hours, Seconds

        print*, ""
        print*, ""
        call Header
        PRINT*, ""
        PRINT*, "                                  No Liability"
        PRINT*, ""
        PRINT*, "                Analysis Finished                         "

        total = etime(elapsed)
        Minutes = total/60
        Seconds = Total - (INT(Minutes) * 60)
        Hours = Minutes/60
        Minutes = INT(Minutes)-(INT(Hours) * 60)

        PRINT '(A107,A7,I3,A9,I3,A9,F6.2)', "Time Elapsed", "Hours", INT(Hours), "Minutes", INT(Minutes), "Seconds", Seconds
        call writeTimer(INT(Hours),INT(Minutes),Seconds,params)
    end subroutine PrintTimerTitles

    subroutine WriteHapLib(library, currentcore, params)
        use OutputParametersModule
        use CoreModule
        use HaplotypeModule
        use HaplotypeLibraryModule

        class(HaplotypeLibrary), intent(in) :: library
        integer, intent(in) :: currentcore
        type(OutputParameters), intent(in) :: params

        integer :: i, SizeCore, nHaps
        integer :: haplibUnitBin, haplibunit
        character(len = 300) :: filout
        type(Haplotype) :: hap
        character(len=100) :: fmt

        SizeCore = library%getNumSnps()

        nHaps = library%getSize()

        if (params%outputHaplotypeLibraryText) then
            write (filout, '(a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"HapLib",i0,".txt")') DASH, DASH, DASH, currentcore
            open (newunit = haplibunit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')
            write (haplibunit,*) nHaps, SizeCore
        endif
        if (params%outputHaplotypeLibraryBinary) then
            write (filout, '(a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"HapLib",i0,".bin")') DASH, DASH, DASH, currentcore
            open (newunit = haplibunitbin, FILE = trim(params%outputDirectory)//filout, form = "unformatted", status = 'unknown')

            write (haplibunitbin) nHaps, SizeCore
        end if

        write(fmt, '(a,i10,a)') '(2i10,a2,', library%nSnps, 'i1)'

        do i = 1, nHaps
            hap = library%getHap(i)
            if (params%outputHaplotypeLibraryText) then   
                write (haplibunit, fmt) i, library%getHapFreq(i), " ", hap%toIntegerArray()

            end if
            if (params%outputHaplotypeLibraryBinary) then
                write (haplibunitbin) hap%toIntegerArray()
            end if
        end do

        if (params%outputHaplotypeLibraryText) then
            close(haplibunit)
        end if
        if (params%outputHaplotypeLibraryBinary) then
            close(haplibunitbin)
        end if
    end subroutine WriteHapLib

    subroutine writeAlphaPhaseResults(results,p,params)
        use AlphaPhaseResultsModule
        use PedigreeModule
        use OutputParametersModule

        class(AlphaPhaseResults), intent(in) :: results
        class(PedigreeHolder), intent(in) :: p
        class(OutputParameters), intent(in) :: params

        integer :: i, id

        call MakeDirectories(params)
        ! write out core index - 
        do i = 1, results%nCores
            id = results%ids(i)
            call WriteHapLib(results%libraries(i), id, params)
            if (params%outputHapCommonality) then
                call HapCommonality(results%libraries(i), id, params)
            end if
            if (params%outputPerCore) then
                call WriteOutCore(results%cores(i), id, results%startIndexes(i), p, params)
            end if
            if (params%outputSurrogates .or. params%outputSurrogatesSummary) then
                call writeSurrogates(results%surrogates(i), id, p, params)
            end if
            if (params%outputCombined) then
                call WriteOutResults(results%cores,results%startIndexes,results%endIndexes,p,params)
            end if
            if (params%outputIndivMistakes .or. params%outputIndivMistakesPercent .or. params%outputCoreMistakesPercent) then
                call writeTestResults(results%testResults(i), results%cores(i), p, id, params)
            end if
        enddo
        if (params%outputGlobalCoreMistakesPercent) then
            call makeCoreMistakes(params, results%nCores)
        end if
    end subroutine writeAlphaPhaseResults  

    subroutine makeCoreMistakes(params, nCores)
        use OutputParametersModule

        type(OutputParameters), intent(in) :: params
        integer, intent(in) :: nCores

        integer :: i
        double precision, dimension(6) :: single, sums
        character(len = 4096) :: filin, filout
        integer :: coreMissUnit, inunit
        write (filout, '(a1,"Simulation",a1,"CoreMistakesPercent.txt")') DASH, DASH
        open (newunit = coreMissUnit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')

        do i = 1, nCores
            write (filin, '(a1,"Simulation",a1,"CoreMistakesPercent",i0,".txt")') DASH, DASH, i
            open (newunit = inunit, FILE = trim(params%outputDirectory)//filin, status = 'unknown') 
            read(inunit,'(6f9.4)') single
            close(inunit)

            sums = sums + single

            write (coreMissUnit, '(6f9.4)') single
        end do

        write(coreMissUnit,*)    
        write (coreMissUnit, '(6f9.4)') sums / nCores

        close(coreMissUnit)
    end subroutine makeCoreMistakes

    subroutine readHapLib(library, currentcore, params)
        use OutputParametersModule
        use CoreModule
        use HaplotypeModule
        use HaplotypeLibraryModule

        type(HaplotypeLibrary), intent(inout) :: library
        type(HaplotypeLibrary) :: tmp
        integer, intent(in) :: currentcore
        type(OutputParameters), intent(in) :: params

        integer :: i, SizeCore, nHaps,freq, dumI
        character(len = 300) :: filout
        character(len=1) :: dumC
        type(Haplotype) :: hap
        integer(kind=1), dimension(:), allocatable :: hapArray
        integer :: haplibunit, haplibunitbin
        character(len=100) :: fmt

        if (params%outputHaplotypeLibraryText) then
            write (filout, '(a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"HapLib",i0,".txt")') DASH, DASH, DASH, currentcore
            open (newunit = haplibunit, FILE = trim(params%outputDirectory)//filout, status = 'unknown')
            read (haplibunit,*) nHaps, SizeCore
        endif
        if (params%outputHaplotypeLibraryBinary) then
            write (filout, '(a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"HapLib",i0,".bin")') DASH, DASH, DASH, currentcore
            open (newunit = haplibunitbin, FILE = trim(params%outputDirectory)//filout, form = "unformatted", status = 'unknown')
            read (haplibunitbin) nHaps, SizeCore
        endif

        library = HaplotypeLibrary(sizeCore,nHaps,1)

        ! library%setSize(nHaps) 

        allocate(hapArray(SizeCore))
        do i = 1, nHaps
            

            if (params%outputHaplotypeLibraryBinary) then
                read (haplibunitbin) hapArray
                call hap%newHaplotypeInt(hapArray)
                dumI = library%addHap(hap)
            end if
            if (params%outputHaplotypeLibraryText) then
                write(fmt, '(a,i10,a)') '(2i10,a2,', sizeCore, 'i1)'
                read (haplibunit, fmt) dumI, freq,dumC, hapArray
                call library%setHapFreq(i,freq)
                call hap%newHaplotypeInt(hapArray)
                dumI = library%addHap(hap)
            end if
        end do

        if (params%outputHaplotypeLibraryText) then
            close(haplibunit)
        end if
        if (params%outputHaplotypeLibraryBinary) then
            close(haplibunitbin)
        end if
    end subroutine readHapLib




    subroutine readInResults(results, params, p)
        use PedigreeModule
        use CoreModule
        use OutputParametersModule
        use AlphaPhaseResultsModule
        use HaplotypeModule

        type(CoreType), dimension(:),allocatable :: allCores
        integer, dimension(:), allocatable :: startIndex, endIndex
        type(AlphaPhaseResults), intent(out) :: results
        type(OutputParameters), intent(in) :: params
        type(pedigreeHolder), intent(in) :: p
        integer(kind=1), dimension(:), allocatable :: tempPhase

        integer :: i, j, k, nAnisG, nSnp, nCores,dumI
        integer, allocatable, dimension(:) :: WorkOut
        integer(kind=1), allocatable, dimension(:) :: TempSwap
        character(len=100) :: fmt
        character(len=IDLENGTH) :: dumC

        type(Haplotype) :: hap1, hap2
        integer :: unit
        nSnp = 0



        if (params%outputCoreIndex) then
             open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"information.txt", status = "unknown")
            read (unit, *) dumC,nCores
            read (unit, *) dumC, nAnisG
            read (unit, *) dumC, nSnp
            close(unit)
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"CoreIndex.txt", status = "unknown")
            results = AlphaPhaseResults(nCores, .true., .true.)
            allocate(startIndex(nCores))
            allocate(endIndex(nCores))
            allocate(allCores(nCores))
            do i = 1, nCores
                read (unit, *) dumI , startIndex(i), endIndex(i)
                allcores(i) = newCore(p,startIndex(i),startIndex(i),endIndex(i),endIndex(i))
            end do
            close(unit)
        end if
        results%startIndexes = startIndex
        results%endIndexes = endIndex


        if (params%outputFinalPhase) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalPhase.txt", status = "unknown")
            allocate(tempPhase(nSnp))
            write(fmt, '(a,i10,a)') '(a20,', nSnp, 'i2)'

            do i = 1, nAnisG
                read(unit, *) dumC, TempPhase
                do j = 1, nCores
                    call hap1%newHaplotypeInt(TempPhase(startIndex(j):endIndex(j)))
                    allCores(j)%phase(i,1) = hap1
                end do
                read(unit, *) dumC, TempPhase
                do j = 1, nCores
                    call hap2%newHaplotypeInt(TempPhase(startIndex(j):endIndex(j)))
                    allCores(j)%phase(i,2) = hap2
                end do
            end do
            deallocate(tempPhase)
            close(unit)
        end if

        if (params%outputHapIndex) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalHapIndCarry.txt", status = "unknown")
            allocate(WorkOut(nCores * 2))
            write(fmt, '(a,i10,a)') '(a20,', nCores*2, 'i8)'
            do i = 1, nAnisG
                k = 0
                read (unit, *) dumC, WorkOut
                do j = 1, nCores
                    k = k + 2
                    AllCores(j)%hapAnis(i,1) = WorkOut(k - 1)
                    AllCores(j)%hapAnis(i,2) = WorkOut(k)
                end do

            end do
            deallocate(WorkOut)
            close(unit)
        end if

        if (params%outputSwappable) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"SwapPatMat.txt", status = "unknown")
            allocate(TempSwap(nCores))
            write(fmt, '(a,i10,a)') '(a20,', nCores, 'i2)'
            do i = 1, nAnisG
                read (unit, fmt) dumC, TempSwap
                do j = 1, nCores
                    AllCores(j)%swappable(i) =TempSwap(j)
                end do

            end do
            deallocate(TempSwap)
            close(unit)
        end if

        results%cores = allCores
    end subroutine readInResults

    subroutine readInPerCoreResults(results, params, p)
        use PedigreeModule
        use CoreModule
        use OutputParametersModule
        use AlphaPhaseResultsModule
        use HaplotypeModule

        type(CoreType), dimension(:),allocatable :: allCores
        integer, dimension(:), allocatable :: startIndex, endIndex
        type(AlphaPhaseResults), intent(out) :: results
        type(OutputParameters), intent(in) :: params
        type(pedigreeHolder), intent(in) :: p
        integer(kind=1), dimension(:), allocatable :: tempPhase

        integer :: i, j, nAnisG, nSnp, nCores,dumI
        integer, allocatable, dimension(:) :: WorkOut
        integer(kind=1), allocatable, dimension(:) :: TempSwap
        character(len=100) :: fmt
        character(len=IDLENGTH) :: dumC

        type(Haplotype) :: hap1, hap2
        integer :: unit
        nSnp = 0



        if (params%outputCoreIndex) then
            open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"CoreIndex.txt", status = "unknown")
            read (unit, *) nCores
            read (unit, *) nAnisG
            read (unit, *) nSnp
            results = AlphaPhaseResults(nCores, .true., .true.)
            allocate(startIndex(nCores))
            allocate(endIndex(nCores))
            allocate(allCores(nCores))
            do i = 1, nCores
                read (unit, *) dumI , startIndex(i), endIndex(i)
                allcores(i) = newCore(p,startIndex(i),startIndex(i),endIndex(i),endIndex(i))
            end do
            close(unit)
        end if
        results%startIndexes = startIndex
        results%endIndexes = endIndex


        do j = 1, nCores
            if (params%outputFinalPhase) then
                open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalPhase" // itoa(j) // ".txt", status = "unknown")

                write(fmt, '(a,i10,a)') '(a20,', endIndex(j) - startIndex(j) + 1, 'i2)'

                allocate(TempPhase(endIndex(j) - startIndex(j) + 1))
                do i = 1, nAnisG
                    read(unit, *) dumC, TempPhase
                    call hap1%newHaplotypeInt(TempPhase)
                    allCores(j)%phase(i,1) = hap1
                    read(unit, *) dumC, TempPhase
                    call hap2%newHaplotypeInt(TempPhase)
                    allCores(j)%phase(i,2) = hap2
                end do
                deallocate(tempPhase)
                close(unit)
            end if

            if (params%outputHapIndex) then
                open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"FinalHapIndCarry" // itoa(j) // ".txt", status = "unknown")
                allocate(WorkOut(2))
                write(fmt, '(a,i10,a)') '(a20,', 2, 'i8)'
                do i = 1, nAnisG
                    read (unit, *) dumC, WorkOut
                    AllCores(j)%hapAnis(i,1) = WorkOut(1)
                    AllCores(j)%hapAnis(i,2) = WorkOut(2)
                end do
                deallocate(WorkOut)
                close(unit)
            end if

            if (params%outputSwappable) then
                open (newunit = unit, file = trim(params%outputDirectory)//DASH//"PhasingResults"//DASH//"SwapPatMat" // itoa(j) // ".txt", status = "unknown")
                allocate(TempSwap(nCores))
                write(fmt, '(a,i10,a)') '(a20,', nCores, 'i2)'
                do i = 1, nAnisG
                    read (unit, fmt) dumC, AllCores(j)%swappable(i)
                end do
                deallocate(TempSwap)
                close(unit)
            end if
        end do

        results%cores = allCores
    end subroutine readInPerCoreResults

    subroutine readAlphaPhaseResults(results,params,p)
        use AlphaPhaseResultsModule
        use PedigreeModule
        use OutputParametersModule

        class(AlphaPhaseResults), intent(out) :: results
        class(OutputParameters), intent(in) :: params
        type(PedigreeHolder), intent(in) :: p
        integer :: i


        call readInResults(results,params, p)

        do i = 1, results%nCores   
            call readHapLib(results%libraries(i), i, params)
        end do


    end subroutine readAlphaPhaseResults

    subroutine printCoreInfo(coreIndex)
        integer, dimension(:,:), pointer, intent(in) :: CoreIndex

        integer :: i

        print *, "Number of cores:", size(CoreIndex,1)
        print *
        print '(3a10)', "ID", "Start", "End"
        do i = 1, size(CoreIndex,1)
            print '(3i10)', i, CoreIndex(i,1), CoreIndex(i,2)
        end do
    end subroutine printCoreInfo


#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

    subroutine PrintVersion
        call Header
        print *
        print *, "                              Commit:   "//TOSTRING(COMMIT)
        print *, "                              Compiled: "//__DATE__//", "//__TIME__  
    end subroutine PrintVersion
end module InputOutput
