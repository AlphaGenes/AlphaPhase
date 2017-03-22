#ifdef OS_UNIX

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"
#DEFINE SH "sh"
#DEFINE EXE ""
#DEFINE NULL ""

#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S /Q"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#DEFINE SH "BAT"
#DEFINE EXE ".exe"
#DEFINE NULL " >NUL"
#endif

module InputOutput
  use ConstantModule
  implicit none

  integer, parameter, private :: lengan = 20
  
contains

  subroutine WriteOutResults(allCores, startIndex, endIndex, p, writeSwappable, params)
    use PedigreeModule
    use CoreDefinition
    use ProgramParametersDefinition
    use HaplotypeModule

    type(Core), dimension(:), intent(in) :: allCores
    integer, dimension(:), intent(in) :: startIndex, endIndex
    type(PedigreeHolder), intent(in) :: p
    logical :: writeSwappable
    type(ProgramParameters), intent(in) :: params

    integer(kind=1), dimension(:), allocatable :: tempPhase

    integer :: i, j, k, l, counter, CounterM, CounterP, nAnisG, nSnp, nCores
    integer, allocatable, dimension(:) :: WorkOut
    double precision, allocatable, dimension(:) :: CoreCount
    integer(kind=1), allocatable, dimension(:) :: TempSwap
    character(len=100) :: fmt

    type(Haplotype), pointer :: hap1, hap2

    nAnisG = allCores(1)%getNAnisG()
    nCores = size(allCores)
    nSnp = 0
    do i = 1, nCores
      nSnp = nSnp + allCores(i)%getNCoreSnp()
    end do

    if (params%outputParams%outputFinalPhase) then
      open (unit = 15, file = "."//DASH//"PhasingResults"//DASH//"FinalPhase.txt", status = "unknown")
      allocate(tempPhase(nSnp))
      write(fmt, '(a,i10,a)') '(a20,', nSnp, 'i2)'
      do i = 1, nAnisG
        do j = 1, nCores
          hap1 => allCores(j)%phase(i,1)
          TempPhase(startIndex(j):endIndex(j)) = hap1%toIntegerArray()
        end do
        write(15, fmt) p%pedigree(p%hdMap(i))%originalId, &
          TempPhase
        do j = 1, nCores
          hap2 => allCores(j)%phase(i,2)
          TempPhase(startIndex(j):endIndex(j)) = hap2%toIntegerArray()
        end do
        write(15, fmt) p%pedigree(p%hdMap(i))%originalId, &
          TempPhase
      end do
      deallocate(tempPhase)
      close(15)
    end if

    if (params%outputParams%outputCoreIndex) then
      open (unit = 25, file = "."//DASH//"PhasingResults"//DASH//"CoreIndex.txt", status = "unknown")
      do i = 1, nCores
        write (25, *) i, startIndex(i), endIndex(i)
      end do
      close(25)
    end if

    if (params%outputParams%outputSnpPhaseRate) then
      open (unit = 28, file = "."//DASH//"PhasingResults"//DASH//"SnpPhaseRate.txt", status = "unknown")
      do i = 1, nCores
        do j = 1, allCores(i)%getNCoreSnp()
          hap1 => allCores(i)%phase(j,1)
          hap2 => allCores(i)%phase(j,2)
          counter = 0
          do k = 1, nAnisG
            if ((hap1%getPhaseMod(j) == 0).or.(hap1%getPhaseMod(j) == 1)) counter = counter + 1
            if ((hap2%getPhaseMod(j) == 0).or.(hap2%getPhaseMod(j) == 1)) counter = counter + 1
          end do
          write (28, '(i10,f7.2)') startIndex(i) + j - 1, (100 * (float(counter)/(2 * nAnisG)))
        end do
      end do
      close(28)
    end if

    if (params%outputParams%outputIndivPhaseRate) then
      open (unit = 30, file = "."//DASH//"PhasingResults"//DASH//"IndivPhaseRate.txt", status = "unknown")
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
        write (30, fmt) p%pedigree(p%hdMap(i))%originalId, CoreCount(:)
      end do
      deallocate(CoreCount)
      close(30)
    end if

    if (params%outputParams%outputHapIndex) then
      open (unit = 33, file = "."//DASH//"PhasingResults"//DASH//"FinalHapIndCarry.txt", status = "unknown")
      allocate(WorkOut(nCores * 2))
      write(fmt, '(a,i10,a)') '(a20,', nCores*2, 'i8)'
      do i = 1, nAnisG
        k = 0
        do j = 1, nCores
          k = k + 2
          WorkOut(k - 1) = AllCores(j)%getHapAnis(i,1)
          WorkOut(k) = AllCores(j)%getHapAnis(i,2)
        end do
        write (33, fmt) p%pedigree(p%hdMap(i))%originalId, WorkOut(:)
      end do
      deallocate(WorkOut)
      close(33)
    end if

    if ((params%outputParams%outputSwappable) .and. (writeSwappable)) then
      open (unit = 44, file = "."//DASH//"PhasingResults"//DASH//"SwapPatMat.txt", status = "unknown")
      allocate(TempSwap(nCores))
      write(fmt, '(a,i10,a)') '(a20,', nCores, 'i2)'
      do i = 1, nAnisG
        do j = 1, nCores
          TempSwap(j) = AllCores(j)%getSwappable(i)
        end do
        write (44, fmt) p%pedigree(p%hdMap(i))%originalId, TempSwap
      end do
      deallocate(TempSwap)
      close(44)
    end if
    
    if (params%outputParams%outputPhasingYield) then
      open (unit = 29, file = "."//DASH//"PhasingResults"//DASH//"PhasingYield.txt", status = "unknown")
      do j = 1,nCores
	write (29, '(i10,f7.2)') j, AllCores(j)%getTotalYield()
      end do
      close(29)
    end if
  end subroutine WriteOutResults

  subroutine writeOutCore(c, coreID, coreStart, p, writeSwappable, params)
    use PedigreeModule
    use CoreDefinition
    use ProgramParametersDefinition
    use HaplotypeModule

    type(Core), intent(in) :: c
    integer, intent(in) :: coreID
    integer, intent(in) :: coreStart
    type(PedigreeHolder), intent(in) :: p
    logical :: writeSwappable
    class(ProgramParameters) :: params

    integer :: i, j, counter, CounterM, CounterP, nAnisG, nSnp
    integer, allocatable, dimension(:) :: WorkOut
    double precision, allocatable, dimension(:) :: CoreCount
    character(len=100) :: fmt

    character(:), allocatable :: coreIDtxt

    type(Haplotype), pointer :: hap1, hap2

    nAnisG = c%getNAnisG()
    nSnp = c%getNCoreSnp()

    allocate(WorkOut(2))
    allocate(CoreCount(2))

    coreIDtxt = itoa(coreID)

    if (params%outputParams%outputFinalPhase) then
      open (unit = 15, file = "."//DASH//"PhasingResults"//DASH//"FinalPhase" // coreIDtxt // ".txt", status = "unknown")
      write(fmt, '(a,i10,a)') '(a20,', c%getNSnp(), 'i2)'
      do i = 1, nAnisG
	hap1 => c%phase(j, 1)
	hap2 => c%phase(j, 2)
	write(15, fmt) p%pedigree(p%hdMap(i))%originalId, &
	hap1%toIntegerArray()
	write(15, fmt) p%pedigree(p%hdMap(i))%originalId, &
	hap2%toIntegerArray()
      end do
      close(15)
    end if

    if (params%outputParams%outputSnpPhaseRate) then
      open (unit = 28, file = "."//DASH//"PhasingResults"//DASH//"SnpPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      do i = 1, nSnp
        counter = 0
        do j = 1, nAnisG
          hap1 => c%phase(j, 1)
          hap2 => c%phase(j, 2)
          if ((hap1%getPhaseMod(i) == 0).or.(hap1%getPhaseMod(i) == 1)) counter = counter + 1
          if ((hap2%getPhaseMod(i) == 0).or.(hap2%getPhaseMod(i) == 1)) counter = counter + 1
        end do
        write (28, '(i10,f7.2)') i + CoreStart - 1, (100 * (float(counter)/(2 * nAnisG)))
      end do
      close(28)
    end if

    if (params%outputParams%outputIndivPhaseRate) then
      open (unit = 30, file = "."//DASH//"PhasingResults"//DASH//"IndivPhaseRate" // coreIDtxt // ".txt", status = "unknown")
      do i = 1, nAnisG
        hap1 => c%phase(j, 1)
        hap2 => c%phase(j, 2)
        CounterP = c%getNCoreSnp() - hap1%numberMissing()
        CounterM = c%getNCoreSnp() - hap2%numberMissing()
        CoreCount(1) = (float(counterP)/(nSnp) * 100)
        CoreCount(2) = (float(counterM)/(nSnp) * 100)
        write(30, '(a20,2f7.2)') p%pedigree(p%hdMap(i))%originalId, CoreCount(:)
      end do
      close(30)
    end if

    if (params%outputParams%outputHapIndex) then
      open (unit = 33, file = "."//DASH//"PhasingResults"//DASH//"FinalHapIndCarry" // coreIDtxt // ".txt", status = "unknown")
      do i = 1, nAnisG
        WorkOut(1) = c%getHapAnis(i, 1)
        WorkOut(2) = c%getHapAnis(i, 2)
        write (33, '(a20,2i8)') p%pedigree(p%hdMap(i))%originalId, WorkOut(:)
      end do
      close(33)
    end if

    if ((params%outputParams%outputSwappable) .and. (writeSwappable)) then
      open (unit = 44, file = "."//DASH//"PhasingResults"//DASH//"SwapPatMat" // coreIDtxt // ".txt", status = "unknown")
      do i = 1, nAnisG
        write (44, '(a20,i2)') p%pedigree(p%hdMap(i))%originalId, c%getSwappable(i)
      end do
      close(44)
    end if
  end subroutine writeOutCore

  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function

  subroutine CombineResults(nAnisG, writeSwappable, params)
    use PedigreeModule
    use ProgramParametersDefinition
    use CoreUtils

    integer, intent(in) :: nAnisG
    logical :: writeSwappable
    type(ProgramParameters) :: params

    integer, dimension(:,:), pointer :: CoreIndex
    integer :: nCores

    integer, dimension(:), allocatable :: inUnits
    integer :: i, j, coreLength, inUnit
    character(:), allocatable :: coreIDtxt
    integer(kind=1), dimension(:), allocatable :: tempPhase
    integer, dimension(2) :: tempHapInd, tempIndivSwap
    double precision, dimension(2) :: tempIndivPhase
    double precision :: tempSnpPhase
    character(len=20) :: id
    
    if (params%library .eq. "None") then
      CoreIndex => CalculateCores(params%nSnp, params%params%Jump, params%params%offset)    
    else
      CoreIndex => getCoresFromHapLib(params%library)
    end if

    nCores = size(CoreIndex,1)

    if (params%outputParams%outputCoreIndex) then
      open (unit = 25, file = "."//DASH//"PhasingResults"//DASH//"CoreIndex.txt", status = "unknown")
      do i = 1, nCores
        write (25, *) i, CoreIndex(i,:)
      end do
      close(25)
    end if

    if (params%outputParams%outputFinalPhase) then
      !!! FINAL PHASE !!!
      open (unit = 15, file = "."//DASH//"PhasingResults"//DASH//"FinalPhase.txt", status = "unknown")
      allocate(inUnits(nCores))
      do i = 1, nCores
        coreIDtxt = itoa(i)
        open (newunit = inUnits(i), file = "."//DASH//"PhasingResults"//DASH//"FinalPhase" // coreIDtxt // ".txt", status = "old")
      end do

      do i = 1, nAnisG * 2
        do j = 1, nCores
          coreLength  = CoreIndex(j,2) - CoreIndex(j,1) + 1
          allocate(tempPhase(coreLength))
          read(inUnits(j),'(a20,' // itoA(coreLength) // 'i2)') id, tempPhase
          if (j == 1) then
            write(15, '(a20)', advance = 'no') id
          end if
          if (j == nCores) then
            write(15, '(' // itoA(coreLength) // 'i2)', advance='yes') tempPhase
          else
            write(15, '(' // itoA(coreLength) // 'i2)', advance='no') tempPhase
          end if
          deallocate(tempPhase)
        end do
      end do

      do i = 1, nCores
        close(inUnits(i))
      end do

      deallocate(inUnits)
      close(15)
    end if

    if (params%outputParams%outputHapIndex) then
      !!! HAPINDCARRY !!!
      open (unit = 33, file = "."//DASH//"PhasingResults"//DASH//"FinalHapIndCarry.txt", status = "unknown")
      allocate(inUnits(nCores))
      do i = 1, nCores
        coreIDtxt = itoa(i)
        open (newunit = inUnits(i), file = "."//DASH//"PhasingResults"//DASH//"FinalHapIndCarry" // coreIDtxt // ".txt", status = "old")
      end do

      do i = 1, nAnisG
        do j = 1, nCores
          read(inUnits(j),'(a20,2i8)') id, tempHapInd
          if (j == 1) then
            write(33, '(a20)', advance = 'no') id
          end if
          if (j == nCores) then
            write(33, '(2i8)', advance='yes') tempHapInd
          else
            write(33, '(2i8)', advance='no') tempHapInd
          end if
        end do
      end do

      do i = 1, nCores
        close(inUnits(i))
      end do

      deallocate(inUnits)
      close(33)
    end if

    if (params%outputParams%outputIndivPhaseRate) then
      !!! INDIVPHASE !!!
      open (unit = 30, file = "."//DASH//"PhasingResults"//DASH//"IndivPhaseRate.txt", status = "unknown")
      allocate(inUnits(nCores))
      do i = 1, nCores
        coreIDtxt = itoa(i)
        open (newunit = inUnits(i), file = "."//DASH//"PhasingResults"//DASH//"IndivPhaseRate" // coreIDtxt // ".txt", status = "old")
      end do

      do i = 1, nAnisG
        do j = 1, nCores
          read(inUnits(j),'(a20,2f7.2)') id, tempIndivPhase
          if (j == 1) then
            write(30, '(a20)', advance = 'no') id
          end if
          if (j == nCores) then
            write(30, '(2f7.2)', advance='yes') tempIndivPhase
          else
            write(30, '(2f7.2)', advance='no') tempIndivPhase
          end if
        end do
      end do

      do i = 1, nCores
        close(inUnits(i))
      end do

      deallocate(inUnits)
      close(30)
    end if

    if (params%outputParams%outputSnpPhaseRate) then
      !!! SNPPHASE !!!
      open (unit = 28, file = "."//DASH//"PhasingResults"//DASH//"SnpPhaseRate.txt", status = "unknown")
      do i = 1, nCores
        coreLength  = CoreIndex(i,2) - CoreIndex(i,1) + 1
        coreIDtxt = itoa(i)
        open (newunit = inUnit, file = "."//DASH//"PhasingResults"//DASH//"SnpPhaseRate" // coreIDtxt // ".txt", status = "old")
        do j = 1, coreLength
          read(inUnit,'(a10,f7.2)') id, tempSnpPhase
          write(28,'(a10,f7.2)') id, tempSnpPhase
        end do
        close(inUnit)
      end do
      close(28)
    end if

    !!! SWAPHAPMAT !!!
    if (params%outputParams%outputSwappable .and. writeSwappable) then
      open (unit = 44, file = "."//DASH//"PhasingResults"//DASH//"SwapPatMat.txt", status = "unknown")
      allocate(inUnits(nCores))
      do i = 1, nCores
        coreIDtxt = itoa(i)
        open (newunit = inUnits(i), file = "."//DASH//"PhasingResults"//DASH//"SwapPatMat" // coreIDtxt // ".txt", status = "old")
      end do

      do i = 1, nAnisG
        do j = 1, nCores
          read(inUnits(j),'(a20,2i2)') id, tempIndivSwap
          if (j == 1) then
            write(44, '(a20)', advance = 'no') id
          end if
          if (j == nCores) then
            write(44, '(2i2)', advance='yes') tempIndivSwap
          else
            write(44, '(2i2)', advance='no') tempIndivSwap
          end if
        end do
      end do

      do i = 1, nCores
        close(inUnits(i))
      end do

      deallocate(inUnits)
      close(44)
    end if

  end subroutine CombineResults

  function ParsePedigreeAndGenotypeData(params) result(p)
    use ProgramParametersDefinition
    use PedigreeModule
    use SortingModule
    use PedigreeModule

    type(ProgramParameters), intent(in) :: params
    type(PedigreeHolder) :: p

    integer :: nAnisG
    integer :: nAnisRawPedigree
    
    call CountInData(nAnisRawPedigree, nAnisG, params)

    if (trim(params%PedigreeFile) /= "NoPedigree") then
       p = PedigreeHolder(params%GenotypeFile,nAnisRawPedigree,params%nSnp,params%GenotypeFileFormat , params%PedigreeFile )

    else
      p = PedigreeHolder(params%GenotypeFile,nAnisRawPedigree,params%nSnp, params%GenotypeFileFormat )
    endif

     if (p%nHd == 0) then
    ! TODO check if this is wanted behaviour
      p%nHd = p%nGenotyped
      p%hdMap = p%genotypeMap
      p%hdDictionary = p%genotypeDictionary
    endif
    
  end function ParsePedigreeAndGenotypeData

  function ParsePhaseData(PhaseFile, nAnisG, nSnp) result(Phase)
    use HaplotypeModule
    
    character(len=300) :: PhaseFile
    integer, intent(in) :: nAnisG, nSnp
    type(Haplotype), pointer, dimension(:,:) :: Phase

    integer :: i
    integer(kind=1), dimension (nSnp) :: ReadingVector
    character(lengan) :: dummy

    open (unit = 3, file = trim(PhaseFile), status = "old")


    allocate(Phase(nAnisG, 2))

    do i = 1, nAnisG
      read (3, *) dummy, ReadingVector(:)
      Phase(i,1) = Haplotype(ReadingVector)
      read (3, *) dummy, ReadingVector(:)
      Phase(i,2) = Haplotype(ReadingVector)
    enddo

    close(3)
  end function ParsePhaseData

  function ReadInParameterFile(filename) result (params)
    use ProgramParametersDefinition
    use AlphaHouseMod, only: parseToFirstWhitespace,splitLineIntoTwoParts,toLower
    character(*), intent(in) :: filename
    type(ProgramParameters) :: params

    double precision :: PercSurrDisagree
    integer :: TempInt, status, cl
    integer :: unit
    character (len = 300) :: OffsetVariable, hold, outputoption
    character(len=300) :: first, line
    character(len=:), allocatable::tag
    character(len=300),dimension(:),allocatable :: second

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
          endif
        endif

      case("numberofsnp")
        read(second(1), *) params%nSnp

      case("generalcoreandtaillength")
        read(second(1), *) params%params%CoreAndTailLength

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


      case("percentagegenohaplodisagree")
        read(second(1),*) params%params%PercGenoHaploDisagree

      case("genotypemissingerrorpercentage")
        read(second(1), *) params%params%GenotypeMissingErrorPercentage


      case("fulloutput")
        read(second(1),*) outputoption


      case("simulation")
        read(second(1), *) tempInt
        params%Simulation = (TempInt == 1)


      case("truephasefile")
        read(second(1), *) params%TruePhaseFile


      case("iteratemethod")
        if (allocated(second)) then
          read(second(1), *) params%params%itterateType
        endif

      case("iteratesubsetsize")
        if(allocated(second)) then
          read(second(1), *) hold
          if (hold(1:1) == "*") then
            read(hold,"(X,I2)") cl
            call get_command_argument(cl,hold)
            read(hold,*) params%params%itterateNumber
          else
            read(hold,*) params%params%itterateNumber
          end if
        end if

      case("iterateiterations")
        if (allocated(second)) then
          read(second(1), *) hold
          if (hold(1:1) == "*") then
            read(hold,"(X,I2)") cl
            call get_command_argument(cl,hold)
            read(hold,*) params%params%numIter
          else
            read(hold,*) params%params%numIter
          end if
        end if
        if (params%params%itterateType .eq. "Off") then
          params%params%numIter = 1
        end if


      case("cores")
        if (size(second) == 2) then
          read(second(1), *) params%params%startCoreChar
          read(second(2), *) params%params%endCoreChar
        endif
	
      case("minhapfreq")
        if(allocated(second)) then
          read (second, *) hold
          if (hold(1:1) == "*") then
            read(hold,"(X,I2)") cl
            call get_command_argument(cl,hold)
            read(hold,*) params%params%minHapFreq
          else
            read(hold,*) params%params%minHapFreq
          end if
        end if

      case("library")
        if (allocated(second)) then
          read(second(1), *) params%library
       endif
       
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
  if (params%params%CoreAndTailLength < params%params%Jump) then
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

  if (outputoption .eq. "Impute") then
    params%outputParams%outputFinalPhase = .false.
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
    params%outputParams%outputMistakes = .false.
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
    params%outputParams%outputMistakes = .true.
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
    params%outputParams%outputMistakes = .false.
  end if

  PercSurrDisagree = PercSurrDisagree/100
  params%params%NumSurrDisagree = int(params%params%UseSurrsN * PercSurrDisagree)
  params%params%PercGenoHaploDisagree = params%params%PercGenoHaploDisagree/100
  params%params%GenotypeMissingErrorPercentage = params%params%GenotypeMissingErrorPercentage/100

end function ReadInParameterFile


  subroutine HapCommonality(library, OutputPoint, params)
    use HaplotypeLibraryDefinition
    use ProgramParametersDefinition
    
    type(HaplotypeLibrary), intent(in) :: library
    integer, intent(in) :: OutputPoint
    class(ProgramParameters), intent(in) :: params
    
    integer :: i, SizeCore, nHaps
    character(len = 300) :: filout

    integer, allocatable, dimension (:,:) :: HapRel

    SizeCore = library%getNumSnps()
    nHaps = library%getSize()

    if (params%outputParams%outputHapCommonality) then
      write (filout, '(".",a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"Extras",a1,"HapCommonality",i0,".txt")') DASH, DASH, DASH, DASH, OutputPoint
      open (unit = 27, FILE = filout, status = 'unknown')
      
      HapRel = library%getHapRel()

      do i = 1, nHaps
	write (27, '(i10,20000F5.2,20000F5.2,20000F5.2,20000F5.2)') i, float(HapRel(i,:))/SizeCore
      enddo
    endif

    close(27)

  end subroutine HapCommonality

  subroutine WriteSurrogates(definition, OutputPoint, p, params)
    use SurrogateDefinition
    use PedigreeModule
    use ProgramParametersDefinition

    character(len = 300) :: filout
    integer :: i, j, nSurrogates

    type(Surrogate), intent(in) :: definition
    integer, intent(in) :: OutputPoint
    type(PedigreeHolder), intent(in) :: p
    type(ProgramParameters), intent(in) :: params

    integer :: nAnisG
    character(len=100) :: fmt

    nAnisG = size(definition%numOppose,1)

    if (params%outputParams%outputSurrogates) then
      write (filout, '(".",a1,"Miscellaneous",a1,"Surrogates",i0,".txt")') DASH, DASH, OutputPoint
      open (unit = 13, FILE = filout, status = 'unknown')
      write(fmt, '(a,i10,a)') '(a20,', size(definition%partition,2), 'i6)'
      do i = 1, nAnisG
	write (13, fmt) p%pedigree(p%hdMap(i))%originalId, definition%partition(i,:)
      end do
      close(13)
    end if
    
    
    if (params%outputParams%outputSurrogatesSummary) then
      write (filout, '(".",a1,"Miscellaneous",a1,"SurrogatesSummary",i0,".txt")') DASH, DASH, OutputPoint
      open (unit = 19, FILE = filout, status = 'unknown')
      do i = 1, nAnisG
	nSurrogates = 0
	do j = 1, nAnisG
	  if ((definition%numOppose(i, j) <= definition%threshold) .and. &
	   definition%enoughIncommon(i,j)) then
	    nSurrogates = nSurrogates + 1
	  end if
	enddo
	write (19, '(a20,5i8)') &
	p%pedigree(p%hdMap(i))%originalId, count(definition%partition(i,:) == 1), count(definition%partition(i,:) == 2)&
	, count(definition%partition(i,:) == 3), nSurrogates, definition%method(i)
      enddo
      close(19)
    end if

  end subroutine WriteSurrogates

  subroutine CountInData(nAnisRawPedigree, nAnisG, params)
    use ProgramParametersDefinition
    
    type(ProgramParameters), intent(in) :: params
    integer, intent(out) :: nAnisRawPedigree, nAnisG

    integer :: k
    character (len = 300) :: dumC

    nAnisRawPedigree = 0
    if (trim(params%PedigreeFile) /= "NoPedigree") then
      open (unit = 2, file = trim(params%PedigreeFile), status = "old")
      do
	read (2, *, iostat = k) dumC
	nAnisRawPedigree = nAnisRawPedigree + 1
	if (k /= 0) then
	  nAnisRawPedigree = nAnisRawPedigree - 1
	  exit
	endif
      enddo
      close(2)
      print*, " ", nAnisRawPedigree, " individuals in the pedigree file"
    endif

    nAnisG = 0
    open (unit = 3, file = trim(params%GenotypeFile), status = "old")
    do
      read (3, *, iostat = k) dumC
      nAnisG = nAnisG + 1
      if (k /= 0) then
	nAnisG = nAnisG - 1
	exit
      endif
    enddo
    if (params%GenotypeFileFormat == 2) then
      nAnisG = nAnisG /2
    end if
    close(3)
    print*, " ", nAnisG, " individuals in the genotype file"

    if (trim(params%PedigreeFile) == "NoPedigree") nAnisRawPedigree = nAnisG

  end subroutine CountInData

  subroutine MakeDirectories(params)
    use ProgramParametersDefinition
    
    type(ProgramParameters), intent(in) :: params
    
    print*, ""
    call system(RMDIR // "Miscellaneous")
    call system(RMDIR // "PhasingResults")

    call system(MD // "PhasingResults")
    call system(MD // "PhasingResults"//DASH//"HaplotypeLibrary")
    if (params%outputParams%outputHapCommonality) call system(MD // "PhasingResults"//DASH//"HaplotypeLibrary"//DASH//"Extras")
    call system(MD // "Miscellaneous")
    
    if (params%Simulation) then
      call system(RMDIR // "Simulation")
      if (params%outputParams%outputIndivMistakes .or. params%outputParams%outputIndivMistakesPercent .or. params%outputParams%outputCoreMistakesPercent .or. &
	    params%outputParams%outputMistakes) then
	call system(MD // "Simulation")
      endif
    end if

  end subroutine MakeDirectories  
  

  
  subroutine writeTimer(hours, minutes, seconds, params)
    use ProgramParametersDefinition
    integer, intent(in) :: hours, minutes
    real, intent(in) :: seconds
    type(ProgramParameters), intent(in) :: params
         
    if (params%outputParams%outputTimer) then
      open (unit = 32, file = "."//DASH//"PhasingResults"//DASH//"Timer.txt", status = "unknown")
      write(32, '(A27,A7,I3,A9,I3,A9,F6.2)') "Time Elapsed", "Hours", INT(Hours), "Minutes", INT(Minutes), "Seconds", Seconds
      close(32)
    end if
  end subroutine writeTimer
  
  subroutine WriteTestResults(results, c, p, OutputPoint, params)
    use SurrogateDefinition
    use PedigreeModule
    use CoreDefinition
    use TestResultDefinition
    use ProgramParametersDefinition

    type(TestResults), intent(in) :: results
    type(Core), intent(in) :: c
    type(PedigreeHolder), intent(in) :: p
    type(ProgramParameters), intent(in) :: params
    integer, intent(in) :: OutputPoint

    integer :: i
    integer :: nAnisG, nSnp
    character(len = 300) :: filout

    nAnisG = c % getNAnisG()
    nSNp = c % getNCoreSnp()

    if (params%outputParams%outputIndivMistakes) then
      write (filout, '(".",a1,"Simulation",a1,"IndivMistakes",i0,".txt")') DASH, DASH, OutputPoint
      open (unit = 17, FILE = filout, status = 'unknown')
      do i = 1, nAnisG
	write(17, '(a20,6i6,a6,6i6,a6,6i6,a6,6i6)') p%pedigree(p%hdMap(i))%originalId, &
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
      close(17)
    end if
    
    if (params%outputParams%outputIndivMistakesPercent) then
      write (filout, '(".",a1,"Simulation",a1,"IndivMistakesPercent",i0,".txt")') DASH, DASH, OutputPoint
      open (unit = 20, FILE = filout, status = 'unknown')

      do i = 1, nAnisG
	write (20, '(a20,6f7.1,a6,6f7.1,a6,6f7.1,a6,6f7.1)') p%pedigree(p%hdMap(i))%originalId, &
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
      close(20)
    end if
    
    if (params%outputParams%outputCoreMistakesPercent) then
      write (filout, '(".",a1,"Simulation",a1,"CoreMistakesPercent.txt")') DASH, DASH
      open (unit = 31, FILE = filout, status = 'unknown', position = 'append')
      write (31, '(6f9.4)') &
      (results%percentAll(1,ALL_,CORRECT_) + results%percentAll(2,ALL_,CORRECT_)) / 2, &
      (results%percentAll(1,HET_,CORRECT_) + results%percentAll(2,HET_,CORRECT_)) / 2, &
      (results%percentAll(1,ALL_,NOTPHASED_) + results%percentAll(2,ALL_,NOTPHASED_)) / 2, &
      (results%percentAll(1,HET_,NOTPHASED_) + results%percentAll(2,HET_,NOTPHASED_)) / 2, &
      (results%percentAll(1,ALL_,INCORRECT_) + results%percentAll(2,ALL_,INCORRECT_)) / 2, &
      (results%percentAll(1,HET_,INCORRECT_) + results%percentAll(2,HET_,INCORRECT_)) / 2
      close(31)
    end if
  end subroutine WriteTestResults

  subroutine CombineTestResults(nCores, params)
    use ProgramParametersDefinition

    integer, intent(in) :: nCores
    type(ProgramParameters) :: params

    character(len = 300) :: filout
    double precision, allocatable, dimension(:,:) :: AverageMatrix
    
    integer :: i

    if (params%outputParams%outputCoreMistakesPercent) then
      allocate(AverageMatrix(nCores, 6))
      write (filout, '(".",a1,"Simulation",a1,"CoreMistakesPercent.txt")') DASH, DASH
      open (unit = 31, FILE = filout, status = 'unknown')
      do i = 1, nCores
	read (31, *) AverageMatrix(i,:)
      end do
      write (31, *) " "
      write (31, '(6f9.4)') sum(AverageMatrix(:, 1))/nCores, sum(AverageMatrix(:, 2))/nCores, sum(AverageMatrix(:, 3))/nCores, &
      sum(AverageMatrix(:, 4))/nCores, sum(AverageMatrix(:, 5))/nCores, sum(AverageMatrix(:, 6))/nCores
      deallocate(AverageMatrix)
      close(31)
    end if
  end subroutine CombineTestResults
  
  function getHaplotypeLibraries(directory) result (libraries)
    use HaplotypeLibraryDefinition
    
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
    print*, "                              **********************                         "
    print*, "                              *                    *                         "
    print*, "                              *   AlphaPhase 1.1   *                         "
    print*, "                              *                    *                         "
    print*, "                              **********************                         "
    print*, "                                                                              "
    print*, "                    Software For Phasing and Imputing Genotypes               "
  end subroutine Header

   subroutine Titles

    call Header
    print*, ""
    print*, ""
    print*, ""
    print*, ""

  end subroutine Titles

  subroutine PrintTimerTitles(params)
    use ProgramParametersDefinition

    implicit none

    class(ProgramParameters), intent(in) :: params

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
    use ProgramParametersDefinition
    use CoreDefinition
    use HaplotypeModule
    use HaplotypeLibraryDefinition
    
    class(HaplotypeLibrary), intent(in) :: library
    integer, intent(in) :: currentcore
    type(ProgramParameters), intent(in) :: params

    integer :: i, SizeCore, nHaps
    character(len = 300) :: filout
    type(Haplotype) :: hap

    SizeCore = library%getNumSnps()

    nHaps = library%getSize()

    if (params%outputParams%outputHaplotypeLibraryText) then
      write (filout, '(".",a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"HapLib",i0,".txt")') DASH, DASH, DASH, currentcore
      open (unit = 24, FILE = filout, status = 'unknown')
    endif
    if (params%outputParams%outputHaplotypeLibraryBinary) then
      write (filout, '(".",a1,"PhasingResults",a1,"HaplotypeLibrary",a1,"HapLib",i0,".bin")') DASH, DASH, DASH, currentcore
      open (unit = 34, FILE = filout, form = "unformatted", status = 'unknown')
      
      write (34) nHaps, SizeCore
    end if
    
    do i = 1, nHaps
      hap = library%getHap(i)
      if (params%outputParams%outputHaplotypeLibraryText) then	
	write (24, '(2i6,a2,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)') &
	i, library%getHapFreq(i), " ", hap%toIntegerArray()
      end if
      if (params%outputParams%outputHaplotypeLibraryBinary) then
	write (34) hap%toIntegerArray()
      end if
    end do

    if (params%outputParams%outputHaplotypeLibraryText) then
      close(24)
    end if
    if (params%outputParams%outputHaplotypeLibraryBinary) then
      close(34)
    end if
  end subroutine WriteHapLib
  

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

  subroutine PrintVersion
    call Header
    print *
    print *, "                              Commit:   "//TOSTRING(COMMIT)
    print *, "                              Compiled: "//__DATE__//", "//__TIME__  
  end subroutine PrintVersion
end module InputOutput