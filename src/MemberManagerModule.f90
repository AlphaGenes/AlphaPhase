module MemberManagerModule
    use CoreModule

    implicit none

    type :: MemberManager

        type(CoreType), pointer :: c
        integer :: number
        integer :: curPos
        integer, dimension(:), allocatable :: order
        logical :: noneLeft
        logical :: random
    contains
        procedure :: hasNext
        procedure :: getNext
        final :: destroy      
    end type MemberManager

    interface MemberManager
        module procedure newMemberManager
    end interface MemberManager

contains

    function newMemberManager(c, iterateType, iterateNumber) result(manager)
        class(CoreType), intent(in) :: c
        type(MemberManager) :: manager
        character (len = 300) :: iterateType
        integer, intent(in) :: iterateNumber

        if (iterateType .eq. "Off") then
            call createAll(manager, c)
        end if
        if (iterateType .eq. "RandomOrder") then
            call createRandomOrder(manager, c, iterateNumber)
        end if
        if (iterateType .eq. "Cluster") then
            call createCluster(manager, c, iterateNumber)
        end if
    end function newMemberManager

    subroutine destroy(manager)
        type(MemberManager) :: manager

        if (allocated(manager%order)) then
            deallocate(manager%order)
        end if
    end subroutine destroy

    subroutine createPedigreeOrder(manager, c, number)
        class(MemberManager) :: manager
        class(CoreType), intent(in), target :: c
        integer, intent(in) :: number

        integer :: nAnisG, i

        manager%c => c

        nAnisG = c%getNAnisG()

        allocate(manager%order(nAnisG))
        do i = 1, nAnisG
            manager%order(i) = i
        end do

        manager%noneleft = .false.
        manager%number = number
        manager%curPos = 1
        manager%random = .false.
    end subroutine createInputOrder

    subroutine createRandomOrder(manager, c, number)
        use Random

        class(MemberManager) :: manager
        class(CoreType), intent(in), target :: c
        integer, intent(in) :: number

        integer :: nAnisG, secs, nCount

        manager%c => c

        call system_clock(nCount)
        secs = mod(nCount, int(1e6))

        nAnisG = c%getNAnisG()
        allocate(manager%order(nAnisG))
        call RandomOrder(manager%order, nAnisG, 1, -abs(secs))

        manager%noneleft = .false.
        manager%number = number
        manager%curPos = 1
        manager%random = .true.
    end subroutine createRandomOrder

    subroutine createCluster(manager, c, number)
        use GenotypeModule
        class(MemberManager) :: manager
        class(CoreType), intent(in), target :: c
        integer, intent(in) :: number

        logical, dimension(c%getNAnisG()) :: used
        integer :: nAnisG, numUsed, curMax, curOrder, seed, i, curIndiv, curSize
        type(Genotype), pointer :: g1, g2

        manager%c => c

        nAnisG = c%getNAnisG()
        allocate(manager%order(nAnisG))
        used = .false.
        numUsed = 0
        curOrder = 0

        do while (numUsed < nAnisG)    
            curSize = 0      
            do i = 1, nAnisG
                if (.not. used(i)) then
                    seed = i
                    used(i) = .true.
                    curOrder = curOrder + 1
                    manager%order(curOrder) = i
                    numUsed = numUsed + 1
                    curSize = 1
                    exit
                end if
            end do

            curMax = 0
            curIndiv = 1
            do while ((curSize < number) .and. (numUsed < nAnisG))
                if (.not. used(curIndiv)) then
                    g1 = c%coreAndTailGenos(seed)
                    g2 = c%coreAndTailGenos(curIndiv)
                    if (dist(g1%toIntegerArray(),g2%toIntegerArray()) <= curMax) then
                        used(curIndiv) = .true.
                        curOrder = curOrder + 1
                        manager%order(curOrder) = curIndiv
                        numUsed = numUsed + 1
                        curSize = curSize + 1
                    end if
                end if
                curIndiv = curIndiv + 1
                if (curIndiv > nAnisG) then
                    curIndiv = 1
                    curMax = curMax + 1
                end if
            end do
        end do

        manager%noneleft = .false.
        manager%number = number
        manager%curPos = 1
        manager%random = .false.

    end subroutine createCluster

    function dist(input1, input2) result (d)
        integer(kind=1), dimension(:), intent(in) :: input1, input2
        integer :: d

        integer :: i

        d = 0
        do i = 1, size(input1)
            if (input1(i) == 0 .and. input2(i) == 2) then
                d = d + 1
            end if
            if (input1(i) == 2 .and. input1(i) == 0) then
                d = d + 1
            end if
        end do
    end function dist

    subroutine createAll(manager, c)
        class(MemberManager) :: manager
        class(CoreType), intent(in), target :: c

        manager%noneleft = .false.
        call createPedigreeOrder(manager, c, c%getNAnisG())
    end subroutine createAll

    function hasNext(manager) result (has)
        class(MemberManager) :: manager
        logical :: has

        has = .not. manager%noneLeft
    end function hasNext

    function getNext(manager) result(members)
        use Random
        class(MemberManager) :: manager
        integer, dimension(:), allocatable :: members

        integer :: maxNumber, unPhased
        integer :: c, i

        integer :: nCount, secs

        unPhased = 0
        do i = 1, manager%c%getNAnisG()
            if (.not. (manager%c%phase(i,1)%fullyPhased() .and. manager%c%phase(i,2)%fullyPhased())) then
                unPhased = unPhased + 1
            end if
        end do
        if (unPhased < manager%number) then
            manager%noneLeft = .true.
            maxNumber = unPhased
        else
            maxNumber = manager%number
        end if
        allocate(members(maxNumber))
        c = 0
        do while (c < maxNumber)
            if (.not. (manager%c%phase(manager%curPos,1)%fullyPhased() .and. manager%c%phase(manager%curPos,2)%fullyPhased())) then
                c  = c + 1
                members(c) = manager%order(manager%curPos)
            end if
            if (manager%curPos == size(manager%order,1)) then
                manager%curPos = 1
                if (manager%random) then
                    call system_clock(nCount)
                    secs = mod(nCount, int(1e6))
                    call RandomOrder(manager%order, size(manager%order,1), 1, -abs(secs))
                end if
                manager%noneLeft = .true.
            else
                manager%curPos = manager%curPos + 1
            end if
        end do

    end function getNext

end module MemberManagerModule