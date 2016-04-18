module NRMCode
  integer, parameter :: NRMmemTemp = 20000
  integer, parameter :: shellmax = 50000
  
  !! Only marginally better than global but as the NRM code is going to be ripped out soon it doesn't seem worth fixing.
  real(kind = 4), allocatable :: xnumrelmatHold(:)
  integer, allocatable :: RecSire(:), RecDam(:)
  integer :: NRMmem, shell
contains
  
  function createNRM result(PseudoNRM)
    use Global
    use GlobalPedigree
    
    integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM
    
    allocate(PseudoNRM(nAnisG, nAnisG))
    allocate(RecSire(0:nAnisP))
    allocate(RecDam(0:nAnisP))
    RecSire(0) = 0
    RecDam(0) = 0
    do i = 1, nAnisP
      RecSire(i) = seqsire(i)
      RecDam(i) = seqdam(i)
    enddo
    
    !!!!! START NRM !!!!!
    NRMmem = NRMmemTemp
    if (NRMmem > nAnisP) NRMmem = nAnisP
!    shellmax = 50000
    allocate(xnumrelmatHold(-1 * NRMmem: NRMmem * (NRMmem + 1)/2))
    xnumrelmatHold = -9.
    xnumrelmatHold(-1 * NRMmem: 0) = 0.

    if (FullFileOutput == 1) then
!      allocate(NRM(nAnisG, nAnisG))
      if (WindowsLinux == 1) then
	!open (unit = 8, file = ".\Miscellaneous\GenotypedNRM.txt", status = "unknown")
	open (unit = 9, file = ".\Miscellaneous\GenotypedPseudoNRM.txt", status = "unknown")
	!open (unit = 10, file = ".\Miscellaneous\SummaryOfGenotypedNRM.txt", status = "unknown")
	!open (unit = 11, file = ".\Miscellaneous\AlleleFrequency.txt", status = "unknown")
	!open (unit = 12, file = ".\Miscellaneous\GenotypedMarkerNRM.txt", status = "unknown")
      else
	!open (unit = 8, file = "./Miscellaneous/GenotypedNRM.txt", status = "unknown")
	open (unit = 9, file = "./Miscellaneous/GenotypedPseudoNRM.txt", status = "unknown")
	!open (unit = 10, file = "./Miscellaneous/SummaryOfGenotypedNRM.txt", status = "unknown")
	!open (unit = 11, file = "./Miscellaneous/AlleleFrequency.txt", status = "unknown")
	!open (unit = 12, file = "./Miscellaneous/GenotypedMarkerNRM.txt", status = "unknown")
      endif
    end if

    print*, " "
    print*, " Making NRM"
    PseudoNRM = 0
    do i = 1, nAnisG
      if (mod(i, 400) == 0) print*, "   NRM done for genotyped individual --- ", i
      if (FullFileOutput == 1) then
!	shellwarning = 0
	do j = i, nAnisG
	  shell = 0
	  value = xnumrelmat(RecodeGenotypeId(i), RecodeGenotypeId(j))
!	  NRM(i, j) = value
!	  NRM(j, i) = value
	enddo
      endif
!      shellwarning = 0
      do j = i, nAnisG
	shell = 0
	valueS = xnumrelmat(seqsire(RecodeGenotypeId(i)), RecodeGenotypeId(j))
	shell = 0
	valueD = xnumrelmat(seqdam(RecodeGenotypeId(i)), RecodeGenotypeId(j))
	if ((valueS > NrmThresh).and.(valueD <= NrmThresh)) PseudoNRM(i, j) = 1
	if ((valueS <= NrmThresh).and.(valueD > NrmThresh)) PseudoNRM(i, j) = 2
	PseudoNRM(j, i) = PseudoNRM(i, j)
      enddo
      if (FullFileOutput == 1) then
!	if (nAnisG < 20000) then
!	  write (8, "(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)") GenotypeId(i), NRM(i,:)
!	else
!	  write (8, *) GenotypeId(i), NRM(i,:)
!	end if
	write (9, *) GenotypeId(i), PseudoNRM(i,:)
      endif
    enddo
    
    close(9)
    close(10)
    
    return

    !!!!! END NRM
  end function createNRM
  
  recursive function xnumrelmat(i, j) result (xA)
    use GlobalPedigree
    use Global
    implicit none

    integer :: i, j
    real :: xA

    shell = shell + 1
    if (i .eq. 0 .or. j .eq. 0 .or. shell > shellmax) then
      xA = 0.0
  !    if (shell > shellmax) shellWarning = shellWarning + 1
      return
    endif

    If (i .le. NRMmem .And. j .le. NRMmem) Then
      xA = xnumrelmat_mem(i, j)
    else
      IF (i .eq. j)then
      xA = 1 + .5 * xnumrelmat(RecSire(i), RecDam(i))
    elseIF (i .lt. j)then
      xA = .5 * (xnumrelmat(i, RecSire(j)) + xnumrelmat(i, RecDam(j)))
    elseIF (j .lt. i)then
      xA = .5 * (xnumrelmat(j, RecSire(i)) + xnumrelmat(j, RecDam(i)))
      endif
    endif

  end function xnumrelmat

  !########################################################################################################################################################################################################

  RECURSIVE function xnumrelmat_mem(i, j) RESULT (xA)
    use GlobalPedigree
    use Global
    implicit none


    INTEGER :: i, j, k
    REAL :: xA

    ! Addressing vector xnumrelmatHold:  k=(i-1)*(float(NRMmem)-float(i)/2)+j
    shell = shell + 1
    !IF (shell>4000 .AND. MOD(shell,1000).eq.0 )   PRINT*, shell
    if (i .eq. 0 .or. j .eq. 0 .or. shell > shellmax) then
      xA = 0.0
  !    if (shell > shellmax) shellWarning = shellWarning + 1
      return
    endif

    if (i .le. j) then
      k = (i - 1) * (float(NRMmem) - float(i) / 2) + j
    else
      k = (j - 1) * (float(NRMmem) - float(j) / 2) + i
    endif

    if (xnumrelmatHold(k) > -1) then
      xA = xnumrelmatHold(k)
      return
    endif


    IF (i .eq. j)then
      xA = 1 + .5 * xnumrelmat_mem(RecSire(i), RecDam(i))
    elseIF (i .lt. j)then
      xA = .5 * (xnumrelmat_mem(i, RecSire(j)) + xnumrelmat_mem(i, RecDam(j)))
    elseIF (j .lt. i)then
      xA = .5 * (xnumrelmat_mem(j, RecSire(i)) + xnumrelmat_mem(j, RecDam(i)))
    endif

    xnumrelmatHold(k) = xA
  end function xnumrelmat_mem
end module NRMcode