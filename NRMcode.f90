module NRMCode
  integer, parameter :: NRMmemTemp = 20000
contains
  
  function createNRM result(PseudoNRM)
    use Global
    use GlobalPedigree
    
    integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM
    
    allocate(PseudoNRM(nAnisG, nAnisG))
    
    !!!!! START NRM !!!!!
    NRMmem = NRMmemTemp
    if (NRMmem > nAnisP) NRMmem = nAnisP
    shellmax = 50000
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
	shellwarning = 0
	do j = i, nAnisG
	  shell = 0
	  value = xnumrelmat(RecodeGenotypeId(i), RecodeGenotypeId(j))
!	  NRM(i, j) = value
!	  NRM(j, i) = value
	enddo
      endif
      shellwarning = 0
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
end module NRMcode