module NRMCode
  implicit none
  
  integer, parameter :: NRMmemTemp = 20000
  integer, parameter :: shellmax = 50000
  
  !! Only marginally better than global but as the NRM code is going to be ripped out soon it doesn't seem worth fixing.
  real(kind = 4), allocatable :: xnumrelmatHold(:)
  integer, allocatable :: RecSire(:), RecDam(:)
  integer :: NRMmem, shell
  character(20), allocatable :: nrmped(:,:)
  character(20), allocatable :: Id(:)
  integer, allocatable :: seqsire(:), seqdam(:), RecodeGenotypeId(:)
contains

  function createNRM(p) result(PseudoNRM)
    use Parameters, only: FullFileOutput, NrmThresh
    use Constants
    use PedigreeDefinition
    
    type(Pedigree), intent(in) :: p
    
    integer(kind = 1), allocatable, dimension (:,:) :: PseudoNRM
    integer :: i, j, nAnisP
    double precision :: value, valueS, valueD
    integer :: nAnisG
    
    nAnisG = p%getNAnis()
    allocate(RecodeGenotypeId(nAnisG))
    call PedigreeViewerRecode(size(nrmped,1), nAnisP, nrmped, id)
    
    do i = 1, nAnisG
      do j = 1, nAnisP
	if (Id(j) == p%getID(i)) then
	!if (ped(j,1) == GenotypeId(i)) then
	  RecodeGenotypeId(i) = j
	  exit
	end if
      end do
    end do
    
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
	open (unit = 9, file = ".\Miscellaneous\GenotypedPseudoNRM.txt", status = "unknown")
      else
	open (unit = 9, file = "./Miscellaneous/GenotypedPseudoNRM.txt", status = "unknown")
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
	write (9, *) p%getID(i), PseudoNRM(i,:)
      endif
    enddo
    
    close(9)
    
    return

    !!!!! END NRM
  end function createNRM
  
  recursive function xnumrelmat(i, j) result (xA)
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
  
  subroutine PedigreeViewerRecode(nobs, nAnisPedigree, ped, id)
    implicit none

    integer :: nobs, nAnisPedigree
    character(20), dimension(:,:), intent(in) :: ped
    character(20), allocatable :: Id(:)

    integer, allocatable :: passedorder(:)
    character (20), allocatable :: holdid(:), holdsireid(:), holddamid(:)
    character (20), allocatable :: Sortedid(:), Sortedsire(:), Sorteddam(:)
    character(20), allocatable :: sire(:), dam(:)
    character (20) :: IDhold
    integer, allocatable :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
    integer, allocatable :: OldN(:), NewN(:), holdsire(:), holddam(:)
    integer :: mode ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
    integer :: i, j, k, kk, newid, itth, itho, ihun, iten, iunit
    integer :: nsires, ndams, newsires, newdams, nbisexuals, flag
    integer :: ns, nd, iextra, oldnobs, kn, kb, oldkn, ks, kd
    integer :: Noffset, Limit, Switch, ihold, ipoint

    mode = 1
  !  allocate(id(0:nobs), sire(nobs), dam(nobs), seqid(nobs), seqsire(nobs), seqdam(nobs))
    allocate(id(0:nobs), sire(nobs), dam(nobs), seqsire(nobs), seqdam(nobs))

    do i = 1, nobs
      id(i) = ped(i, 1)
      sire(i) = ped(i, 2)
      dam(i) = ped(i, 3)
    end do

    !!!!!!!!!!!!!!!!!!
    ! INSERT DUMMY IDS FOR UNKNWON PARENTS
    !!!!!!!!!!!!!!!!!!

    do j = 1, nobs
      if (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') then
	dam(j) = '0'
	seqdam(j) = 0
      endif
      if (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') then
	sire(j) = '0'
	seqsire(j) = 0
      endif
    enddo
    if (mode .eq. 1) then
      !print*,  ' Inserting dummy IDs into Pedigree ... '
      newid = 0
      do j = 1, nobs
	if (((sire(j) == '0').and.(dam(j) .ne. '0')) .or. ((sire(j) .ne. '0').and.(dam(j) == '0'))) then
	  newid = newid + 1
	  if (newid .gt. 99999) then
	    PRINT*, newid, ' ...'
	    stop 'too many dummy single parent IDs'
	  endif
	  itth = int(newid/10000)
	  itho = int(newid/1000) - 10 * itth
	  ihun = int(newid/100) - 10 * itho - 100 * itth
	  iten = int(newid/10) - 10 * ihun - 100 * itho - 1000 * itth
	  iunit = newid - 10 * iten - 100 * ihun - 1000 * itho - 10000 * itth
	  if (sire(j) == '0') sire(j) = 'dum' // achar(48 + itth) // achar(48 + itho) // achar(48 + ihun) // achar(48 + iten) // achar(48 + iunit)
	  if (dam(j) == '0') dam(j) = 'dum' // achar(48 + itth) // achar(48 + itho) // achar(48 + ihun) // achar(48 + iten) // achar(48 + iunit)
	endif
      enddo
    endif  

    !!!!!!!!!!!!!!!!!!
    ! SORT SIRES BY ID AND CREATE LIST OF UNIQUE SIRES (sortedsire)
    !!!!!!!!!!!!!!!!!!

    !print*,  ' Sorting Sires ... '
    allocate(Sortedid(nobs), SortedIdIndex(nobs))
    Sortedid(1:nobs) = sire(1:nobs)
    Noffset = int(nobs/2)
    do while (Noffset > 0)
      Limit = nobs - Noffset
      switch = 1
      do while (Switch .ne. 0)
	Switch = 0
	do i = 1, Limit
	  if (Sortedid(i) .gt. Sortedid(i + Noffset)) then
	    IDhold = Sortedid(i)
	    Sortedid(i) = Sortedid(i + Noffset)
	    Sortedid(i + Noffset) = IDhold
	    Switch = i
	  endif
	enddo
	Limit = Switch - Noffset
      enddo
      Noffset = int(Noffset/2)
    enddo
    nsires = 0
    if (Sortedid(1) /= '0') nsires = 1
    do i = 2, nobs
      if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') nsires = nsires + 1
    end do
    allocate(Sortedsire(0:nsires), SortedSireIndex(nsires))
    Sortedsire(0) = '0'
    nsires = 0
    if (Sortedid(1) /= '0') then
      nsires = 1
      Sortedsire(1) = Sortedid(1)
    endif
    do i = 2, nobs
      if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') then
	nsires = nsires + 1
	Sortedsire(nsires) = Sortedid(i)
      endif
    end do

    !!!!!!!!!!!!!!!!!!
    ! SORT DAMS BY ID AND CREATE LIST OF UNIQUE SIRES (sorteddam)
    !!!!!!!!!!!!!!!!!!

    !print*,  ' Sorting Dams ... '
    Sortedid(1:nobs) = dam(1:nobs)
    Noffset = int(nobs/2)
    do while (Noffset > 0)
      Limit = nobs - Noffset
      switch = 1
      do while (Switch .ne. 0)
	Switch = 0
	do i = 1, Limit
	  if (Sortedid(i) .gt. Sortedid(i + Noffset)) then
	    IDhold = Sortedid(i)
	    Sortedid(i) = Sortedid(i + Noffset)
	    Sortedid(i + Noffset) = IDhold
	    Switch = i
	  endif
	enddo
	Limit = Switch - Noffset
      enddo
      Noffset = int(Noffset/2)
    enddo
    nDams = 0
    if (Sortedid(1) /= '0') nDams = 1
    do i = 2, nobs
      if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') nDams = nDams + 1
    end do
    allocate(Sorteddam(0:nDams), SortedDamIndex(ndams))
    Sorteddam(0) = '0'
    nDams = 0
    if (Sortedid(1) /= '0') then
      nDams = 1
      Sorteddam(1) = Sortedid(1)
    endif
    do i = 2, nobs
      if (Sortedid(i) /= Sortedid(i - 1).and.Sortedid(i) /= '0') then
	nDams = nDams + 1
	Sorteddam(nDams) = Sortedid(i)
      endif
    end do

    !!!!!!!!!!!!!!!!!!
    ! SORT ANIMALS BY ID
    !!!!!!!!!!!!!!!!!!

    !print*,  ' Sorting IDs ... '
    Sortedid(1:nobs) = id(1:nobs)
    do i = 1, nobs
      SortedIdIndex(i) = i
    end do
    Noffset = int(nobs/2)
    do while (Noffset > 0)
      Limit = nobs - Noffset
      switch = 1
      do while (Switch .ne. 0)
	Switch = 0
	do i = 1, Limit
	  if (Sortedid(i) .gt. Sortedid(i + Noffset)) then
	    IDhold = Sortedid(i)
	    Sortedid(i) = Sortedid(i + Noffset)
	    Sortedid(i + Noffset) = IDhold
	    ihold = SortedIdIndex(i)
	    SortedIdIndex(i) = SortedIdIndex(i + Noffset)
	    SortedIdIndex(i + Noffset) = ihold
	    Switch = i
	  endif
	enddo
	Limit = Switch - Noffset
      enddo
      Noffset = int(Noffset/2)
    enddo

    !!!!!!!!!!!!!!!!!!
    ! CHECK FOR DUPLICATE IDS
    !!!!!!!!!!!!!!!!!!

    !print*,  ' Check for duplicate IDs ... '
    flag = -1
    do i = 2, nobs
      if (Sortedid(i) == Sortedid(i - 1)) then
	if (flag == -1) then
	  open (4, file = 'ID_err.txt', status = 'unknown')
	  write(4, *) 'Duplicated IDs ...'
	  flag = 0
	end if
	write(4, *) Sortedid(i)
	flag = flag + 1
      end If
    enddo
    if (flag > -1) then
      close (4)
      print*, flag, ' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
      !pause
    end if

    !!!!!!!!!!!!!!!!!!
    ! SET SORTEDSIREINDEX to position in id list or new position at end if not already there
    !!!!!!!!!!!!!!!!!!

    !print*,  ' Males ... '
    !print*,  '  Find or set sire indices ... '
    newsires = 0
    do j = 1, nsires
      !check if already listed as an individual
      ipoint = int(nobs/2)
      Noffset = int(ipoint/2)
      do while (Noffset > 1)
	if (Sortedsire(j) .lt. Sortedid(ipoint)) then
	  ipoint = ipoint - Noffset
	  Noffset = int(Noffset/2)
	else
	  ipoint = ipoint + Noffset
	  Noffset = int(Noffset/2)
	endif
      enddo
      kn = 0
      if (Sortedsire(j) == Sortedid(ipoint)) kn = 1
      do while (ipoint < nobs .and. kn == 0 .and. Sortedsire(j) > Sortedid(ipoint))
	ipoint = ipoint + 1
      enddo
      if (Sortedsire(j) == Sortedid(ipoint)) kn = 1
      do while (ipoint > 1 .and. kn == 0 .and. Sortedsire(j) < Sortedid(ipoint))
	ipoint = ipoint - 1
      enddo
      if (Sortedsire(j) == Sortedid(ipoint)) kn = 1
      if (kn == 1) then
	SortedSireIndex(j) = SortedIdIndex(ipoint)
      else ! sire is unlisted base sire
	newsires = newsires + 1
	SortedSireIndex(j) = nobs + newsires ! for now
      endif
    end do
    allocate(holdsireid(newsires))
    kn = 0
    do j = 1, nsires
      if (SortedSireIndex(j) > nobs) then
	kn = kn + 1
	holdsireid(SortedSireIndex(j) - nobs) = Sortedsire(j)
      end if
    enddo
    if (kn /= newsires) stop'newsires error'
    !print*,  '  Find seqsire ... '

    !!!!!!!!!!!!!!!!!!
    ! SET SEQSIRE to the position of the sire in the sorted sire list?
    !!!!!!!!!!!!!!!!!!

    do j = 1, nobs
      If (sire(j) == '0') Then
	seqsire(j) = 0
      else
	ipoint = int(nsires/2)
	Noffset = int(ipoint/2)
	do while (Noffset > 1)
	  if (sire(j) .lt. Sortedsire(ipoint)) then
	    ipoint = ipoint - Noffset
	    Noffset = int(Noffset/2)
	  else
	    ipoint = ipoint + Noffset
	    Noffset = int(Noffset/2)
	  endif
	enddo
	kn = 0
	if (sire(j) == Sortedsire(ipoint)) kn = 1
	do while (ipoint < nsires .and. kn == 0 .and. sire(j) > Sortedsire(ipoint))
	  ipoint = ipoint + 1
	enddo
	if (sire(j) == Sortedsire(ipoint)) kn = 1
	do while (ipoint > 1 .and. kn == 0 .and. sire(j) < Sortedsire(ipoint))
	  ipoint = ipoint - 1
	enddo
	if (sire(j) == Sortedsire(ipoint)) kn = 1
	if (kn == 1) then
	  seqsire(j) = SortedSireIndex(ipoint)
	else
	  print*, ' Error: Sire missing: ', sire(j)
	  stop
	endif
      endif
    enddo

    !!!!!!!!!!!!!!!!!!
    ! REPEAT FOR DAMS with extra bisexual check
    !!!!!!!!!!!!!!!!!!


    !print*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'
    !print*,  ' Females ... '
    !print*,  '  Find or set dam indices ... '
    newdams = 0
    nbisexuals = 0
    do j = 1, ndams
      !check if already listed as an individual
      ipoint = int(nobs/2)
      Noffset = int(ipoint/2)
      do while (Noffset > 1)
	if (Sorteddam(j) .lt. Sortedid(ipoint)) then
	  ipoint = ipoint - Noffset
	  Noffset = int(Noffset/2)
	else
	  ipoint = ipoint + Noffset
	  Noffset = int(Noffset/2)
	endif
      enddo
      kn = 0
      if (Sorteddam(j) == Sortedid(ipoint)) kn = ipoint ! store ipoint here as ipoint can change with bisexuals
      do while (ipoint < nobs .and. kn == 0 .and. Sorteddam(j) > Sortedid(ipoint))
	ipoint = ipoint + 1
      enddo
      if (Sorteddam(j) == Sortedid(ipoint)) kn = ipoint
      do while (ipoint > 1 .and. kn == 0 .and. Sorteddam(j) < Sortedid(ipoint))
	ipoint = ipoint - 1
      enddo
      if (Sorteddam(j) == Sortedid(ipoint)) kn = ipoint
      !check if already listed as a sire (and therefore bisexual)
      ipoint = int(nsires/2)
      Noffset = int(ipoint/2)
      do while (Noffset > 1)
	if (Sorteddam(j) .lt. Sortedsire(ipoint)) then
	  ipoint = ipoint - Noffset
	  Noffset = int(Noffset/2)
	else
	  ipoint = ipoint + Noffset
	  Noffset = int(Noffset/2)
	endif
      enddo
      kb = 0
      if (Sorteddam(j) == Sortedsire(ipoint)) kb = 1
      do while (ipoint < nsires .and. kb == 0 .and. Sorteddam(j) > Sortedsire(ipoint))
	ipoint = ipoint + 1
      enddo
      if (Sorteddam(j) == Sortedsire(ipoint)) kb = 1
      do while (ipoint > 1 .and. kb == 0 .and. Sorteddam(j) < Sortedsire(ipoint))
	ipoint = ipoint - 1
      enddo
      if (Sorteddam(j) == Sortedsire(ipoint)) kb = 1
      if (kb == 1) then
	nbisexuals = nbisexuals + 1
	open (5, file = 'bisex.txt', position = 'append')
	write(5, *) Sorteddam(j)
	close(5)
      endif
      if (kb == 1) then
	SorteddamIndex(j) = SortedSireIndex(ipoint)
      elseif (kn >= 1) then
	SorteddamIndex(j) = SortedIdIndex(kn)
      else ! dam is unlisted base dam
	newdams = newdams + 1
	SorteddamIndex(j) = nobs + newsires + newdams ! for now
      endif
    end do
    if (nbisexuals > 0) then
      print*, nbisexuals, ' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'
      !pause
    endif
    allocate(holddamid(newdams))
    kn = 0
    do j = 1, ndams
      if (SortedDamIndex(j) > nobs + newsires) then
	kn = kn + 1
	holddamid(SortedDamIndex(j) - nobs - newsires) = Sorteddam(j)
      end if
    enddo
    if (kn /= newdams) stop'newdams error'
    !print*,  '  Find seqdam ... '
    do j = 1, nobs
      if (dam(j) == '0') then
	seqdam(j) = 0
      else
	ipoint = int(ndams/2)
	Noffset = int(ipoint/2)
	do while (Noffset > 1)
	  if (dam(j) .lt. Sorteddam(ipoint))then
	    ipoint = ipoint - Noffset
	    Noffset = int(Noffset/2)
	  else
	    ipoint = ipoint + Noffset
	    Noffset = int(Noffset/2)
	  endif
	enddo
	kn = 0
	if (dam(j) == Sorteddam(ipoint)) kn = 1
	do while (ipoint < ndams .and. kn == 0 .and. dam(j) > Sorteddam(ipoint))
	  ipoint = ipoint + 1
	enddo
	if (dam(j) == Sorteddam(ipoint)) kn = 1
	do while (ipoint > 1 .and. kn == 0 .and. dam(j) < Sorteddam(ipoint))
	  ipoint = ipoint - 1
	enddo
	if (dam(j) == Sorteddam(ipoint)) kn = 1
	if (kn == 1) then
	seqdam(j) = SorteddamIndex(ipoint)
      else
	print*, ' Error: dam missing: ', dam(j)
	stop
	endif

      endif

    enddo
    !print*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'
    !print*,  ' Arranging unlisted base parents ... '
    iextra = newsires + newdams
    if (iextra > 0) then
      !print*, ' ', iextra, ' unlisted base parents found.'
      ! SortedId and SortedIdIndex just used as a holder while redimensioning

      ! New animals will go at start
      Sortedid(1:nobs) = id(1:nobs)
      deallocate (id)
      allocate(id(nobs + iextra))
      id(1 + iextra:nobs + iextra) = Sortedid(1:nobs)

      Sortedid(1:nobs) = sire(1:nobs)
      deallocate (sire)
      allocate(sire(nobs + iextra))
      sire(1 + iextra:nobs + iextra) = Sortedid(1:nobs)

      Sortedid(1:nobs) = dam(1:nobs)
      deallocate (dam)
      allocate(dam(nobs + iextra))
      dam(1 + iextra:nobs + iextra) = Sortedid(1:nobs)

      SortedIdIndex(1:nobs) = seqsire(1:nobs)
      deallocate (seqsire)
      allocate(seqsire(nobs + iextra))
      seqsire(1 + iextra:nobs + iextra) = SortedIdIndex(1:nobs)

      SortedIdIndex(1:nobs) = seqdam(1:nobs)
      deallocate (seqdam)
      allocate(seqdam(nobs + iextra))
      seqdam(1 + iextra:nobs + iextra) = SortedIdIndex(1:nobs)

  !    deallocate (seqid)
  !    allocate(seqid(nobs + iextra))

    endif

    !print*, ' Inserting unlisted base parents ...'

    oldnobs = nobs
    nobs = nobs + iextra
    !print*, ' Total number of animals = ',nobs

    allocate (passedorder(nobs))
    passedorder = 0

    do i = 1 + iextra, nobs
      passedorder(i) = i - iextra

      if (sire(i) == '0')then
	seqsire(i) = 0
      else
	seqsire(i) = iextra + seqsire(i)
	If (seqsire(i) > nobs) seqsire(i) = seqsire(i) - nobs ! for unlisted sires
      end If

      if (dam(i) == '0') then
	seqdam(i) = 0
      else
	seqdam(i) = iextra + seqdam(i)
	if (seqdam(i) > nobs) seqdam(i) = seqdam(i) - nobs
      end if
    enddo

    do i = 1, newsires
      id(i) = holdsireid(i)
      passedorder(i) = 0
      seqsire(i) = 0
      seqdam(i) = 0
    enddo !i
    do i = newsires + 1, newsires + newdams
      id(i) = holddamid(i - newsires)
      passedorder(i) = 0
      seqsire(i) = 0
      seqdam(i) = 0
    enddo !i
    deallocate(holdsireid, holddamid, SortedIdIndex, SortedId)
    flag = 0
    do i = 1, nobs
      if (i <= seqsire(i) .Or. i <= seqdam(i)) flag = 1
    enddo
    if (flag == 0) then
      goto 8012
    endif
    !print*, ' Re-Ordering pedigree ...'
    allocate ( OldN(0:nobs), NewN(0:nobs))
    allocate ( holdid(0:nobs), holdsire(nobs), holddam(nobs))
    OldN(0) = 0
    NewN = 0
    !seqsire(0) = 0 !not needed !
    !seqdam(0) = 0
    holdid(1:nobs) = id(1:nobs)
    holdsire = seqsire
    holddam = seqdam
    !Find base ancestors ...
    kn = 0
    do i = 1, nobs
      if (seqsire(i) == 0 .And. seqdam(i) == 0) then
	kn = kn + 1
	NewN(i) = kn
	OldN(kn) = i
      end if
    enddo !i
    !Re-order pedigree ...
    NewN(0) = nobs + 1
    flag = 0
    do while (kn < nobs)
      oldkn = kn
      do i = 1, nobs
	if (NewN(i) == 0) then !And id(i) <> 'UniqueNULL' Then
	  Ks = seqsire(i)
	  Kd = seqdam(i)
	  if (NewN(Ks) > 0 .And. NewN(Kd) > 0) then
	    kn = kn + 1
	    NewN(i) = kn
	    OldN(kn) = i
	  end if
	end if
      enddo !i
      ! to avoid hang on unexpected problem ...
      if (kn == oldkn) then
	flag = flag + 1
      else
	flag = 0
      endif

      if (flag > 10) then
	open(6, file = 'ped_err.txt', status = 'unknown')
	write(6, *) 'Pedigree errors found involving two or more of the following relationships ...'
	write(6, *)
	write(6, *) '       Index numbers are followed by names.'
	write(6, *) '       Index number 0 means unknown, whence name is blank.'
	write(6, *)
	do i = 1, nobs
	  if (NewN(i) == 0) then
	    write(6, *) 'Individual:', i, ':  ', id(i)
	    write(6, *) '    Father:', seqsire(i), ':  ', id(seqsire(i))
	    write(6, *) '    Mother:', seqdam(i), ':  ', id(seqdam(i))
	    write(6, *)
	  end if
	enddo !i
	close (6)
	print*, 'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
	stop
      endif
    enddo
    NewN(0) = 0
    do i = 1, nobs
      id(i) = holdid(OldN(i))
    enddo
    do i = 1, nobs
      seqsire(i) = NewN(holdsire(OldN(i)))
      seqdam(i) = NewN(holddam(OldN(i)))
      if (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
	print*, 'out of order'
	stop
      endif
    enddo

    deallocate ( OldN, NewN, holdid, holdsire, holddam)

    8012 do i = 1, nobs
    !  seqid(i) = i
    enddo
    nAnisPedigree = nobs

    deallocate(Sortedsire, Sorteddam, sire, dam, SortedSireIndex, SortedDamIndex)

  end subroutine PedigreeViewerRecode

end module NRMcode