module HaplotypeLibraryPhasing
  use ConstantModule
  implicit none

  integer, parameter, private :: nMaxRounds = 100
    
contains
  subroutine UpdateHapLib(c, library, percminpresent, minoverlap, percGenoHaploDisagree)
    use HaplotypeLibraryDefinition
    use HaplotypeModule
    use CoreDefinition
    use Random
    
    type(Core), intent(in) :: c
    type(HaplotypeLibrary), intent(in) :: library
    integer, intent(in) :: minoverlap
    double precision, intent(in) :: PercGenoHaploDisagree, percminpresent
    
    type(Haplotype) :: hap

    integer :: i, errorallow, minpresent
    
    errorallow = int(percGenoHaploDisagree * c%getNCoreSnp())
    minpresent = max(int(percminpresent * c%getNCoreSnp()),minoverlap)
    
    do i = 1, c % getNAnisG()
      if (mod(i,2000) == 0) then
	print *, "Done", i
      end if
      !Paternal Haps
      hap = c % phase(i, 1)
      if (hap%numberNotMissing() >= minpresent) then
	call newHaplotype(c, i, 1, library, minoverlap, errorallow)
      endif
      
      !Maternal Haps
      hap = c % phase(i, 2)
      if (hap%numberNotMissing() >= minpresent) then
	call newHaplotype(c, i, 2, library, minoverlap, errorallow)
      endif
   enddo
  end subroutine UpdateHapLib

  subroutine newHaplotype(c, animal, phase, library, minoverlap, errorallow)
    use HaplotypeModule
    use CoreDefinition
    use HaplotypeLibraryDefinition

    class(Core) :: c
    integer, intent(in) :: animal, phase, minoverlap, errorallow
    class(HaplotypeLibrary) :: library

    integer :: id
    integer, dimension(:), allocatable :: ids
    integer, dimension(:), allocatable :: newkeys, oldkeys
    type(Haplotype) :: hap, merged
    
    integer :: i, j, k, n
    
    hap = c % phase(animal, phase)
      ids = library%matchWithErrorAndMinOverlap(hap,0,minoverlap)
      if (size(ids) == 0) then
	id = library%addHap(hap)
      end if
      if (size(ids) == 1) then
	id = ids(1)
	if (hap%equalhap(library%newstore(id))) then
	  library%hapfreq(id) = library%hapfreq(id) + 1
	else
	  oldkeys = getKeys(library%newstore(id),library%key)
	  merged = hap%mergeMod(library%newstore(id))
	  do i = 1, c%getNAnisG()
	    do j = 1, 2
	      if (c%hapanis(i,j) == id) then
		c%phase(i,j) = merged
	      end if
	    end do
	  end do
	  
	  library%newstore(id) = merged	  
	  
	  do i = 1, size(oldkeys)
	    call library%partialMap(oldkeys(i))%list_remove(id)
	  end do
	  newkeys = getKeys(merged,library%key)
	  do i = 1, size(newkeys)
	    call library%partialMap(newkeys(i))%list_add(id)
	  end do
	  
	  library%hapfreq(id) = library%hapfreq(id) + 1
	  c%phase(animal,phase) = merged
	end if
      end if
      if (size(ids) > 1) then
	id = ids(1) 
	merged = hap
	do i = 1, size(ids)
	  merged = merged%mergeMod(library%newstore(ids(i)))
	end do
	
	if (merged%numberError() <= ErrorAllow) then
	  oldkeys = getKeys(library%newstore(id),library%key)
	  library%newstore(id) = merged
	  do n = 1, size(ids)
	    do i = 1, c%getNAnisG()
	      do j = 1, 2
		if (c%hapanis(i,j) == ids(n)) then
		  c%phase(i,j) = merged
		  c%hapanis(i,j) = id
		end if
	      end do
	    end do
	  end do
	  c%phase(animal,phase) = merged
	  
	  do i = 1, size(oldkeys)
	    call library%partialMap(oldkeys(i))%list_remove(id)
	  end do
	  newkeys = getKeys(merged,library%key)
	  do i = 1, size(newkeys)
	    call library%partialMap(newkeys(i))%list_add(id)
	  end do

	  do i = 2, size(ids)
	    library%hapfreq(id) = library%hapfreq(id) + library%hapfreq(ids(i))
	    call library%removehap(ids(i))
	    do j = 1, c%getNAnisG()
	      do k = 1, 2
		if (c%hapanis(j,k) > ids(i)) then
		  c%hapanis(j,k) = c%hapanis(j,k) - 1
		end if
	      end do
	    end do
	    do j = i+1, size(ids)
		if (ids(j) > ids(i)) then
		  ids(j) = ids(j) - 1
		end if
	    end do
	  end do	  


	else
	  id = library%addHap(hap)
	end if
      end if

       c%hapAnis(animal, phase) = id
  end subroutine newHaplotype
  
  subroutine imputeFromLib(library, c, PercGenoHaploDisagree, percminpresent, minoverlap, minHapFreq, percMinToKeep, quiet)
    use HaplotypeLibraryDefinition
    use CoreDefinition
    
    type(HaplotypeLibrary), intent(in) :: library
    class(Core) :: c
    double precision, intent(in) :: PercGenoHaploDisagree
    integer, intent(in) :: minHapFreq, minoverlap
    double precision :: percMinToKeep, percminpresent
    logical, intent(in) :: quiet
    
    integer :: nOldHaps, iterations, errorallow, minpresent
    
    errorallow = int(percGenoHaploDisagree * c%getNCoreSnp())
    minpresent = max(int(percminpresent * c%getNCoreSnp()), minoverlap)
    
    call complementStart(library, c, errorallow, minpresent, minoverlap, minhapfreq)
    if (.not. quiet) then
      print '(6x, a21, 5x, f6.2, a8, i7, a11)', " After complementing ", c%getTotalYield(), "% Yield ", &
	library%numberPercentPhased(percMinToKeep), " Haplotypes"
    end if
    
    nOldHaps = 0
    iterations = 0
    do while (nOldHaps /= library%getSize())
      nOldHaps = library%getSize()
      call singleImputationRound(library, c, errorallow, minpresent, minoverlap,  minhapfreq)
      iterations = iterations + 1
      if (.not. quiet) then
	print '(6x, a7, i2, a12, 5x, f6.2, a8, i7, a11)', " After ", iterations, " iterations ", c%getTotalYield(), "% Yield ", &
	  library%numberPercentPhased(percMinToKeep), " Haplotypes"
      end if
    end do
  end subroutine imputeFromLib
  
  subroutine singleImputationRound(library, c, errorAllow, minPresent, minOverlap, minHapFreq)
    use HaplotypeLibraryDefinition
    use CoreDefinition
    use GenotypeModule
    use HaplotypeModule
    
    type(HaplotypeLibrary), intent(in) :: library
    class(Core) :: c
    integer, intent(in) :: errorAllow, minPresent, minOverlap, minHapFreq
    
    integer :: i
    type(Genotype), pointer :: geno
    type(Haplotype) :: hap1, hap2, comp
    type(Haplotype) :: newHap, libHap
    logical :: hap1changed, hap2changed
    integer, allocatable, dimension(:) :: candHapsPat, candHapsMat !, compatHaps
    integer, pointer, dimension(:,:) :: candPairs
    
    do i = 1, c%getNAnisG()
      geno => c%coreGenos(i)
      hap1 = c%phase(i,1)
      hap2 = c%phase(i,2)
      
      hap1changed = .false.
      hap2changed = .false.
      
      if (mod(i,2000) == 0) then
	print *, "Done", i
      end if
      
      if ((.not. hap1%fullyPhased()) .or. (.not. hap2%fullyPhased())) then
!	compatHaps => library % getCompatHapsFreq(geno,minHapFreq,errorAllow)

	!! Can save some time at various points in this function by checking for fully phased
!	candHapsPat => library % limitedMatchWithErrorAndMinOverlap(c % phase(i, 1), ErrorAllow, minOverlap, compatHaps)
!	candHapsMat => library % limitedMatchWithErrorAndMinOverlap(c % phase(i, 2), ErrorAllow, minOverlap, compatHaps)
	candHapsPat = library % matchWithErrorAndMinOverlap(c % phase(i, 1), ErrorAllow, minOverlap)
	candHapsMat = library % matchWithErrorAndMinOverlap(c % phase(i, 2), ErrorAllow, minOverlap)
	
	if ((size(candHapsPat) > 0) .and. (size(candHapsMat) == 0)) then
	  libHap = getLibraryHap(library, candHapsPat)
	  newHap = newHaplotypeHaplotype(hap1)
	  call newHap%setFromOther(libHap)
	  hap1changed = .not. newHap%equalHap(hap1)
	  c%phase(i,1) = newHap
	  if (hap1changed .and. (newHap%numberNotMissing() >= minpresent)) then
	    call newHaplotype(c, i, 1, library, minoverlap, errorallow)
	  end if
	end if
	
	if ((size(candHapsPat) == 0) .and. (size(candHapsMat) > 0)) then
	  libHap = getLibraryHap(library, candHapsMat)
	  newHap = newHaplotypeHaplotype(hap2)
	  call newHap%setFromOther(libHap)
	  hap2changed = .not. newHap%equalHap(hap2)
	  c%phase(i,2) = newHap
	  if (hap2changed .and. (newHap%numberNotMissing() >= minpresent)) then
	    call newHaplotype(c, i, 2, library, minoverlap, errorallow)
	  end if
	end if
	
	if ((size(candHapsPat) > 0) .and. (size(candHapsMat) > 0)) then
	  candPairs => getCompatPairsWithError(library, geno, ErrorAllow, CandHapsPat, CandHapsMat, c%getNAnisG())
	  
	  if (size(CandPairs,1) > 0) then
	    libHap = getLibraryHap(library, candPairs(:,1))
	    newHap = newHaplotypeHaplotype(hap1)
	    call newHap%setFromOther(libHap)
	    hap1changed = .not. newHap%equalHap(hap1)
	    c%phase(i,1) = newHap
	    if (hap1changed .and. (newHap%numberNotMissing() >= minpresent)) then
	      call newHaplotype(c, i, 1, library, minoverlap, errorallow)
	    end if

	    libHap = getLibraryHap(library, candPairs(:,2))
	    newHap = newHaplotypeHaplotype(hap2)
	    call newHap%setFromOther(libHap)
	    hap2changed = .not. newHap%equalHap(hap2)
	    c%phase(i,2) = newHap
	    if (hap2changed .and. (newHap%numberNotMissing() >= minpresent)) then
	      call newHaplotype(c, i, 2, library, minoverlap, errorallow)
	    end if
	  end if
	end if

	hap2 = c%phase(i,2)
	comp = geno%complement(c%phase(i,1))
!	candHapsPat => library % limitedMatchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent, compatHaps)
	candHapsPat = library % matchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent)
	newHap = newHaplotypeHaplotype(hap2)
	if (size(candHapsPat) > 0) then
	  libHap = getLibraryHap(library, candHapsPat)
	  call newHap%setFromOther(libHap)
	else
	  call newHap%setFromOther(comp)
	end if
	c%phase(i,2) = newHap
	if (.not. newHap%equalHap(hap2) .and. (newHap%numberNotMissing() >= minpresent)) then
	  call newHaplotype(c, i, 2, library, minoverlap, errorallow)
	end if

	hap1 = c%phase(i,1)
	comp = geno%complement(c%phase(i,2))
!	candHapsMat => library % limitedMatchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent, compatHaps)
	candHapsMat = library % matchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent)
	newHap = newHaplotypeHaplotype(hap1)
	if (size(candHapsMat) > 0) then
	  libHap = getLibraryHap(library, candHapsMat)
	  call newHap%setFromOther(libHap)
	else
	  call newHap%setFromOther(comp)
	end if
	c%phase(i,1) = newHap
	if (.not. newHap%equalHap(hap1) .and. (newHap%numberNotMissing() >= minpresent)) then
	  call newHaplotype(c, i, 1, library, minoverlap, errorallow)
	end if
      end if
    end do
    
  end subroutine singleImputationRound
  
  subroutine complementStart(library, c, errorAllow, minPresent, minOverlap, minHapFreq)
    use HaplotypeLibraryDefinition
    use CoreDefinition
    use GenotypeModule
    use HaplotypeModule
    
    type(HaplotypeLibrary), intent(in) :: library
    class(Core) :: c
    integer, intent(in) :: errorAllow, minPresent, minOverlap, minHapFreq
    
    type(Haplotype) :: hap1, hap2, comp, libHap, newHap
    type(Genotype), pointer :: geno
    integer :: i
    integer, allocatable, dimension(:) :: candHapsPat, candHapsMat!, compatHaps
    
    do i = 1, c%getNAnisG()
      geno => c%coreGenos(i)
!      compatHaps => library % getCompatHapsFreq(geno,minHapFreq,errorAllow)
      hap2 = c%phase(i,2)
      comp = geno%complement(c%phase(i,1))
!      candHapsPat => library % limitedMatchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent, compatHaps)
      candHapsPat = library % matchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent)
      newHap = newHaplotypeHaplotype(hap2)
      if (size(candHapsPat) > 0) then
	libHap = getLibraryHap(library, candHapsPat)
	call newHap%setFromOther(libHap)
      else
	call newHap%setFromOther(comp)
      end if
      
      c%phase(i,2) = newHap
      if (.not. newHap%equalHap(hap2) .and. (newHap%numberNotMissing() >= minPresent)) then
	call newHaplotype(c, i, 2, library, minoverlap, errorallow)
      end if

      hap1 = c%phase(i,1)
      comp = geno%complement(c%phase(i,2))
!      candHapsMat => library % limitedMatchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent, compatHaps)
      candHapsMat = library % matchWithErrorAndMinOverlap(comp, ErrorAllow, minPresent)
      newHap = newHaplotypeHaplotype(hap1)
      if (size(candHapsMat) > 0) then
	libHap = getLibraryHap(library, candHapsMat)
	call newHap%setFromOther(libHap)
      else
	call newHap%setFromOther(comp)
      end if
      c%phase(i,1) = newHap
      if (.not. newHap%equalHap(hap1) .and. (newHap%numberNotMissing() >= minPresent)) then
	call newHaplotype(c, i, 1, library, minoverlap, errorallow)
      end if
    end do
    
  end subroutine complementStart
    
  function getLibraryHap(library, candHaps) result (libhap)
    use HaplotypeLibraryDefinition
    use HaplotypeModule
    
    type(HaplotypeLibrary), intent(in) :: library
    integer, dimension(:), intent(in) :: candHaps
    
    type(Haplotype) :: libhap
    
    libhap = newHaplotypeMissing(library%nSnps)
    
    ! Here for speed!
    if (size(CandHaps) == 1) then
      libhap = library%newstore(CandHaps(1))
    end if
    
    if (size(CandHaps) > 1) then
      call libhap%setZeroBits(library%oneZeroNoOnes(candHaps))
      call libhap%setOneBits(library%oneOneNoZeros(candHaps))
    end if
    
  end function getLibraryHap    
  
  function getCompatPairsWithError(library, geno, ErrorAllow, patLimit, matLimit, nAnisG) result(pairs)
    use GenotypeModule
    use HaplotypeLibraryDefinition
    
    class(HaplotypeLibrary) :: library
    type(Genotype), intent(in) :: geno
    integer, intent(in) :: ErrorAllow
    integer, dimension(:), intent(in) :: patLimit, matLimit
    integer, intent(in) :: nAnisG
    integer, dimension(:,:), pointer :: pairs
    
    integer, dimension(:,:), pointer :: tempPairs
    integer :: i, j, p, ii, jj

    allocate(tempPairs(nAnisG*2,2))
    
    p = 0
    i = 1
    do while ((i <= size(patLimit)) .and. ((p*p) <= (nAnisG - 1)))
      j = 1
      do while ((j <= size(matLimit)) .and. ((p*p) <= (nAnisG - 1)))
	ii = patLimit(i)
	jj = matLimit(j)
	if (geno%compatibleHaplotypes(library%newstore(ii), library%newstore(jj), ErrorAllow)) then
	  p = p + 1
	  tempPairs(p,1) = ii
	  tempPairs(p,2) = jj
	end if
	j = j + 1
      end do
      i = i + 1
    end do
    
    allocate(pairs(p,2))
    pairs = tempPairs(1:p,:)
    deallocate(tempPairs)
  end function getCompatPairsWithError
    

end module HaplotypeLibraryPhasing