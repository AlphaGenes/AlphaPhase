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
      !Paternal Haps
      hap = c % phase(i, 1)
      if (hap%numberNotMissing() >= minoverlap) then
	call newHaplotype(c, i, 1, library, minoverlap, errorallow)
      endif

      !Maternal Haps
      hap = c % phase(i, 2)
      if (hap%numberNotMissing() >= minoverlap) then
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
    type(HaplotypeLibrary) :: library

    integer :: id
    integer, dimension(:), pointer :: ids
    type(Haplotype) :: hap, merged, libhap
    
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
	  merged = hap%mergeMod(library%newstore(id))
	  do i = 1, c%getNAnisG()
	    do j = 1, 2
	      if (c%hapanis(i,j) == id) then
		c%phase(i,j) = merged
	      end if
	    end do
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

	  do i = size(ids), 2, -1
	    library%hapfreq(id) = library%hapfreq(id) + library%hapfreq(ids(i))
	    call library%removehap(ids(i))
	    do j = 1, c%getNAnisG()
	      do k = 1, 2
		if (c%hapanis(j,k) > ids(i)) then
		  c%hapanis(j,k) = c%hapanis(j,k) - 1
		end if
	      end do
	    end do
	  end do
	else
	  libhap = getLibraryHap(library, ids)
	  call hap%setFromOtherIfMissing(libhap)
	  c%phase(animal,phase) = hap
	  id = library%addHap(hap)
	end if
      end if

      c%hapAnis(animal, phase) = id
  end subroutine newHaplotype
  
  subroutine updateHaplotype(c, library, id, hap, minoverlap, errorallow)
    use HaplotypeModule
    use CoreDefinition
    use HaplotypeLibraryDefinition
    
    class(Core) :: c
    integer, intent(in) :: id, minoverlap, errorallow
    type(HaplotypeLibrary) :: library
    
    type(Haplotype) :: hap, merged, libhap    
    integer :: i, j, k, n    
    integer, dimension(:), pointer :: ids

      library%newstore(id) = hap
      ids = library%matchWithErrorAndMinOverlap(hap,0,minoverlap)
      if (size(ids) > 1) then
	merged = newHaplotypeMissing(hap%length)
	do i = 1, size(ids)
	  merged = merged%mergeMod(library%newstore(ids(i)))
	end do
	
	if (merged%numberError() <= ErrorAllow) then
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

	  do i = size(ids), 1, -1
	    if (ids(i) /= id) then
	      library%hapfreq(id) = library%hapfreq(id) + library%hapfreq(ids(i))
	      call library%removehap(ids(i))
	      do j = 1, c%getNAnisG()
		do k = 1, 2
		  if (c%hapanis(j,k) > ids(i)) then
		    c%hapanis(j,k) = c%hapanis(j,k) - 1
		  end if
		end do
	      end do
	    end if
	  end do
	else
	  libhap = getLibraryHap(library, ids)
	  call hap%setFromOtherIfMissing(libhap)
	  
	  library%newstore(id) = hap
	  
	  do i = 1, c%getNAnisG()
	    do j = 1, 2
	      if (c%hapanis(i,j) == id) then
		c%phase(i,j) = hap
	      end if
	    end do
	  end do
	end if
      else
	do i = 1, c%getNAnisG()
	  do j = 1, 2
	    if (c%hapanis(i,j) == id) then
	      c%phase(i,j) = hap
	    end if
	  end do
	end do	
      end if
  end subroutine updateHaplotype
  
  
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
    
    integer :: i, j, oj
    type(Haplotype) :: newhap, comp
    
    do i = 1, c%getNAnisG()
      do j = 1, 2
	oj = 3 - j
	newhap = newHaplotypeHaplotype(c%phase(i,oj))
	comp = c%coregenos(i)%complement(c%phase(i,j))
	call newhap%setFromOtherIfMissing(comp)
	if (.not. newhap%equalHap(c%phase(i,oj))) then
	  if (c%hapAnis(i,oj) /= MissingHaplotypeCode) then
	    call updateHaplotype(c,library,c%hapAnis(i,oj),newhap,minOverlap,errorallow)
	  else
	    c%phase(i,j) = newhap
	    if (newhap%numberNotMissing() >= minOverlap) then
	      call newHaplotype(c, i, j, library, minOverlap, errorallow)
	    end if
	  end if
	end if
      end do
    end do
    
  end subroutine singleImputationRound
    
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