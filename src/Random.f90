module Random
  
  implicit none
  
contains
  
subroutine RandomOrder(order, n, start, idum)
    !     Generate a random ordering of the integers 1 ... n.

    integer, INTENT(IN) :: n, start
    integer, allocatable, INTENT(inout) :: order(:)
    integer :: idum

    !     Local variables
    integer :: i, j, k
    double precision :: wk


    do i = 1, n
	order(i) = start - 1 + i
    end do

    !     Starting at the end, swap the current last indicator with one
    !     randomly chosen from those preceeding it.

    do i = n, 2, -1
	wk = ran1(idum)
	j = 1 + i * wk
	if (j < i) then
	    k = order(i)
	    order(i) = order(j)
	    order(j) = k
	end if
    end do

    RETURN
end subroutine RandomOrder



FUNCTION ran1(idum)
    ! This Function returns a uniform random deviate between 0.0 and 1.0.
    ! Set IDUM to any negative value to initialize or reinitialize the sequence.
    !MODIFIED FOR REAL
    INTEGER idum, IA, IM, IQ, IR, NTAB, NDIV
    DOUBLE PRECISION ran1, AM, EPS, RNMX
    PARAMETER (IA = 16807, IM = 2147483647, AM = 1./IM, IQ = 127773, IR = 2836, NTAB = 32, NDIV = 1 + (IM - 1)/NTAB, EPS = 1.2e-7, &
    RNMX = 1. - EPS)
    INTEGER j, k, iv(NTAB), iy
    SAVE iv, iy
    DATA iv /NTAB * 0/, iy /0/
    IF (idum .le. 0.or.iy .eq. 0) then
	idum = max(-idum, 1)
	DO 11 j = NTAB + 8, 1, -1
	    k = idum/IQ
	    idum = IA * (idum - k * IQ) - IR * k
	    IF (idum .lt. 0) idum = idum + IM
	    IF (j .le. NTAB) iv(j) = idum

	    11 CONTINUE
	    iy = iv(1)
	END IF
	k = idum/IQ
	idum = IA * (idum - k * IQ) - IR * k
	IF (idum .lt. 0) idum = idum + IM
	j = 1 + iy/NDIV
	iy = iv(j)
	iv(j) = idum
	ran1 = min(AM * iy, RNMX)
	RETURN
END function ran1

end module Random