SUBROUTINE chebypol(x,n,a,b,poly)

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: x, a, b
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, DIMENSION(n+1), INTENT(OUT) :: poly

DOUBLE PRECISION :: xrebased
INTEGER :: i
! n is order 
! [a b] is the range over which the approximation is desired
! x needs to be mapped back into -1 1, the domain of definition
! for Chebyshev polynomials.

! evaluates the n+1 cheby polynomials of nth order
! and stacks them in poly= [order 0 at x 
                           !order 1 at x ...]

xrebased = (x -a)/(b-a)*2.d0-1.d0


poly=0.0d0

poly(1) = 1.d0

IF (n>0) THEN
    poly(2) = xrebased
END IF 

IF (n>1) THEN
    DO i=3,n+1
        poly(i)=2.d0*xrebased*poly(i-1)-poly(i-2)
    END DO
END IF


END SUBROUTINE