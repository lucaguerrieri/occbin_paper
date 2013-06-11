SUBROUTINE linspace(d1,d2,n,grid)

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: d1, d2
DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: grid

INTEGER :: indxi


grid(1) = d1
DO indxi= 0,n-2
   grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
END DO
grid(n) = d2

!MATLAB
!grid = [d1+(0:n-2)*(d2-d1)/(floor(n)-1) d2];


END SUBROUTINE