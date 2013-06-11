SUBROUTINE  SORT(x,n)

IMPLICIT NONE


INTEGER, INTENT(IN) :: n

DOUBLE PRECISION, DIMENSION(n), INTENT(INOUT) :: x
DOUBLE PRECISION :: xtemp
	
INTEGER :: iindx,loc, status

WRITE(*,*) 'positions'
DO iindx = 1,n-1
   loc = MINLOC(x(iindx:n),1)+iindx-1
   WRITE(*,*) loc
	  
   xtemp = x(iindx)
   x(iindx) = x(loc)
   x(loc) = xtemp
   
END DO

WRITE(*,*) 'x'
WRITE(*,*) x

END SUBROUTINE