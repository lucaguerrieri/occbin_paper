SUBROUTINE rfunc_wrapper(n,x,fvec,iflag)

USE wrapper_vars

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, DIMENSION(n) :: x, fvec
INTEGER, INTENT(INOUT) :: iflag
EXTERNAL :: rfunc

!iflag =1

CALL rfunc(x,(order+1)*nstates,nodes,order+1,theta,nstates,P,kstart,kend,fvec)

END SUBROUTINE
