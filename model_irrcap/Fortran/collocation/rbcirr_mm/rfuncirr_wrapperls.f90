SUBROUTINE rfuncirr_wrapperls(m,n,x,fvec,iflag)

USE wrapper_vars

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, m
DOUBLE PRECISION, DIMENSION(n) :: x, fvec
INTEGER, INTENT(INOUT) :: iflag
EXTERNAL :: rfunc

!iflag =1

CALL rfuncirr(x,(order+1)*nstates,nodes,(order+1),theta,nstates,P,kstart,kend,fvec)

END SUBROUTINE
