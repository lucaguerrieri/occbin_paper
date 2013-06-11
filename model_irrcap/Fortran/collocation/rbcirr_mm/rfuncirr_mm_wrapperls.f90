SUBROUTINE rfuncirr_mm_wrapperls(m,n,x,fvec,iflag)

USE wrapper_vars

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, m
DOUBLE PRECISION, DIMENSION(n) :: x, fvec
INTEGER, INTENT(INOUT) :: iflag
EXTERNAL :: rfunc

!iflag =1

CALL rfuncirr_mm(x,(order+1)*nstates,nodes_mm,(order_mm+1),theta,nstates,P,kstart,kend,fvec,order)


!CALL rfuncirr(x,(order+1)*nstates,nodes_mm,(order_mm+1),theta,nstates,P,kstart,kend,fvec)


END SUBROUTINE
