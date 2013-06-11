SUBROUTINE rfunc_wrapperls(m,n,x,fvec,iflag)

USE wrapper_vars

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, m
DOUBLE PRECISION, DIMENSION(n) :: x, fvec
DOUBLE PRECISION, DIMENSION(1,1) :: Pd
DOUBLE PRECISION, DIMENSION(1) :: thetad
INTEGER, INTENT(INOUT) :: iflag
EXTERNAL :: rfunc

Pd = 1.d0
thetad = 0.d0
!iflag =1

CALL rfuncirr(x,(order+1),nodes,order+1,thetad,1,Pd,kstart,kend,fvec)

END SUBROUTINE
