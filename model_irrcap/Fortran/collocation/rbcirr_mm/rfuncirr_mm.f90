SUBROUTINE rfuncirr_mm(avec,nparams,nodes_mm,npoints_mm,theta,nstates,P,kstart,kend,integral,order)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nparams, npoints_mm, nstates, order

DOUBLE PRECISION, DIMENSION(nparams), INTENT(IN) :: avec
DOUBLE PRECISION, DIMENSION(npoints_mm), INTENT(IN) :: nodes_mm
DOUBLE PRECISION, DIMENSION(nstates), INTENT(IN) :: theta
DOUBLE PRECISION, DIMENSION(nstates,nstates), INTENT(IN) :: P
DOUBLE PRECISION, INTENT(IN) :: kstart, kend
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Rvec, polmat
DOUBLE PRECISION, DIMENSION((order+1)*nstates), INTENT(OUT) :: integral

INTEGER :: i, j, k, t, status


ALLOCATE(Rvec(npoints_mm*nstates),polmat(npoints_mm), STAT=status)


CALL rfuncirr(avec,nparams,nodes_mm,npoints_mm,theta,nstates,P,kstart,kend,Rvec)


polmat = 1.d0
integral = 0.d0

j = 0; 
DO i = 1,order+1
	DO k = 1,nstates
		j = j+1;
		DO t = 1,npoints_mm
			integral(j) = integral(j) + Rvec(npoints_mm*(k-1)+t)*polmat(t)
		END DO
	END DO
	DO t = 1,npoints_mm
		polmat(t) = polmat(t)*nodes_mm(t)
	END DO
END DO



DEALLOCATE(Rvec, polmat, STAT=status)


END SUBROUTINE