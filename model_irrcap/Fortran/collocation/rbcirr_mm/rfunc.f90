SUBROUTINE rfunc(avec,nparams,kap,npoints,theta,nstates,P,kstart,kend,Rvec)

! model parameters
! contains BETA ALPHA DELTAK SIGMAC GAMMA SIGMAZERO SIGMAD
USE model_parameters

IMPLICIT NONE

INTEGER, INTENT(IN) :: nparams, npoints, nstates

DOUBLE PRECISION, DIMENSION(nparams), INTENT(IN) :: avec
DOUBLE PRECISION, DIMENSION(npoints), INTENT(IN) :: kap
DOUBLE PRECISION, DIMENSION(nstates), INTENT(IN) :: theta
DOUBLE PRECISION, DIMENSION(nstates,nstates), INTENT(IN) :: P
DOUBLE PRECISION, INTENT(IN) :: kstart, kend

DOUBLE PRECISION, DIMENSION(npoints*nstates), INTENT(OUT) :: Rvec

INTEGER :: order, state, nextstate, status, inode, posindx
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: a_forstate, a_fornextstate, poly
DOUBLE PRECISION :: kap_prime, kap_2prime, c, c_prime, expectation_term, this_term
EXTERNAL :: chebypol, cfunc


order = nparams/nstates-1

ALLOCATE( a_forstate(order+1), a_fornextstate(order+1), poly(order+1), STAT=status)

IF (status .NE. 0) THEN
   WRITE(*, *) 
   WRITE(*, *) '***In Call_....: Could not allocate memeory***'
END IF


posindx =0
DO state = 1,nstates
    a_forstate = avec((order+1)*(state-1)+1:(order+1)*state)

    DO inode = 1,npoints
       posindx = posindx+1
       CALL chebypol(kap(inode),order,kstart,kend,poly)
       CALL constructpol(a_forstate,poly,order,kap_prime)


       CALL cfunc(theta(state),kap(inode),kap_prime,c)

       ! compute exptectation term
       expectation_term = 0.0d0
       DO nextstate = 1,nstates
          
           a_fornextstate =  avec((order+1)*(nextstate-1)+1:(order+1)*nextstate)
           call chebypol(kap_prime,order,kstart,kend,poly)
           call constructpol(a_fornextstate,poly,order,kap_2prime)
           
           call cfunc(theta(nextstate),kap_prime,kap_2prime,c_prime)
           
           this_term = c_prime**(-GAMMAC)*((1-DELTAK)-2.d0*PSI*(kap_2prime/kap_prime-1.d0)*(-kap_2prime/kap_prime**2.d0)+ALPHA*exp(theta(nextstate))*kap_prime**(ALPHA-1.d0))
           expectation_term = expectation_term+P(state,nextstate)*this_term;
       END DO
       Rvec(posindx) = - c + ((BETA*expectation_term)/(1.d0+2.d0*PSI*(kap_prime/kap(inode)-1)/kap(inode)))**(-1/GAMMAC)

       
   END DO
END DO



DEALLOCATE(a_forstate, a_fornextstate, poly, STAT=status)

END SUBROUTINE