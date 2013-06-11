SUBROUTINE rouwenhorst(rho,sigma_eps,n,zgrid,P)

!rouwenhorst.f90
!
![zgrid, P] = rouwenhorst(rho, sigma_eps, n)
!
! rho is the 1st order autocorrelation
! sigma_eps is the standard deviation of the error term
! n is the number of points in the discrete approximation
!
! see "Finite State Markov-chain Approximations to Highly Persistent
! Processes"


IMPLICIT NONE

INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: rho,sigma_eps
DOUBLE PRECISION, DIMENSION(n,n), INTENT(OUT) :: P
DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: zgrid

INTEGER :: i, status
DOUBLE PRECISION :: mu_eps, q, nu
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: P11, P12, P21, P22
EXTERNAL :: linspace

mu_eps = 0.d0

q = (rho+1.d0)/2.d0
nu = SQRT((DBLE(n)-1.d0)/(1.d0-rho**2.d0)) * sigma_eps
P = 0.d0

P(1,1) = q
P(1,2) = 1 - q
P(2,1) = 1 - q
P(2,2) = q

DO i=2,n-1
   WRITE(*,*) i
   ALLOCATE(P11(i+1,i+1), P12(i+1,i+1), P21(i+1,i+1), P22(i+1,i+1), STAT=status)
   
   P11 = 0.d0
   P12 = 0.d0
   P21 = 0.d0
   P22 = 0.d0
 
   P11(1:i,1:i) = P(1:i,1:i)
   P12(1:i,2:i+1) = P(1:i,1:i)
   P21(2:i+1,1:i) = P(1:i,1:i)
   P22(2:i+1,2:i+1) = P(1:i,1:i)
   WRITE(*,*) 'Assigned partial matrices'   

   P(1:i+1,1:i+1) = q*P11 + (1-q)*P12 + (1-q)*P21 + q*P22

   P(2:i,:) = P(2:i,:)/2.d0

   DEALLOCATE(P11,P12,P21,P22, STAT = status)
END DO

WRITE(*,*) 'About to call linspace'
CALL linspace(mu_eps/(1.d0-rho)-nu,mu_eps/(1.d0-rho)+nu,n,zgrid)

END SUBROUTINE