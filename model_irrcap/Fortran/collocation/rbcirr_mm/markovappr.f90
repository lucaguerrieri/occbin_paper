SUBROUTINE markovappr(lambda,sigma,m,n,s,tran)

! translation of matlab routine
! function [Tran,s]=markovappr(lambda,sigma,m,N)
!
! the simple case of approximating first-order 
! autoregressive process with Markov chain
!
! y_t = lambda * y_(t-1) + u_t
!
! u_t is a Gaussian white noise process with standard deviation sigma.
!
! m determines the width of discretized state space, Tauchen uses m=3
! ymax=m*vary,ymin=-m*vary, ymax and ymin are two boundary points
!
! n is the number of possible states chosen to approximate
! the y_t process, usually N=9 should be fine
!
! n2 is the subset of nodes placed between the interval [-sigma sigma]
!
! Tran is the transition matrix of the Markov chain
!
! s is the discretized state space of y_t



IMPLICIT NONE

INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: m, lambda,sigma 

DOUBLE PRECISION, DIMENSION(n,n), INTENT(OUT) :: tran
DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: s

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: s1,s2

DOUBLE PRECISION :: stvy, ymax, ymin, w
INTEGER :: ymin2_pos, ymax2_pos, j, k, status

DOUBLE PRECISION, EXTERNAL :: normcdf


! discretize the state space


stvy = sqrt(sigma**2.d0/(1.d0-lambda**2.d0)); ! standard deviation of y_t

ymax = m*stvy                  ! upper boundary of state space
ymin = -ymax                  ! lower boundary of state space

WRITE(*,*) 'ymax'
WRITE(*,*) ymax
WRITE(*,*) 'ymin'
WRITE(*,*) ymin

w = (ymax-ymin)/DBLE(n-1)            ! length of interval 

CALL linspace(ymin,ymax,n,s)       



! calculate the transition matrix

DO j=1,n
   
   DO k=2,n-1
      
      tran(j,k)= normcdf(s(k)-lambda*s(j)+(s(k+1)-s(k))/2.d0,0.d0,sigma) &
         - normcdf(s(k)-lambda*s(j)-(s(k)-s(k-1))/2.d0,0.d0,sigma)
      
   END DO
   
   Tran(j,1) = normcdf(s(1)-lambda*s(j)+w/2.d0,0.d0,sigma)      ! only subtract half the interval on the right
   Tran(j,n) = 1.d0 - normcdf(s(N)-lambda*s(j)-w/2.d0,0.d0,sigma)  ! only subtract half the interval on the left
   
END DO


   
END SUBROUTINE   
   




