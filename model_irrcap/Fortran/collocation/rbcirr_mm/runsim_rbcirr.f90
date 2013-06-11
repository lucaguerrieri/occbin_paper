PROGRAM runsim_rbcirr
! solves the RBC model by orthogonal collocation with 
! Chebyshev polynomials
! Remember to update the path for saving the results file at the end of the program -- the path is hardwired.



USE model_parameters
USE wrapper_vars


IMPLICIT NONE


DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: poly
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:) :: Pdat
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: avec, Rvec
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:) :: thetadat, nodesdat, nodes_mmdat
DOUBLE PRECISION :: Pd, thetad
INTEGER :: status
INTEGER :: info, m, n, indxi, nstates2
INTEGER :: version

EXTERNAL ::  rfuncirr, rfuncirr_wrapper, rfuncirr_wrapperls, rfuncirr_mm_wrapperls, rfunc, rfunc_wrapper, rfunc_wrapperls



BETA=0.96d0
ALPHA=0.33d0
DELTAK=0.10d0
GAMMAC=2.d0
RHOA = 0.9d0
PHI = 0.975d0
PSI = 0.d0        ! adjustment cost for capital
SIGMA_EPS = 0.013d0! deterministic steady state



kss = ((1.d0/BETA-1.d0+DELTAK)/ALPHA)**(1.d0/(ALPHA-1.d0))
css = -DELTAK*kss +kss**ALPHA
iss = DELTAK*kss
!uss = (css**(1.d0-GAMMAC)-1.d0)/(1.d0-GAMMAC)
!vss = uss/(1.d0-BETA)

version = 7

IF (version == 1) THEN
kstart = 0.9d0*kss
kend = 1.5d0*kss
order = 6

nstates = 51
PHI = 0.975d0
GAMMAC = 2.d0

ELSEIF (version == 2) THEN
kstart = 0.9d0*kss
kend = 1.5d0*kss
order = 6

nstates = 51
PHI = 0.975d0
GAMMAC = 3.d0

ELSEIF (version == 3) THEN
kstart = 0.9d0*kss
kend = 1.5d0*kss
order = 6

nstates = 51
PHI = 0.975d0
GAMMAC = 4.d0

ELSEIF (version == 4) THEN
kstart = 0.9d0*kss
kend = 1.5d0*kss
order = 6

nstates = 51
PHI = 0.975d0
GAMMAC = 5.d0

ELSEIF (version == 5) THEN
kstart = 0.7d0*kss
kend = 1.3d0*kss
order = 6

nstates = 51
PHI = 0.d0
GAMMAC = 2.d0

ELSEIF (version == 6) THEN
kstart = 0.9d0*kss
kend = 1.5d0*kss
order = 6

ELSEIF (version == 7) THEN
kstart = 0.9d0*kss
kend = 1.5d0*kss
order = 6
nstates = 51
PHI = 0.975d0
GAMMAC = 1.01d0

END IF



order_mm = order*1


nstates2 = 0



ALLOCATE(nodesdat(order+1),nodes_mmdat(order_mm+1), poly(order+1,1),Pdat(nstates,nstates), &
avec((order+1)*nstates), thetadat(nstates),Rvec((order+1)*nstates), STAT=status)

nodes => nodesdat
nodes_mm => nodes_mmdat
P => Pdat
theta => thetadat


CALL chebyroots((order+1),kstart,kend,nodes)
nodes_mm = nodes
CALL chebyroots((order_mm-order),kstart,kend,nodes_mm)
CALL linspace(kss*.975,kss*1.025,order_mm-order,nodes_mm(1:order_mm-order))
!nodes_mm(order_mm-order+1:order_mm+1)=nodes


WRITE(*,*) 'nodes'
CALL writevec_screen(nodes,(order+1),1,(order+1))

! get markov approximation to AR process
!CALL rouwenhorst(RHOA,SIGMA_EPS,nstates,theta,P)
!CALL markovappr(RHOA,SIGMA_EPS*SQRT(1.d0 - RHOA**2),2.d0,nstates,theta,P)

CALL markovappr(RHOA,SIGMA_EPS,3.d0,nstates,theta,P)


WRITE(*,*) 'P'
CALL writemat_screen(P,nstates,nstates,1,nstates,1,nstates)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SOLVE MODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first guess -- keep K constant at its SS value, regardless of theta
! the zeroth order coefficient is fixed at the SS value
! fix the first order coefficient at a small negative value -- if capital
! is above SS, capital is expected to be lower next period.

! first solve the deterministic model, then apply the solution from the
! deterministic model to the stochastic model


avec = 0.d0  ! first guess for Chebyshev coefficients

! first guesses for the deterministic model 


	avec(1) = css**(-GAMMAC)
	avec(2) = -0.1d0



Pd = 1.d0
thetad = 0.d0  ! one state of nature, set at steady state for technology
n = (order+1)
m = n
!CALL hybrd1(rfunc_wrapper,n,avec(1:(order+1)*nstates),Rvec(1:(order+1)*nstates),0.000000000001d0,info)

CALL rfuncirr(avec(1:order+1),n,nodes,order+1,thetad,1,Pd,kstart,kend,Rvec(1:(order+1)))
WRITE(*,*) 'Residuals For RBCIRR  Model at First Guess for SS'
WRITE(*,*) Rvec(1:order+1)


CALL lmdif1(rfunc_wrapperls,m,n,avec(1:(order+1)),Rvec(1:(order+1)),0.00000000001d0,info)


WRITE(*,*) 'info from optimizer'
WRITE(*,*) info


DO indxi =2,nstates
   avec((indxi-1)*(order+1)+1:(indxi)*(order+1))=avec(1:order+1)
END DO


WRITE(*,*) 'Residuals For RBC Model Witout Capital Constraint After  Call'
WRITE(*,*) Rvec(1:(order+1))

WRITE(*,*) 'avec'
WRITE(*,*) avec(1:(order+1))

! check output of rfuncirr
n = (order+1)*nstates
CALL rfuncirr(avec,n,nodes,order+1,theta,nstates,P,kstart,kend,Rvec)
WRITE(*,*) 'Residuals For RBC Model at First Guess'
WRITE(*,*) Rvec




!CALL hybrd1(rfuncirr_wrapper,n,avec,Rvec,0.0000000001d0,info)

n = (order+1)*nstates
m = n


CALL lmdif1(rfuncirr_wrapperls,m,n,avec,Rvec,0.000000001d0,info)


WRITE(*,*) 'info from optimizer'
WRITE(*,*) info


WRITE(*,*) 'avec'
WRITE(*,*) avec

WRITE(*,*) 'Residuals For Full Model'
WRITE(*,*) Rvec

!CALL rfuncirr_mm(avec,n,nodes_mm,order_mm+1,theta,nstates,P,kstart,kend,Rvec,order)

!WRITE(*,*) 'Integral'
!WRITE(*,*) Rvec



!n = (order+1)*nstates
!m = n

!CALL lmdif1(rfuncirr_mm_wrapperls,m,n,avec,Rvec,0.000000001d0,info)




CALL writevec(avec, (order+1)*nstates, 1, (order+1)*nstates,'/Users/Jason/Documents/MATLAB/consumption/occbin_work/occbin_work/model_irrcap/Fortran/collocation/rbcirr_matlab_1poly/load_avec_gammac_1.m')


DEALLOCATE(nodesdat, nodes_mmdat, poly,Pdat, avec, thetadat, Rvec, STAT=status)


END PROGRAM
