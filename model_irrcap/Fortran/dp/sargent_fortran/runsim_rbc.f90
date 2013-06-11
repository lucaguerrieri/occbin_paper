PROGRAM runsim_rbc
! solves the RBC model with irreversible capital 
! using dynamic programming


IMPLICIT NONE

DOUBLE PRECISION :: kss, css, iss, uss, vss
DOUBLE PRECISION :: lbound, ubound
INTEGER :: nstates, npoints, state, pointi, pointj, maxiter, i,j,k,kk, ipoint, jstate, howard_iter
INTEGER :: nlinvars, status
DOUBLE PRECISION :: stepsize, vdiff, maxdiff, shock, candidate
DOUBLE PRECISION :: inf_level
DOUBLE PRECISION,  DIMENSION(:), ALLOCATABLE :: kap, theta, maximand,expectation

DOUBLE PRECISION,  DIMENSION(:,:), ALLOCATABLE :: P, v, vprime, vnew, decrulea, decruleb, init, decisions
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: R
DOUBLE PRECISION, DIMENSION(1) :: tmp_1_1
INTEGER :: inttmp_1

DOUBLE PRECISION :: BETA, ALPHA, DELTAK, GAMMAC, RHOA, PSI, PHI, SIGMA_EPS

INTEGER, DIMENSION(:,:), ALLOCATABLE :: u
INTEGER :: kpos, vpos, apos

CHARACTER(len=98) :: path

BETA=0.96d0
ALPHA=0.33d0
DELTAK=0.10d0
GAMMAC=1.0001d0
RHOA = 0.9d0
PHI = 0.975d0
PSI = 0.d0        ! adjustment cost for capital
SIGMA_EPS = 0.01d0! deterministic steady state




kss = ((1.d0/BETA-1.d0+DELTAK)/ALPHA)**(1.d0/(ALPHA-1.d0))
css = -DELTAK*kss +kss**ALPHA
iss = DELTAK*kss
uss = (css**(1.d0-GAMMAC)-1.d0)/(1.d0-GAMMAC)
vss = uss/(1.d0-BETA)


lbound = 0.9d0*kss
ubound = 1.2d0*kss
npoints = 2001
nstates = 51
inf_level = 1000000000000.D0
howard_iter = 50  ! set to 0 to turn of Howard improvement algorithm


maxiter = 1000
vdiff = 100000
maxdiff = 10.d0**-7.d0

! at work, from windows -- change size to 46
!path = 'g:\ofs\research2\guerrieri\matlab\skander_dyn\'
! or from linux -- change size to 44
!path = '/ofs/research2/Guerrieri/matlab/skander_dyn/'

! on mac, 98
path = '/Users/Jason/Documents/MATLAB/consumption/occbin_20121229/model_irrcap/Fortran/dp/sargent_fortran/'

! number of variable in linear decision rule saved from matlab
nlinvars = 7 

! assume there is only one shock in the linear model
! then the dimensions of decrulea are nlinvars by nlinvars. the dimensions of decruleb are nlinvars by 1.
kpos = 2
apos = 6
vpos = 7

stepsize = (ubound-lbound)/(DBLE(npoints-1))


ALLOCATE(kap(npoints),P(nstates,nstates),theta(nstates),R(npoints,npoints,nstates), STAT=status)

WRITE(*,*) status
ALLOCATE(v(npoints,nstates), vprime(nstates,npoints), vnew(npoints,nstates))
ALLOCATE(decrulea(nlinvars,nlinvars),decruleb(nlinvars,1),maximand(npoints))
ALLOCATE(init(nlinvars,1), decisions(nlinvars,1), expectation(npoints), u(npoints,nstates))



!kap = (lbound:stepsize:ubound)
CALL linspace(lbound,ubound,npoints,kap)

!WRITE(*,*) 'kap'
!WRITE(*,*) (kap(pointi),pointi=1,npoints)

! get markov approximation to AR process

CALL rouwenhorst(RHOA,SIGMA_EPS,nstates,theta,P)

!WRITE(*,*) 'P'
!DO pointi = 1,nstates
!WRITE(*,*) (P(pointi,pointj),pointj=1,nstates)
!END DO


DO state = 1,nstates;
  DO pointj = 1,npoints
!$OMP PARALLEL DO
    DO pointi = 1,npoints
      R(pointi,pointj,state) = -kap(pointj)+(1.d0-DELTAK)*kap(pointi) - PSI*(kap(pointj)/kap(pointi)-1.d0)**2.d0 + EXP(theta(state))*kap(pointi)**ALPHA;
      IF ((candidate>0.d0) .AND. (kap(pointj)-(1.d0-DELTAK)*kap(pointi)>= PHI*iss)) THEN
        R(pointi,pointj,state) = ((R(pointi,pointj,state))**(GAMMAC-1.d0)-1.d0)/(GAMMAC-1.d0)
		!WRITE(*,*) R(pointi,pointj,state)
      ELSE
        R(pointi,pointj,state) = -inf_level
      END IF
    END DO
  END DO
END DO

!WRITE(*,*) 'First row of R for state 1'
!WRITE(*,*) (R(1,pointj,1),pointj=1,npoints)


! this is a mixture of the algorithm on page 100 of Sargent and Ljungvist
! and of the stuff on page 413 of Judd.

! choose initial guess for v based on the linear solution
! load linear solution
CALL loaddata(decrulea,nlinvars,nlinvars,path//'decrulea.dat')
CALL loaddata(decruleb,nlinvars,1,path//'decruleb.dat')

!WRITE(*,*) 'Decrulea'
!DO pointi = 1,nlinvars
!  WRITE(*,*) (decrulea(pointi,pointj),pointj=1,nlinvars)
!END DO
!WRITE(*,*) 'Decruleb'
!DO pointi = 1,nlinvars
!   WRITE(*,*) (decruleb(pointi,1))
!END DO


init = 0.d0


DO jstate = 1,nstates
  init(apos,1) = theta(jstate);
  shock = theta(jstate)-RHOA*theta(jstate)
  DO ipoint = 1,npoints
     init(kpos,1) = log(kap(ipoint))-log(kss)
	 	 
     decisions = matmul(decrulea,init)+decruleb*shock
     v(ipoint,jstate) = EXP(decisions(vpos,1)+log(vss))
  END DO
END DO

!WRITE (*,*) (v(1,ipoint), ipoint=1,nstates)

vnew = v;





i = 0;
DO WHILE ((vdiff>maxdiff) .AND. (i<=maxiter))
  i = i+1
  
  WRITE(*,*) 'iteration'
  WRITE(*,*) i
  
  vprime = TRANSPOSE(v);

  DO j = 1,nstates
    expectation = BETA*MATMUL(P(j,:),vprime)

    DO k=1,npoints
          maximand = R(k,:,j) + expectation;
	  !WRITE(*,*) 'Maximand'
	  !wRITE(*,*) (maximand(ipoint), ipoint=1,npoints)
	  vnew(k,j)= MAXVAL(maximand,1)
	  !inttmp_1 = MAXLOC(maximand,1)
	  u(k,j) = MAXLOC(maximand,1)
     END DO
  END DO
  v = vnew;
  !WRITE(*,*) (v(3,ipoint), ipoint=1,nstates)
  
  !WRITE(*,*) 'u'
  !WRITE(*,*) (u(3,ipoint), ipoint=1,nstates)
  
  DO kk=1,howard_iter

        DO j = 1,nstates
!$OMP PARALLEL DO       
            DO k = 1,npoints	     
                 vnew(k,j) = R(k,u(k,j),j)+ SUM(BETA*MATMUL(P(j,:),RESHAPE(v(u(k,j),:),(/nstates,1/))))
            END DO
        END DO
        v = vnew;
  END DO

  vdiff = MAXVAL(MAXVAL(ABS(v-TRANSPOSE(vprime)),1),1)
  WRITE(*,*) 'Max change in V is ' 
  WRITE(*,*) vdiff
END DO


  !WRITE(*,*) 'u'
  !WRITE(*,*) (u(3,ipoint), ipoint=1,nstates)


CALL writemat(v, npoints, nstates, 1, npoints, 1, nstates, path//'vmat1.m')
CALL writemat_int(u, npoints, nstates, 1, npoints, 1, nstates, path//'umat1.m')

DEALLOCATE(kap,P,theta,R,v, vprime, vnew, decrulea, decruleb, init, decisions, expectation, u, maximand)

END PROGRAM
