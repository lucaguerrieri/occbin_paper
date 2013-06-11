

setpathdynare4

global M_ oo_ 
global BETA ALPHA RHOA DELTAK GAMMAC PSI PHI kss iss

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
modnam = 'dynrbc';
modnamstar = 'dynrbcirr';


if ispc                             
!copy paramfile_gammac_2.m paramfile.m
else
!cp paramfile_gammac_2.m paramfile.m    
end
paramfile
def_parm
SIGMA_EPS = 0.013;              
% express the occasionally binding constraint
% in linearized form
% one can use any combination of endogenous variables and parameters 
% declared in the the dynare .mod files
% constraint1 defines the first constraint
% if constraint1 is true, solution switches to model2
% but if constraint_relax is true, solution reverts to model1


constraint = 'i<log(PHI)';
constraint_relax ='lambdak<0';

% Pick innovation for IRFs
irfshock =char('erra');      % label for innovation for IRFs
                             % needs to be an exogenous variable in the
                             % dynare .mod files

maxiter = 10;     
tol0 = 1e-8;


nperiods_out = 1;

                          

 kstart = .9*kss;
 kend = 1.5*kss;
 nkgrid = 100;
 kgrid = linspace(kstart,kend,nkgrid);
 
 nstates = 51;
 [P, theta] =markovappr(RHOA,SIGMA_EPS,3,nstates);

 
 
 % Solve model, generate model IRFs
[zdata zdataconcatenated oobase_ Mbase_ oostar_ Mstar_] = ...
         solve_one_constraint_temp1(modnam,modnamstar,...
                              constraint, constraint_relax,...
                              0,irfshock,nperiods_out,tol0,maxiter);

% median technology 
 kpos = strmatch('k',M_.endo_names);
 init = 0*zdata(1,:)';
 Rvec = zeros(nkgrid);
 Cvec = zeros(nkgrid);
 
 for thisk = 1:nkgrid
    
     init(kpos) = log(kgrid(thisk))-log(kss);
     [Rvec(thisk) dummy1 dummy2 Cvec(thisk)]= rfuncirr_i(init,0,irfshock,theta,P,modnam,modnamstar,constraint,constraint_relax,oobase_,Mbase_,oostar_,Mstar_);
 end

l1=plot((log(kgrid)-log(kss))*100,log10(abs(Rvec./Cvec)),'k'); 
set(l1,'LineWidth',2);

%%%%%%%%%%%%
if ispc                             
!copy paramfile_gammac_2_phi_0.m paramfile.m
else
!cp paramfile_gammac_2_phi_0.m paramfile.m    
end
paramfile
def_parm



[zdata zdataconcatenated oobase_ Mbase_ oostar_ Mstar_] = ...
         solve_one_constraint_temp1(modnam,modnamstar,...
                              constraint, constraint_relax,...
                              0,irfshock,nperiods_out,tol0,maxiter);

 
 % median technology 
 kpos = strmatch('k',M_.endo_names);
 init = 0*zdata(1,:)';
 Rvec = zeros(nkgrid);
 Cvec = zeros(nkgrid);
 
 for thisk = 1:nkgrid
    
     init(kpos) = log(kgrid(thisk))-log(kss);
     [Rvec(thisk) dummy1 dummy2 Cvec(thisk)]= rfuncirr_i(init,0,irfshock,theta,P,modnam,modnamstar,constraint,constraint_relax,oobase_,Mbase_,oostar_,Mstar_);
 end

l1=plot((log(kgrid)-log(kss))*100,log10(abs(Rvec./Cvec)),'b-.');
set(l1,'LineWidth',2);
legend('Collocation Solution','Piece-Wise Linear Solution','Log-Linear Solution (model without investment constraint)')