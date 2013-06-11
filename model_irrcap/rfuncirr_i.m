function [R kap_difference a_difference c] = rfuncirr_i(init,shock,irfshock,theta,P,modnam,modnamstar,constraint,constraint_relax,oobase_,Mbase_,oostar_,Mstar_)
        
% model parameters
global BETA ALPHA RHOA DELTAK GAMMAC PSI PHI kss iss
nperiods = 50;
tol0 = 1e-8;
maxiter = 10;

kpos = strmatch('k',Mbase_.endo_names,'exact');
ipos = strmatch('i',Mbase_.endo_names,'exact');
apos = strmatch('a',Mbase_.endo_names,'exact');
lambdakpos = strmatch('lambdak',Mbase_.endo_names,'exact');

kap = exp(init(kpos)+log(kss));  

nstates = length(theta);
 

%         [zdata zdataconcatenated] = ...
%          solve_one_constraint(modnam,modnamstar,...
%                               constraint, constraint_relax,...
%                               shock,irfshock,nperiods,tol0,maxiter);
        
[zdata zdataconcatenated] = ...
    solve_one_constraint_temp2(modnam,modnamstar,...
    constraint, constraint_relax,...
    shock,irfshock,nperiods,tol0,maxiter,init,oobase_,Mbase_,oostar_,Mstar_);


kap_prime=(1-DELTAK)*kap+exp(zdataconcatenated(1,ipos)+log(iss));
kap_difference = log(kap_prime)-log(kss);

lambdak = zdataconcatenated(1,lambdakpos);
a_difference = zdataconcatenated(1,apos);

c = cfunc(a_difference,kap,kap_prime);

newinit = zdataconcatenated(1,:)';


% compute exptectation term
expectation_term = 0;
state = max(find(min(abs(theta-newinit(apos)))==abs(theta-newinit(apos))));
for nextstate = 1:nstates
    
    next_shock = theta(nextstate)-RHOA*theta(state);
    
    [zdata zdataconcatenated] = ...
        solve_one_constraint_temp2(modnam,modnamstar,...
        constraint, constraint_relax,...
        next_shock,irfshock,nperiods,tol0,maxiter,newinit,oobase_,Mbase_,oostar_,Mstar_);
    
    % extract kap_2prime and lambdak_prime
    kap_2prime=(1-DELTAK)*kap_prime + exp(zdataconcatenated(1,ipos)+log(iss));
    lambdak_prime= zdataconcatenated(1,lambdakpos);
    
    c_prime = cfunc(theta(nextstate),kap_prime,kap_2prime);
    
    
    
    this_term = c_prime^(-GAMMAC)*((1-DELTAK)-2*PSI*(kap_2prime/kap_prime-1)*(-kap_2prime/kap_prime^2)...
        +ALPHA*exp(theta(nextstate))*kap_prime^(ALPHA-1))-(1-DELTAK)*lambdak_prime;
    expectation_term = expectation_term+P(state,nextstate)*this_term;
    
end
R = - c + ((lambdak+BETA*expectation_term)/(1+2*PSI*(kap_prime/kap-1)/kap))^(-1/GAMMAC);




