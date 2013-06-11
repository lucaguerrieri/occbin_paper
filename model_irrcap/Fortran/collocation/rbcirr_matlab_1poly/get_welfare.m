% solves the RBC model by orthogonal collocation with 
% Chebyshev polynomials


clear

global BETA ALPHA DELTAK GAMMAC PSI PHI iss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BETA=0.96;
ALPHA=0.33;
DELTAK=0.10;
GAMMAC= 5;
RHOA = 0.9;
PHI = 0.975;
PSI = 0;        % adjustment cost for capital
SIGMA_EPS = 0.01;

% deterministic steady state
kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;
uss = (css^(1-GAMMAC)-1)/(1-GAMMAC);
vss = uss/(1-BETA);

% discretize the asset space
kstart = .8*kss;
kend = 1.2*kss;

% order = 6;  % this is the order of the approximating polynomial
%             % choose an even number so that the constant for the
%             % polynomial family matches the steady state
% 
% % get markov approximation to AR process
% nstates = 61; % this is the number of states in the Markov approximation
%               % to the AR(1) process for technology
%               
%               % choose an odd number to include the SS
%               
% nstates2 = 12;  

order = 5;  % this is the order of the approximating polynomial
            % choose an even number so that the constant for the
            % polynomial family matches the steady state

% get markov approximation to AR process
nstates = 51; % this is the number of states in the Markov approximation
              % to the AR(1) process for technology
              
              % choose an odd number to include the SS
              
nstates2 = 0;  

nnegs = nstates;
              
load_from_fortran_program = 1;  
             
if nstates>1
[theta, P] = rouwenhorst(RHOA,SIGMA_EPS,nstates);
%CALL markovappr(RHOA,SIGMA_EPS,3.d0,nstates,4,theta,P)

%[P,theta]=markovappr2(RHOA,SIGMA_EPS,3,nstates,nstates2);


else
P = 1; theta=0;
end

if ~load_from_fortran_program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SOLVE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first guess -- keep K constant at its SS value, regardless of theta
% the zeroth order coefficient is fixed at the SS value
% fix the first order coefficient at a small negative value -- if capital
% is above SS, capital is expected to be lower next period.
%
% first solve the deterministic model without capital irreversibility,
% then apply the solution from the deterministic model to the stochastic model
% without capital irreversibility
% and finally solve the model with capital irreversibility.


avec = zeros(order+1,1);  % first guess for Chebyshev coefficients

% first guesses for the deterministic model 
avec(1) = kss;
avec(2) = -0.001;

% set Markov prob matrix and state vector for deterministic case
Pd = 1; thetad = 0;  % one state of nature, set at steady state for technology

% options for fsolve
options = optimset('Display','Iter','MaxFunEvals',1e10,'MaxIter',1e5,'TolFun',1e-7,'Algorithm','trust-region-reflective');

nodes = chebyroots(order+1,kstart,kend);

% test residual function
% rfunc(avec,nodes,theta,P,kstart,kend)

% find deterministic solution of model without capital irreversibility constraint
A_sol_deterministic = ...
     fsolve(@(avec) rfunc(avec,nodes,thetad,Pd,kstart,kend), avec,options);

nnegs = length(find(theta<0));  % number of states in which technology is below SS
                                % assume and then verify that 
                                % capital constraint will not be binding 
                                % when tehcnology is above SS


avec = repmat(A_sol_deterministic,nstates,1); % stacks nstates copy of the deterministic solution
% find stochastic solution
A_sol_stoch = fsolve(@(avec) rfunc(avec,nodes,theta,P,kstart,kend), avec,options);


alambdafirstguess = zeros(order+1,1);
alambdafirstguess(1) =0.0 ;
alambdafirstguess(2) = -0.001;

nnegs = nstates;

% use deterministic solution as first guess for the stochastic solution
avec = [A_sol_stoch 
        repmat(alambdafirstguess,nstates,1)]; 
   
nodes = chebyroots(2*(order+1),kstart,kend);    

% find stochastic solution

A_sol = fsolve(@(avec) rfuncirr2(avec,nodes,theta,P,kstart,kend), avec,options);

% reshape A_sol so that each column contains the Chebyshev coefficients for
% a particular state of nature
else

   mat = zeros(2*(order+1)*nstates,1);
   load_avec;
   A_sol = mat;

end

avec = A_sol;
A1_sol = reshape(A_sol(1:(order+1)*nstates),order+1,nstates);
A2_sol = reshape(A_sol((order+1)*nstates+1:end),order+1,nstates);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot model variables for a given markov chain 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% generate path according to markov chain defined below
%chain_pos1 = markov_match_ar(RHOA,2*SIGMA_EPS,theta,1);
%chain_pos2 = markov_match_ar(RHOA,-2*SIGMA_EPS,theta,1);

%chain = [ones(9,1)*ceil(nstates/2) ; chain_pos1*ones(10,1); chain_pos2*ones(10,1); ceil(nstates/2)*ones(10,1)];

rand('state',1);
%chain = markov(P,10000,floor(nstates/2))';

%chain=markov_luca(P,10000,floor(nstates/2));



[chain1]=markov_match_ar(RHOA,theta(5),theta,50);
[chain2]=markov_match_ar(RHOA,theta(nstates-5+1),theta,50);
chain = [ceil(nstates/2)*ones(1,9) chain1 chain2];

%[chain1]=markov_match_ar(RHOA,theta(5),theta,50);
%chain = [chain1];




nperiods = length(chain);
k_path = zeros(nperiods+1,1);
c_path = zeros(nperiods+1,1);
lambdak_path = zeros(nperiods+1,1);
resid_path = zeros(nperiods+1,1);

k_path(1) = kss;
lambdak_path(1) = 0;
for period = 2:nperiods+1;
    waitbar(period/(nperiods+1))
    a1_state = [A1_sol(:,chain(period-1))]';
    k_path(period) = a1_state*chebypol(k_path(period-1),order,kstart,kend);
    if (k_path(period)-(1-DELTAK)*k_path(period-1))<PHI*iss
           k_path(period) = PHI*iss+(1-DELTAK)*k_path(period-1);
       
           a2_state = [A2_sol(:,chain(period-1))]';
           lambdak_path(period-1) = a2_state*chebypol(k_path(period-1),order,kstart,kend);
       
    end
    c_path(period-1) = cfunc(theta(chain(period-1)),k_path(period-1),k_path(period));
    expectation_term = 0;
    kap_prime = k_path(period);
    resids = rfuncirr2(avec,k_path(period-1),theta,P,kstart,kend);
    resid_path(period-1) = resids(chain(period-1))/css;
end




figure
subplot(3,2,1)
plot((k_path(2:end)-kss)/kss*100)
ylabel('Percent Dev. from SS')
title('Capital')

subplot(3,2,2)
i_path = k_path(2:end)-(1-DELTAK)*k_path(1:end-1);
plot(100*(i_path-iss)/iss)
ylabel('Percent Dev. from SS')
title('Investment')

subplot(3,2,3)
plot(lambdak_path)
ylabel('Level')
title('Lagrangian Multiplier on Kuhn Tucker Condition')

subplot(3,2,4)
theta_path = theta(chain)';
%c_path = cfunc(theta_path,k_path(1:end-1),k_path(2:end));
plot((c_path(1:end-1)-css)/css)
ylabel('Percent Dev. from SS')
title('Consumption')

subplot(3,2,5)
plot(theta_path)
ylabel('Percent Dev. from SS')
title('Technology')

subplot(3,2,6)
plot(abs(resid_path(1:end-1))*100)
title('Euler Residual (percent of SS consumption)')


% bridge between Luca and Matteo's notation

simC = c_path(1:end-1)';
simK = k_path(2:end)';
simI = i_path';
simZ = exp(theta_path(1:end)');
simLAMBDAK = lambdak_path(1:end-1)';

save ../../../luca_nl_sim simC simK simI simZ simLAMBDAK resid_path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot decision rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure
kgrid = [kstart*(1.025):.005:kend*(0.975)];
npoints = length(kgrid);

for ia=[ 1:nstates]
a1_state = [A1_sol(:,ia)]';

   for i=1:npoints
        kprime(i) =  a1_state*chebypol(kgrid(i),order,kstart,kend);
       
        if (kprime(i) - (1-DELTAK)*kgrid(i))<PHI*iss
            kprime(i) = (1-DELTAK)*kgrid(i)+PHI*iss;    
        end
   end
   h=plot((kgrid-kss)/kss*100,(kprime-kgrid)/kss*100); hold on; 
   set(h,'Color',[ia/nstates 0.5 1 ])
end
ylabel('Percent Dev. from Previous Period')
xlabel('Percent Dev. from SS')
title('Decision Rules')