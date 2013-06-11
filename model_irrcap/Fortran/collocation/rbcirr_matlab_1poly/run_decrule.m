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

RHOA = 0.9;

PSI = 0;        % adjustment cost for capital
SIGMA_EPS = 0.013;

% deterministic steady state
kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;




version =1;
if version ==1
    
    
    % also choose k =3.44 z=-6%
    % and choos k +10% z = 4%
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 51; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .9*kss;
    kend = 1.5*kss;
    
    GAMMAC = 2;
    PHI = 0.975;
    
   
    
elseif version==2
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 51; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .9*kss;
    kend = 1.5*kss;
    
    GAMMAC = 3;
    PHI = 0.975;
    
   
elseif version == 3
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 51; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .9*kss;
    kend = 1.5*kss;
    
    GAMMAC = 4;
    PHI = 0.975;
    
   
    
elseif version == 4
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 51; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .9*kss;
    kend = 1.5*kss;
    
    GAMMAC = 5;
    PHI = 0.975;
    
  
    
elseif version == 5
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 51; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .7*kss;
    kend = 1.3*kss;
    
    GAMMAC = 2;
    PHI = 0;
    
      
elseif version == 6
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 3; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .9*kss;
    kend = 1.5*kss;
    
    GAMMAC = 10;
    PHI = 0.975;
    SIGMA_EPS = 0.03;
    
   
elseif version == 7
    order = 6;  % this is the order of the approximating polynomial
    % choose an even number so that the constant for the
    % polynomial family matches the steady state
    
    
    nstates = 51; % this is the number of states in the Markov approximation
    % to the AR(1) process for technology
    
    % choose an odd number to include the SS
    
    kstart = .9*kss;
    kend = 1.5*kss;
    
    GAMMAC = 1.01;
    PHI = 0.975;
    
end

nstates2 = 0;

nnegs = nstates;



if nstates>1
    %[theta, P] = rouwenhorst(RHOA,SIGMA_EPS,nstates);
    
    [P,theta]=markovappr(RHOA,SIGMA_EPS,3,nstates);
    
    
else
    P = 1; theta=0;
end



mat = zeros((order+1)*nstates,1);
if version == 1
    load_avec_gammac_2
elseif version == 2
    load_avec_gammac_3
elseif version ==3
    load_avec_gammac_4
elseif version == 4
    load_avec_gammac_5
elseif version == 5
    load_avec_gammac_2_phi_0
    
elseif version == 6
    load_avec_test
end
A_sol = mat;



avec = A_sol;
A_sol = reshape(A_sol(1:(order+1)*nstates),order+1,nstates);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot model variables for a given markov chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

kgrid = linspace(kstart,kend,800);

[kdec cdec]=decrule(kgrid, theta, A_sol, order, kstart, kend);

crit = 1e-8;



v=welfare(P,kgrid,theta,kdec,cdec,GAMMAC,BETA,crit)
