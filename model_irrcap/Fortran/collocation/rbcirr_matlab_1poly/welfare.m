% Computes value of using a particular policy forever 

% Inputs
% P: transition matrix for Z, nz by nz
% K: grid of values for which Kdec is defined
% Cdec: consumption function
% Z: grid of values for Z
% GAMMAC and BETA: model parameters
% crit: stopping criterion for convergence

function v=welfare(P,K,Z,Kdec,Cdec,GAMMAC,BETA,crit)

nz = numel(Z);
nk = numel(K);

if GAMMAC==1
U = log(Cdec);
else
U = (Cdec.^(1-GAMMAC)/(1-GAMMAC));
end

newV=U/(1-BETA);
EVPRIME=newV;

iHOWARD=1;
maxdiff=10;

while maxdiff>crit && iHOWARD<500
  
  disp([ maxdiff ])

  oldV = newV;
  EV=P*newV ;
  for iz=1:nz
     EVPRIME(iz,:)=interp1(K',EV(iz,:),Kdec(iz,:)','*spline')';
%     EVPRIME(iz,:)=qinterp1(K',EV(iz,:),Kdec(iz,:)')';
  end
  newV = U + BETA * EVPRIME ;
  
  iHOWARD=iHOWARD+1;
  maxdiff=max(abs(newV(:)-oldV(:)));
  
end

v=newV;