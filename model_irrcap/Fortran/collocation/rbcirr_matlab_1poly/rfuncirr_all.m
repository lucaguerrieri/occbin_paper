
function Rvec = rfuncirr(avec,nodes,theta,P,kstart,kend)
%
% Rvec = rfuncirr(avec,nodes1,nodes2,theta,P,kstart,kend)
%
%
% The function returns the residuals of the Euler equation
% for the RBC model with capital irreversibility
%
% INPUTS
% avec
%
% now use two functions; one for the decision rule, the other for the
% Lagrangian multiplier. 
%
% The parameters for the chebychev polynomial functions aproximating the
% decision rule and the multiplier are stacked in avec.
% avec first stacks the parameters for the decision rule, state by state,
% then it stacks the parameters for the multiplier, state by state.
% The function is built to economize on parameters for the multiplier
% It is conjectured that the multipllier is never positive when the states
% are such that technology is non-negative.
% If the conjecture is not verified, the function issues a warning.
%
% nodes1, nodes2:
% In states when techonology is above steady state, 
% the Lagrange multiplier is conjectured to be zero
% hence, for collocation, we need just as many nodes as number of
% parameters for the function approximating the decisiont rule for capital
% Those nodes are stored in nodes1.
% In states when technology is below steady state, 
% the Lagrange multiplier may be positive, hence we need to parametrize
% both the decision rule and the multiplier. In those states we need
% twice as many nodes as when technology is above steady state. In this
% latter case, the nodes are stored in nodes2
%
% P
% transition probability matrix for the Markov process for technology
% 
% kstart, kend
% range over which the approximating polynomials will be defined


% model parameters
global BETA ALPHA DELTAK GAMMAC PSI PHI iss



nstates = length(theta);


order = length(avec)/(2*nstates)-1;


if length(nodes)==1
   Rvec = zeros(nstates,1);
else
   Rvec = zeros(length(avec),1); 
end
a1 = reshape(avec(1:((order+1)*nstates)),order+1,nstates);
a2 = reshape(avec(((order+1)*nstates+1):end),order+1,nstates);


rindex = 0;
for state = 1:nstates;
    
    a1_forstate = a1(:,state)';
    
    
    kap = nodes;
    
    
    npoints = length(kap);
    for inode = 1:npoints
        
        kap_prime = a1_forstate*chebypol(kap(inode),order,kstart,kend);
        lambdak = 0;
        if (kap_prime - (1-DELTAK)*kap(inode)) < PHI*iss
            kap_prime = PHI*iss+(1-DELTAK)*kap(inode);
            
            a2_forstate = a2(:,state)';
            lambdak = a2_forstate*chebypol(kap(inode),order,kstart,kend);
            
        end
        
        c = cfunc(theta(state),kap(inode),kap_prime);
        
        % compute exptectation term
        expectation_term = 0;
        for nextstate = 1:nstates
            a1_fornextstate = a1(:,nextstate)';
            
            lambdak_prime = 0;
            kap_2prime = a1_fornextstate*chebypol(kap_prime,order,kstart,kend);
            if (kap_2prime - (1-DELTAK)*kap_prime) < PHI*iss
                kap_2prime = PHI*iss+(1-DELTAK)*kap_prime;
                
                    a2_fornextstate = a2(:,nextstate)';
                    lambdak_prime = a2_fornextstate*chebypol(kap_prime,order,kstart,kend);
            end
           c_prime = cfunc(theta(nextstate),kap_prime,kap_2prime);
           this_term = c_prime^(-GAMMAC)*((1-DELTAK)-2*PSI*(kap_2prime/kap_prime-1)*(-kap_2prime/kap_prime^2)...
                       +ALPHA*exp(theta(nextstate))*kap_prime^(ALPHA-1))-(1-DELTAK)*lambdak_prime;
           expectation_term = expectation_term+P(state,nextstate)*this_term;
           
       end
       rindex = rindex+1;
       Rvec(rindex) = - c+((lambdak + (BETA*expectation_term))/(1+2*PSI*(kap_prime/kap(inode)-1)/kap(inode)))^(-1/GAMMAC);
       
   end
end


