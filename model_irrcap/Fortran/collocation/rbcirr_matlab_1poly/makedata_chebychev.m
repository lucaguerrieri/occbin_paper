function [simC simK simY simZ euler_err_nl ] = makedata_chebychev(A1_sol,A2_sol,chain,nperiods,order);

global PHI DELTAK PSI BETA ALPHA GAMMAC RHOA kss kstart kend iss theta nstates nnegs P avec

nperiods = length(chain);
k_path = zeros(nperiods+1,1);
c_path = zeros(nperiods,1);
lambdak_path = zeros(nperiods,1);
resid_path = zeros(nperiods,1);

k_path(1) = kss;
lambdak_path(1) = 0;
for period = 2:nperiods+1;
    a1_state = [A1_sol(:,chain(period-1))]';
    k_path(period) = a1_state*chebypol(k_path(period-1),order,kstart,kend);
    if (k_path(period)-(1-DELTAK)*k_path(period-1))<PHI*iss
       k_path(period) = PHI*iss+(1-DELTAK)*k_path(period-1);
       if theta(chain(period-1))<0
           a2_state = [A2_sol(:,chain(period-1))]';
           lambdak_path(period-1) = a2_state*chebypol(k_path(period-1),order,kstart,kend);
       else
           warning('Conjecture that lambdak is 0 when technology is above SS is not verified')
       end
    end
    c_path(period-1) = cfunc(theta(chain(period-1)),k_path(period-1),k_path(period));
    expectation_term = 0;
    kap_prime = k_path(period);
    for nextstate = 1:nstates
        a1_fornextstate = A1_sol(:,nextstate)';
        
        lambdak_prime = 0;
        kap_2prime = a1_fornextstate*chebypol(kap_prime,order,kstart,kend);
        if (kap_2prime - (1-DELTAK)*kap_prime) < PHI*iss
            kap_2prime = PHI*iss+(1-DELTAK)*kap_prime;
            if nextstate<=nnegs
                a2_fornextstate = A2_sol(:,nextstate)';
                lambdak_prime = a2_fornextstate*chebypol(kap_prime,order,kstart,kend);
            else
                warning('set lambdak_prime')
                lambdak_prime = 0.1;
            end
        end
        c_prime = cfunc(theta(nextstate),kap_prime,kap_2prime);
        this_term = c_prime^(-GAMMAC)*((1-DELTAK)-2*PSI*(kap_2prime/kap_prime-1)*(-kap_2prime/kap_prime^2)...
            +ALPHA*exp(theta(nextstate))*kap_prime^(ALPHA-1))-(1-DELTAK)*lambdak_prime;
        expectation_term = expectation_term+P(chain(period-1),nextstate)*this_term;
        
    end
    %resid_path(period-1) = (- c_path(period-1) ...
    %     +((lambdak_path(period-1) + (BETA*expectation_term))...
    %     /(1+2*PSI*(kap_prime/k_path(period-1)-1)/k_path(period-1)))^(-1/GAMMAC))/css;
    resids = rfuncirr(avec,k_path(period-1),k_path(period-1),theta,P,kstart,kend);
    resid_path(period-1) = resids(chain(period-1));
end


theta_path = theta(chain(1:end))';

simC=c_path(1:end)';
simK=k_path(2:end)';
simZ=exp(theta_path(1:end))';

simKlag = [ k_path(1:end-1) ]';

simY=simZ.*simKlag.^ALPHA;

euler_err_nl = resid_path(1:end);


% figure
% subplot(3,2,1)
% plot((k_path(2:end)-kss)/kss*100)
% ylabel('Percent Dev. from SS')
% title('Capital')
% 
% subplot(3,2,2)
% i_path = k_path(2:end)-(1-DELTAK)*k_path(1:end-1);
% plot(100*(i_path-iss)/iss)
% ylabel('Percent Dev. from SS')
% title('Investment')
% 
% subplot(3,2,3)
% plot(lambdak_path)
% ylabel('Level')
% title('Lagrangian Multiplier on Kuhn Tucker Condition')
% 
% subplot(3,2,4)
% theta_path = theta(chain)';
% %c_path = cfunc(theta_path,k_path(1:end-1),k_path(2:end));
% plot((c_path-css)/css*100)
% ylabel('Percent Dev. from SS')
% title('Consumption')
% 
% subplot(3,2,5)
% plot(theta_path*100)
% ylabel('Percent Dev. from SS')
% title('Technology')
% 
% subplot(3,2,6)
% plot(resid_path(1:end-1)*100)
% title('Euler Residual (consumption equivalent)')
% 
% save chebypaths k_path i_path c_path lambdak_path theta_path
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Plot decision rules
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% kgrid = [kstart:.005:kend];
% npoints = length(kgrid);
% 
% for ia=[ 1:nstates]
% a1_state = [A1_sol(:,ia)]';
% 
%    for i=1:npoints
%         kprime(i) =  a1_state*chebypol(kgrid(i),order,kstart,kend);
%        
%         if (kprime(i) - (1-DELTAK)*kgrid(i))<PHI*iss
%             kprime(i) = (1-DELTAK)*kgrid(i)+PHI*iss;    
%         end
%    end
%    h=plot((kgrid-kss)/kss*100,(kprime-kgrid)/kss*100); hold on; 
%    set(h,'Color',[ia/nstates 0.5 1 ])
% end
% ylabel('Percent Dev. from Previous Period')
% xlabel('Percent Dev. from SS')
% title('Decision Rules')




