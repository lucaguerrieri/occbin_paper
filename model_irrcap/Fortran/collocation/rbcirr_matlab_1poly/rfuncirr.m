function Rvec = rfuncirr(avec,nodes,theta,P,kstart,kend)

% model parameters
global BETA ALPHA DELTAK GAMMAC PSI PHI kss iss



npoints = length(nodes);
nstates = length(theta);

order = length(avec)/nstates-1;


kap = nodes;

R = zeros(npoints,nstates);
a = reshape(avec,order+1,nstates);


for state = 1:nstates;
    a_forstate = [a(:,state)]';
    for inode = 1:npoints
       efunc = a_forstate*chebypol(kap(inode),order,kstart,kend);
       
       lambdak = 0;
        kap_prime = invcfunc(theta(state),kap(inode),efunc^(-1/GAMMAC));
        c = cfunc(theta(state),kap(inode),kap_prime);      
        % next verify guess for lambdak
        if kap_prime -(1-DELTAK)*kap(inode)<PHI*iss   
            kap_prime = (1-DELTAK)*kap(inode)+PHI*iss;
            c=cfunc(theta(state),kap(inode),kap_prime);
            lambdak = c^(-GAMMAC)-efunc;
        end
       
       
       % compute exptectation term
       expectation_term = 0;
       for nextstate = 1:nstates
           a_fornextstate = [a(:,nextstate)]';
           efunc = a_fornextstate*chebypol(kap_prime,order,kstart,kend);
           
           lambdak_prime = 0;
           kap_2prime = invcfunc(theta(nextstate),kap_prime,efunc^(-1/GAMMAC));
           c_prime = cfunc(theta(nextstate),kap_prime,kap_2prime);
           if kap_2prime-(1-DELTAK)*kap_prime<PHI*iss
                 kap_2prime = (1-DELTAK)*kap_prime+PHI*iss;
                c_prime = cfunc(theta(nextstate),kap_prime,kap_2prime);
                lambdak_prime = c_prime^(-GAMMAC)-efunc;
            end
           
           this_term = c_prime^(-GAMMAC)*((1-DELTAK)-2*PSI*(kap_2prime/kap_prime-1)*(-kap_2prime/kap_prime^2)...
                       +ALPHA*exp(theta(nextstate))*kap_prime^(ALPHA-1))-(1-DELTAK)*lambdak_prime;
           expectation_term = expectation_term+P(state,nextstate)*this_term;
           
       end
       R(inode,state) = - c + ((lambdak+BETA*expectation_term)/(1+2*PSI*(kap_prime/kap(inode)-1)/kap(inode)))^(-1/GAMMAC);
       
   end
end

Rvec = R(:);
