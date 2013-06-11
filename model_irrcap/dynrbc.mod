// variables
var c, k, klag, i, lambdak, elambdak, a, v;    
 
// innovations to shock processes
varexo erra;

// parameters
parameters ALPHA, DELTAK, BETA, GAMMAC, RHOA, PHI, PSI, kss, css, iss,
           uss, vss;

// parameter values
// defined in external files
// paramfile
// and def_parm
// called form steady state file

model;
/////////////////////////////////////////////////////////////////
// 1.
-exp(c)^(-GAMMAC)*(1+2*PSI*(exp(k)/exp(k(-1)) -1)/exp(k(-1)))
+ BETA*exp(c(1))^(-GAMMAC)*((1-DELTAK)-2*PSI*(exp(k(1))/exp(k)-1)*
  (-exp(k(1))/exp(k)^2)+ALPHA*exp(a(1))*exp(k)^(ALPHA-1))= 
  -lambdak+BETA*(1-DELTAK)*elambdak;

// 2.
exp(c)+exp(k)-(1-DELTAK)*exp(k(-1))+
PSI*(exp(k)/exp(k(-1))-1)^2=exp(a)*exp(k(-1))^(ALPHA);

// 3.
exp(i) = exp(k)-(1-DELTAK)*exp(k(-1));

// 4.
lambdak = 0;

// 5. 
a = RHOA*a(-1)+erra;

// 6. 
exp(v) =  (exp(c)^(1-GAMMAC)-1)/(1-GAMMAC) + BETA*exp(v(1));

// 7
elambdak = lambdak(1);

// 8
klag = k(-1);

end;




shocks;
  var erra; stderr 0.015;
end;



stoch_simul(order=1,nocorr,nomoments,noprint,irf=0);
