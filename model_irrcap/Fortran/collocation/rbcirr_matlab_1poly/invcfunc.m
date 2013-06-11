function kap_prime=invcfunc(theta,kap,c)

global ALPHA DELTAK 

 kap_prime = exp(theta).*kap.^(ALPHA)-c+(1-DELTAK)*kap;