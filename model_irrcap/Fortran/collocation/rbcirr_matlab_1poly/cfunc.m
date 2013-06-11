function c=cfunc(theta,kap,kap_prime)

global ALPHA DELTAK PSI

c = exp(theta).*kap.^(ALPHA)-kap_prime+(1-DELTAK)*kap...
            -PSI*(kap_prime./kap-1).^2;