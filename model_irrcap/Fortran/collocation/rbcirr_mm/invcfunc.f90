SUBROUTINE invcfunc(theta,kap,c,kap_prime)

USE model_parameters

DOUBLE PRECISION, INTENT(IN) :: theta, kap, c
DOUBLE PRECISION, INTENT(OUT) :: kap_prime

kap_prime = exp(theta)*kap**(ALPHA)+(1.d0-DELTAK)*kap-c

END SUBROUTINE