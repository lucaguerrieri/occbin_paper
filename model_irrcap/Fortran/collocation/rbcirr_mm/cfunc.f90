SUBROUTINE cfunc(theta,kap,kap_prime,c)

USE model_parameters

DOUBLE PRECISION, INTENT(IN) :: theta, kap, kap_prime
DOUBLE PRECISION, INTENT(OUT) :: c

c = exp(theta)*kap**(ALPHA)-kap_prime+(1.d0-DELTAK)*kap-PSI*(kap_prime/kap-1.d0)**2.d0

END SUBROUTINE