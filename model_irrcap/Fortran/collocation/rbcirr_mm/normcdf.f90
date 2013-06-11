DOUBLE PRECISION FUNCTION normcdf(x,mu,sigma)

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: x, mu, sigma


normcdf = (1.d0+erf((x-mu)/sqrt(2.d0*sigma**2)))/2.d0

END FUNCTION