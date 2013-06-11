SUBROUTINE constructpol(polcoefs,poly,order,x)

INTEGER, INTENT(IN) :: order
DOUBLE PRECISION, INTENT(OUT) :: x
DOUBLE PRECISION, DIMENSION(order+1), INTENT(IN) :: polcoefs
DOUBLE PRECISION, DIMENSION(order+1), INTENT(IN) :: poly

INTEGER :: indxi


x =0.d0

DO indxi = 1,order+1
x = x + polcoefs(indxi)*poly(indxi)
END DO



END SUBROUTINE