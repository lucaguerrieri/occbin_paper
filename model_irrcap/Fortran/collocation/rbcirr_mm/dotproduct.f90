SUBROUTINE dotproduct(npoints,x,y,product)

INTEGER, INTENT(IN) :: npoints
DOUBLE PRECISION, DIMENSION(npoints), INTENT(IN) :: x, y
DOUBLE PRECISION, INTENT(OUT) :: product

INTEGER :: i

product = 0.d0
DO i=1,npoints
	product = product + x(i)*y(i)
END DO


END SUBROUTINE