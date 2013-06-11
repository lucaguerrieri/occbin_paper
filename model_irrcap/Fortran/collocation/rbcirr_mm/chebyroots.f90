SUBROUTINE chebyroots(order,a,b,nodes)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: a,b
INTEGER, INTENT(IN) :: order
DOUBLE PRECISION, DIMENSION(order), INTENT(OUT) :: nodes


DOUBLE PRECISION :: pi = 3.141592653589793d0
INTEGER :: k

! see Judd, Algorithm 6.4 step 1



nodes=0.d0
WRITE(*,*) 'nodes'
DO k=1,order
    nodes(k) = cos(DBLE(2*k-1)/DBLE(2*order)*pi)
END DO

nodes = (nodes+1.d0)*(b-a)/2.d0+a

END SUBROUTINE