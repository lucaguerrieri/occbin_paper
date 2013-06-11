MODULE wrapper_vars

IMPLICIT NONE
SAVE

INTEGER :: nstates, order, order_mm
DOUBLE PRECISION :: kstart, kend
DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: P
DOUBLE PRECISION, POINTER, DIMENSION(:) :: theta, nodes, nodes_mm



END MODULE wrapper_vars
