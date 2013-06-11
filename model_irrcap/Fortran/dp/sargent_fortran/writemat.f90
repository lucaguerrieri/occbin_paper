SUBROUTINE writemat(mat, rows, cols, strow, endrow, stcol, endcol, filename)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, stcol, endcol, rows, cols
DOUBLE PRECISION, DIMENSION(rows,cols), INTENT(IN) :: mat
INTEGER :: i, j
CHARACTER(len=200), INTENT(IN) :: filename
DOUBLE PRECISION :: zerotol = 0.000000000000001

OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')

DO j=stcol,endcol
   DO i=strow,endrow
      IF (ABS(mat(i,j))>zerotol) THEN
         WRITE (4,100) i,j,mat(i,j)
      END IF
     END DO
END DO
100 FORMAT ('mat(', I4, ',', I4, ')=', ES30.15, ';')

CLOSE(UNIT=4)
END SUBROUTINE
