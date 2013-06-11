SUBROUTINE writemat_int(mat, rows, cols, strow, endrow, stcol, endcol, filename)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, stcol, endcol, rows, cols
INTEGER, DIMENSION(rows,cols), INTENT(IN) :: mat
INTEGER :: i, j
CHARACTER(len=200), INTENT(IN) :: filename

OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')

DO j=stcol,endcol
   DO i=strow,endrow
      
         WRITE (4,100) i,j,mat(i,j)
      
     END DO
END DO
100 FORMAT ('mat(', I7, ',', I7, ')=', I7, ';')

CLOSE(UNIT=4)
END SUBROUTINE
