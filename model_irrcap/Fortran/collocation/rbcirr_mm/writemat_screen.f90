SUBROUTINE writemat_screen(mat, rows, cols, strow, endrow, stcol, endcol)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, stcol, endcol, rows, cols
DOUBLE PRECISION, DIMENSION(rows,cols), INTENT(IN) :: mat
INTEGER :: i, j
DOUBLE PRECISION :: zerotol = 0.000000000000001
CHARACTER(len=500) :: output


DO j=stcol,endcol
   DO i=strow,endrow
      IF (ABS(mat(i,j))>zerotol) THEN
         WRITE (*,100) i,j,mat(i,j)
         !when using matlab -- then
         !CALL mexPrintf(trim(output)//char(13))
      END IF
     END DO
END DO

100 FORMAT ('mat(', I4, ',', I4, ')=', ES30.15, ';')


END SUBROUTINE
