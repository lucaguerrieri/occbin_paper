SUBROUTINE writevec(mat, rows, strow, endrow, filename)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, rows
DOUBLE PRECISION, DIMENSION(rows), INTENT(IN) :: mat
INTEGER :: i, j
CHARACTER(len=*), INTENT(IN) :: filename
DOUBLE PRECISION :: zerotol = 0.000000000000001

OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')


   DO i=strow,endrow
      IF (ABS(mat(i))>zerotol) THEN
         WRITE (4,100) i,mat(i)
      END IF
     END DO

100 FORMAT ('mat(', I4, ')=', ES30.15, ';')

CLOSE(UNIT=4)
END SUBROUTINE
