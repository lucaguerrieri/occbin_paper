SUBROUTINE writevec_screen(vec, rows, strow, endrow)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, rows
DOUBLE PRECISION, DIMENSION(rows), INTENT(IN) :: vec
INTEGER :: i
DOUBLE PRECISION :: zerotol = 0.000000000000001
CHARACTER(len=500) :: output



DO i=strow,endrow
      IF (ABS(vec(i))>zerotol) THEN
         WRITE (*,100) i,vec(i)
         !when using matlab -- then
         !CALL mexPrintf(trim(output)//char(13))
      END IF
END DO


100 FORMAT ('vec(', I5, ')=',ES30.15,';')


END SUBROUTINE
