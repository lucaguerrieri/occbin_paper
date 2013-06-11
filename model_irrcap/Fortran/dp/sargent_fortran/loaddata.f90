SUBROUTINE loaddata(cof,rows,cols,filename)

IMPLICIT NONE

! define arguments
INTEGER, INTENT(IN) :: rows,cols
CHARACTER(len=200), INTENT(IN) :: filename
DOUBLE PRECISION, DIMENSION(rows,cols), INTENT(OUT) :: cof

! define working variables
DOUBLE PRECISION :: matsum

INTEGER :: skip, ierr
INTEGER :: indxi, indxj


OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

matsum = 0.

! now read the various pieces of the sparse matrix
DO indxj=1,cols
   DO indxi=1,rows
      READ (4,100) cof(indxi,indxj)
      matsum = matsum + cof(indxi,indxj)
   END DO
END DO

CLOSE (UNIT=4)

100 FORMAT (ES30.16)

END SUBROUTINE
