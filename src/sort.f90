MODULE sort
  use mctc_env, only : wp
  implicit none
 
CONTAINS
 
SUBROUTINE Shell_Sort(a)
 
  IMPLICIT NONE
  INTEGER :: i, j, increment
  REAL(wp) :: temp
  REAL(wp), INTENT(in out) :: a(:)
 
  increment = SIZE(a) / 2
  DO WHILE (increment > 0)
      DO i = increment+1, SIZE(a)
         j = i
         temp = a(i)
         DO WHILE (j >= increment+1 .AND. a(j-increment) > temp)
            a(j) = a(j-increment)
            j = j - increment
         END DO
         a(j) = temp
      END DO
      IF (increment == 2) THEN
         increment = 1
      ELSE
         increment = increment * 5 / 11
      END IF      
  END DO
 
END SUBROUTINE Shell_Sort
 
END MODULE sort
  

