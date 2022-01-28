! This file is part of CPCM-X.
! SPDX-Identifier: LGPL-3.0-or-later
!
! CPCM-X is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CPCM-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with CPCM-X.  If not, see <https://www.gnu.org/licenses/>.

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
  

