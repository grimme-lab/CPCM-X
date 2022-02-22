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
  private

  public :: unique, shell_sort
 
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

function unique(char_array) result(unique_array)
   !> Input Character Array
   character(len=*), allocatable, intent(in) :: char_array(:)
   !> Output unique array
   character(len=len(char_array)), allocatable:: unique_array(:)

   character(len=len(char_array)), allocatable :: tmp_array(:)

   character(len=len(char_array)) :: unique_element

   integer :: unique_size, i, j

   unique_size=1

   allocate(unique_array(unique_size))
   unique_array(1)=char_array(1)

   do i=1,size(char_array)
      unique_element=char_array(i)
      if (any(unique_element .eq. unique_array)) cycle
      allocate(tmp_array(unique_size))
      tmp_array=unique_array
      deallocate(unique_array)
      allocate(unique_array(unique_size+1))
      do j=1,unique_size
         unique_array(j)=tmp_array(j)
      end do
      deallocate(tmp_array)
      unique_size=unique_size+1
      unique_array(unique_size)=unique_element
   end do

   end function unique



 
END MODULE sort
  

