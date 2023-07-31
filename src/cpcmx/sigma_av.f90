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

module sigma_av
   use mctc_env, only : wp
   implicit none

   ! Module contains charge averaging and orthogonalizing algorythms
   !
   ! Averaging Uses:
   ! r_av : averaging radius
   ! xyz : Array of the coordinates of the segments
   ! charges : Array of Charge/Area of the Segments
   ! Area : Array of the Area of the Sements
   ! Results:
   ! av_charge : Averaged charge/area over the averaging radius r_av
   !
   ! Orthogonalizing Uses:
   ! v = Averaged charge/area over the averaging radius r_av
   ! v0 = Averaged charge/area over the averaging radius 2*r_av
   ! Results:
   ! vt = Orthogonalized charge/area



   contains

      subroutine average_one_charge (r_av, xyz, charges, area,av_charge)
         use globals
         implicit none
         real(wp), dimension(:), intent(in) :: charges, area
         real(wp), dimension(:,:), intent(in) :: xyz
         real(wp), intent(in) :: r_av
         real(wp), dimension(:), allocatable, intent(out) :: av_charge
         real(wp) :: tmpcharge, tmpcounter, tmpdenominator, r_u2, r_av2
         integer :: num, i, j

         num = size(charges)
         allocate(av_charge(num))
         r_av2=r_av**2.0_wp
         r_u2=0.0_wp
         tmpcharge=0.0_wp
         tmpcounter=0.0_wp
         tmpdenominator=0.0_wp
         do i=1,num
            do j=1,num
               r_u2=(area(j)/pi)
               tmpcounter=tmpcounter+(charges(j)*((r_u2*r_av2)/(r_u2+r_av2))*&
                         &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_wp)/(r_u2+r_av2))))
               tmpdenominator=tmpdenominator+(((r_u2*r_av2)/(r_u2+r_av2))*&
                             &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_wp)/(r_u2+r_av2))))
               r_u2=0.0_wp
            end do
            tmpcharge=tmpcounter/tmpdenominator
            tmpcounter=0.0_wp
            tmpdenominator=0.0_wp
            av_charge(i)=tmpcharge
            !write(*,*) av_charge(i)
         end do
      end subroutine average_one_charge


      subroutine ortho_charge (v,v0,vt)
         implicit none
         real(wp), dimension(:), allocatable, intent(in) :: v,v0
         real(wp), dimension(:), allocatable, intent(out) :: vt

         integer :: i

         allocate(vt(size(v)))

         do i=1,size(v)
            vt(i)=v0(i)-0.816_wp*v(i)
         end do

      end subroutine
end module sigma_av
