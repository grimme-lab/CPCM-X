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

module isodens
   use mctc_env, only : wp
   use globals, only : distance

   PRIVATE

   public:: get_isodens_radii

contains
    subroutine get_isodens_radii(segment_xyz,segment_ident,atom_xyz, isodens_radii)
        !> Coordinates of Segments from COSMO File
        real(wp), allocatable, intent(in) :: segment_xyz(:,:)
        !> Identifier of Segments from COSMO File
        integer, allocatable, intent(in) :: segment_ident(:)
        !> Coordinates of Atoms from COSMO File
        real(wp), allocatable, intent(in) :: atom_xyz(:,:)

        !> Calulcate isodens radii from average distance of the surface segments
        real(wp), intent(out), allocatable :: isodens_radii(:)

        integer :: natom

        integer :: i, j, divider

        natom=maxval(segment_ident)
        allocate(isodens_radii(maxval(segment_ident)))

        isodens_radii=0.0_wp
        divider=0

        do i=1,natom
            do j=1,size(segment_ident)
                if (segment_ident(j) .eq. i) then
                    isodens_radii(i)=&
                    &isodens_radii(i)+distance(segment_xyz(j,:),atom_xyz(i,:))
                    divider=divider+1
                end if
            end do
            if (divider == 0) then
                isodens_radii(i)=0
            else
                isodens_radii(i)=isodens_radii(i)/divider
            end if
            divider=0
        end do

    end subroutine get_isodens_radii

end module isodens