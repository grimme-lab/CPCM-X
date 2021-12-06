! This file is part of COSMO-X.
! SPDX-Identifier: LGPL-3.0-or-later
!
! COSMO-X is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! COSMO-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with COSMO-X.  If not, see <https://www.gnu.org/licenses/>.

module bonding
   use mctc_env, only : wp
   implicit none

   ! This module is used to determine the Bonding Situation of the Compounds

   contains

   ! This Subroutine determines covalent bonds between Atoms/Segments
   ! It uses the Positions of the Segments and the hard coded covalent Radii
   ! A Covalent Bond is assumed when the distance is lower than the sum of the convalent Radii*1.15
   subroutine det_bonds(ident,xyz,elements,is_bonded,oh_count,nh_count)
      use globals
      use element_dict
      implicit none

      integer, dimension(:), allocatable, intent(in) :: ident
      character(2), dimension(:), allocatable, intent(in) :: elements
      real(wp), dimension(:,:), allocatable, intent(in) :: xyz

      logical, dimension(:,:), allocatable, intent(out) :: is_bonded

      integer, intent(out), optional :: oh_count, nh_count

      type(DICT_DATA) :: radius
      integer :: i,j
      real(wp) :: cov_a, cov_b

      allocate(is_bonded(size(ident),size(ident)))
      is_bonded(:,:) = .FALSE.

      cov_a=0.0_wp
      cov_b=0.0_wp
      do i=1,int(maxval(ident))
         radius=dict_get_key(cov_r,elements(i))
         cov_a=radius%param
         do j=i+1,int(maxval(ident))
            radius=dict_get_key(cov_r,elements(i))
           ! write(*,*) cov_a, cov_b
            cov_b=radius%param
            if (distance(xyz(i,:),xyz(j,:)) .LE. ((cov_a+cov_b)*1.15)) then
               is_bonded(i,j)=.TRUE.
               is_bonded(j,i)=.TRUE.
            end if
            cov_b=0.0_wp
         end do
         cov_a=0.0_wp
      end do

      if (present(nh_count)) then
         oh_count=0
         nh_count=0
         do i=1,int(maxval(ident))
            select case (elements(i))
               case default
                  cycle
               case ("o")
                  do j=1,int(maxval(ident))
                     if ((is_bonded(i,j)) .AND. (elements(j) .eq. "h")) then
                        oh_count=oh_count+1
                        exit
                     end if
                  end do
               case ("n")
                  do j=1,int(maxval(ident))
                     if ((is_bonded(i,j)) .AND. (elements(j) .eq. "h")) then
                        nh_count=nh_count+1
                        exit
                     end if
                  end do
            end select
         end do
        ! write(*,*) oh_count, nh_count
      end if


      ! Check points, uncomment the following:
     ! do i=1,int(maxval(ident))
     !    write(*,*) i, elements(i)
     !    do j=1,int(maxval(ident))
     !       write(*,*) "to ",j,is_bonded(i,j)
     !    end do
     ! end do

   end subroutine

   subroutine hb_grouping(ident,elements,is_bonded,hb_group)
      implicit none

      integer, dimension(:), allocatable, intent(in) :: ident
      character(2), dimension(:), allocatable, intent(in) :: elements
      logical, dimension(:,:), allocatable, intent(in) :: is_bonded

      character(2), dimension(:), allocatable, intent(out) :: hb_group


      integer :: i,j


      allocate(hb_group(size(ident)))

      hb_group(:) = 'NH'

      do i=1,size(ident)
         if (elements(int(ident(i))) .EQ. 'h') then
            do j=1,size(ident)
               if (is_bonded(int(ident(i)),int(ident(j)))) then
                  if (elements(int(ident(j))) .EQ. 'o') then
                     hb_group(i)='OH'
                     hb_group(j)='OH'
                     exit
                  else if ((elements(int(ident(j))) .EQ. 'n') .OR. &
                     &(elements(int(ident(j))) .EQ. 'f')) then
                     hb_group(i)='OT'
                     hb_group(j)='OT'
                     exit
                  end if
               end if
            end do
         end if

         if (elements(int(ident(i))) .EQ. 'o') then
            do j=1,size(ident)
               if (is_bonded(int(ident(i)),int(ident(j)))) then
                  if (elements(int(ident(j))) .EQ. 'h') then
                     hb_group(i)='OH'
                     hb_group(j)='OH'
                     exit
                  end if
               end if
            end do
            if (hb_group(i) .EQ. 'NH') then
               hb_group(i) = 'OT'
            end if
         end if

         if (elements(int(ident(i))) .EQ. 'n') hb_group(i)='OT'
         if (elements(int(ident(i))) .EQ. 'f') hb_group(i)='OT'

       ! write(*,*) i, elements(int(ident(i))), hb_group(i)
      end do

!      write(*,*) hb_group(66),size(ident)
!      do i=1,size(ident)
!         write(*,*) i,elements(int(ident(i))), is_bonded(i,66), is_bonded(66,i)
!      end do

   end subroutine hb_grouping

   subroutine det_rings(ident,is_bonded,is_ring,N_ear)
      implicit none

      integer, dimension(:), intent(in) :: ident
      logical, dimension(:,:), intent(in) :: is_bonded
      logical, dimension(:), intent(out), allocatable :: is_ring
      integer, intent(out) :: N_ear ! Number of effective atoms in rings

      logical, dimension(:), allocatable :: visited

      integer :: eles, i, j, cut

      eles=int(maxval(ident))

      allocate(is_ring(eles))
      allocate(visited(eles))

      do i=1,eles
         visited(:)=.false.
         is_ring(i)=.false.
         cut=1
         do j=1,eles
            if (is_bonded(i,j)) then
               Call iter_ring(i,i,j,is_ring(i),eles,is_bonded,cut,visited)
            end if
            if (is_ring(i)) exit
         end do
         !write(*,*) i, is_ring(i)
         !stop
      end do

      N_ear=0

      do i=1,eles
         if (is_ring(i)) then
            do j=1,eles
               if (is_bonded(i,j)) then
                  if (.NOT. is_ring(j)) then
                     N_ear=N_ear+1
                     exit
                  end if
               end if
               if (j .eq. eles) N_ear=N_ear+2
            end do
         end if
      end do



   end subroutine det_rings

   recursive subroutine iter_ring(check,i,j,is_ring,eles,is_bonded,cut,visited)

      integer, intent(in) :: check, i,j,eles
      integer, intent(in) :: cut
      logical, intent(inout) :: is_ring
      logical, dimension(:,:), intent(in) :: is_bonded
      logical, dimension(:) :: visited

      integer :: z

      visited(j)=.true.

      do z=1,eles
         if ((is_bonded(j,z)) .AND. (i .ne. z) .AND. (.NOT. visited(z))) then
            !write(*,*) check,z, cut
            if (check .eq. z) then
               is_ring=.true.
               exit
            else if (is_ring) then
               exit
            else

               if (cut .GE. 11) exit
               Call iter_ring(check,j,z,is_ring,eles,is_bonded,cut+1,visited)
            end if
         end if
      end do

   end subroutine iter_ring





end module bonding
