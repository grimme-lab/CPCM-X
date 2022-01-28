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

module pr
   use mctc_env, only : wp
   implicit none

   !Peng-Robinson Equation of State model for disperive interactions

   contains



      subroutine pr2018(area,elements,ident,oh_count,nh_count,n_ear)
         use globals
         use element_dict
         implicit none
         !New version of PR EOS
         real(wp), dimension(:), allocatable, intent(in) :: area
         integer, allocatable, intent(in) :: ident(:)
         character(2), dimension(:), allocatable, intent(in) :: elements
         integer, intent(in) :: oh_count, nh_count, n_ear

         character(2) :: element
         real(wp), dimension(:), allocatable :: atom_area
         real(wp) :: dG_hb, dG_ring, dG_vdw
         integer :: i
         type(DICT_DATA) :: A, B

         !write(*,*) n_ear
         dG_hb = 0.0_wp
         dG_ring = 0.0_wp
         dG_vdw = 0.0_wp
         allocate(atom_area(int(maxval(ident))))
         atom_area=0.0_wp
         do i=1,size(area)
            atom_area(int(ident(i)))=atom_area(int(ident(i)))+area(i)
         end do

        ! do i=1,int(maxval(ident))
        !    write(*,*) i, atom_area(i)
        ! end do

        do i=1,int(maxval(ident))
         element=elements(i)
         A=dict_get_key(A_dsp,element)
         B=dict_get_key(B_dsp,element)
        ! write(*,*) atom_area(i), A%param, B%param
         dG_vdw=dG_vdw+atom_area(i)*(A%param*log(SysTemp)+B%param)
        end do
       !  write(*,*) oh_count
         dG_hb=oh_count*(pr_param(3)*(log(SysTemp)/SysTemp)+pr_param(4))&
            &+nh_count*(pr_param(5)*(log(SysTemp)/SysTemp)+pr_param(6))
!write(*,*) dG_hb
         dG_ring=n_ear*(pr_param(1)*log(SysTemp)+pr_param(2))

         dG_vdw=dG_vdw*jtokcal

         ! dG_disp = dG_vdw + dG_hb + dG_ring
        ! write(*,*)  dG_vdw, dG_hb, dG_ring, dG_disp
        if (ML) then
           open(5,file='ML.pr')
           write(5,'(F0.2,A,I0,A,I0,A,I0,A,I0)') SysTemp,",",oh_count,",",nh_count,",",n_ear,",",int(maxval(ident))
           close(5)
         end if


      end subroutine pr2018








end module pr
