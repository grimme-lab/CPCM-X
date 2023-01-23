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

module globals
   use mctc_env, only : wp
   use element_dict
   implicit none

   ! Global Parameters
   real(wp), parameter :: autokcal=627.509_wp, R=8.314_wp, jtokcal=0.000239006
   real(wp), parameter :: k_b=1.380649E-23,N_A=6.02214076E23, BtoA=0.529177
   real(wp), parameter :: pi = 4*atan(1.0_wp)
   real(wp), dimension(10) :: param, pr_param !Parameter for several methods
   real(wp) :: SysTemp !System Temperature

   real(wp) :: dG_is, dG_cc, dG_res, dG_disp, dG_shift !contributions to free energy of solvation

   logical :: ML
   ! Covalent Radii of Elements

   type(DICT_STRUCT), pointer :: cov_r

   ! PR Parameters

   type(DICT_STRUCT), pointer :: A_dsp, B_dsp, r_pr

   type :: configuration_type
      character(len=:), allocatable :: input
      character(len=:), allocatable :: smd_solvent
      character(len=:), allocatable :: csm_solvent
      character(len=:), allocatable :: csm_solute
      real(wp) :: T
      real(wp) :: probe
      real(wp) :: z1,z2
      real(wp) :: qc_eps
      character(len=:), allocatable :: sac_param_path
      character(len=:), allocatable :: smd_param_path
      character(len=:), allocatable :: database
      character(len=:), allocatable :: qc_calc
      character(len=:), allocatable :: xyz_input
      logical :: ML, sig_in, prof, smd_default, time, isodens
      character(len=:), allocatable :: model
      character(len=:), allocatable :: config_path
      !> Use internal parameters
      logical :: internal = .false.
   end type configuration_type

   contains
      ! Functions that are used in several modules
      function distance(xyz1,xyz2)
         implicit none
         real(wp), intent(in) :: xyz1(:), xyz2(:)
         real(wp) :: distance

         distance=(sqrt((xyz1(1)-xyz2(1))**2.0_wp+(xyz1(2)&
                  &-xyz2(2))**2.0_wp+(xyz1(3)-xyz2(3))**2.0_wp))

      end function distance

   Pure Function to_upper (str) Result (string)

   !   ==============================
   !   Changes a string to upper case
   !   ==============================

      Implicit None
      Character(*), Intent(In) :: str
      Character(LEN(str))      :: string

      Integer :: ic, i

      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

   !   Capitalize each letter if it is lowecase
       string = str
       do i = 1, LEN_TRIM(str)
         ic = INDEX(low, str(i:i))
           if (ic > 0) string(i:i) = cap(ic:ic)
       end do

   End Function to_upper

   Pure Function to_lower (str) Result (string)

   !   ==============================
   !   Changes a string to upper case
   !   ==============================

      Implicit None
      Character(*), Intent(In) :: str
      Character(LEN(str))      :: string

      Integer :: ic, i

      Character(26), Parameter :: low = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       Character(26), Parameter :: cap = 'abcdefghijklmnopqrstuvwxyz'

   !   Capitalize each letter if it is lowecase
       string = str
       do i = 1, LEN_TRIM(str)
         ic = INDEX(low, str(i:i))
           if (ic > 0) string(i:i) = cap(ic:ic)
       end do

   End Function to_lower

  subroutine rename(src, tgt, stat)
    use, intrinsic :: iso_c_binding, only : c_null_char
    character(len=*), intent(in) :: src
    character(len=*), intent(in) :: tgt
    integer, intent(out) :: stat
    interface
       function sys_rename(src, tgt, lena, lenb) bind(c, name="rename") result(stat)
          use, intrinsic :: iso_c_binding, only : c_char, c_int
          integer(c_int), intent(in) :: lena, lenb
          character(kind=c_char), intent(in) :: src(lena)
          character(kind=c_char), intent(in) :: tgt(lenb)
          integer(c_int) :: stat
       end function sys_rename
    end interface

 
    stat = sys_rename(src//c_null_char, tgt//c_null_char, len(src), len(tgt))
 end subroutine rename

end module globals
