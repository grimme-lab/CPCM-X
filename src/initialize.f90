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

module initialize_cosmo
   use mctc_env, only : wp, error_type, fatal_error
   use, intrinsic :: iso_fortran_env, only : output_unit
   implicit none

contains
   subroutine read_cosmo(compound,elements,ident,xyz,charges,area,pot,volume,c_energy,atom_xyz,database)
      use globals
      character(*), intent(in) :: compound
      character(*), intent(in) :: database
      character(100) :: line, ld1,ld2, ld3, ld4, ld5, ld6,home
      character(:), allocatable :: filen
      real(wp), dimension(:), allocatable,intent(out) :: charges,&
         &area,pot
      real(wp), dimension(:,:), allocatable, intent(out) :: xyz, atom_xyz
      character(2), allocatable, dimension(:),intent(out) :: elements
      real(wp), intent(out) :: volume, c_energy
      integer :: io_error, dummy1, num, ele_num
      real(wp) :: dummy3, dummy4, dummy5
      integer, allocatable :: ident(:)
      character(2) :: element

      logical :: exists

      filen=trim(compound)

      INQUIRE(file=trim(filen),EXIST=exists)

      if (.NOT. exists) then
         if (database .ne. "NONE") then
            filen=trim(database)//"/"//compound
            INQUIRE(file=trim(filen), EXIST=exists)
            write(output_unit,'(5x,a)') "Database specified, reading .cosmo file from", &
               filen
         end if
      end if
      if (.NOT. exists) then
         Call get_environment_variable("CSMHOME", home,dummy1,io_error,.TRUE.)
         if (io_error .EQ. 0) then
            filen=trim(home)//"/"//compound
            INQUIRE(file=trim(filen), EXIST=exists)
               if (exists) then
                  write(*,*) "No COSMO file in working directory, reading COSMO file from ", filen
               else
                  write(*,*) "No COSMO file in working directory or HOME directory for ", compound
                  stop
               end if
         else
            write(*,*) "No COSMO file in working directory for", compound,"No Home directory to check."
            stop
         end if
      end if

      open(1,file=trim(filen),iostat=io_error)

      if (io_error .NE. 0) then
         write(*,*) "Problem while reading COSMO input. &
            &Check your data."
         stop
      end if

      io_error=0
      do while (io_error .GE. 0)
         read(1,*,iostat=io_error) line
      end do
      read(line,*) num
      allocate(charges(num))
      allocate(ident(num))
      allocate(area(num))
      allocate(xyz(num,3))
      allocate(pot(num))
      rewind(1)
      ident(:)=0
      do while (line .NE. "$segment_information")
         read(1,*) line
      end do

      num=1
      dummy4=0
      do while (.TRUE.)
         read(1,'(A1)',advance='no',iostat=io_error) line
         if (line=="#") then
            read(1,*)
            cycle
         else if (io_error .LT. 0) then
            exit
         else
            read(1,*) dummy1,ident(num),xyz(num,1),xyz(num,2),&
               &xyz(num,3),dummy3,area(num),charges(num),pot(num)
         !   charges(num)=anint(charges(num)*1000)/1000
            pot(num)=pot(num)*BtoA
            xyz(num,1)=xyz(num,1)*btoa
            xyz(num,2)=xyz(num,2)*btoa
            xyz(num,3)=xyz(num,3)*btoa
            num=num+1
         end if
      end do
      rewind(1)
      do while (line .NE. "#atom")
         read(1,*) line
      end do

      ele_num=0
      do while (.TRUE.)
         read(1,'(A1)',iostat=io_error) line
         if (index(line,"$") .eq. 0) then
           ele_num=ele_num+1
         else
            exit
         end if
      end do

      allocate(elements(ele_num))
      allocate(atom_xyz(ele_num,3))

      rewind(1)
      do while (line .NE. "#atom")
         read(1,*) line
      end do

      do while (.TRUE.)
         read(1,'(A1)',advance='no',iostat=io_error) line
         if (line .NE. "$") then
            read(1,*) dummy1, dummy3, dummy4, dummy5, element
            elements(dummy1)=to_lower(element)
            !write(*,*) elements(dummy1)
            atom_xyz(dummy1,1)=dummy3*btoa
            atom_xyz(dummy1,2)=dummy4*btoa
            atom_xyz(dummy1,3)=dummy5*btoa
         else
            exit
         end if
      end do
!      write(*,*) elements
      rewind(1)
      do while (line .NE. "area=")
         read(1,*) line
      end do

      read(1,*) line, volume

      INQUIRE(file=filen(:len(filen)-6)//".energy",exist=exists)
      if (exists) then
         write(output_unit,'(5x,a)') "Reading energy for compound "//trim(compound)&
            &//" from "//filen(:len(filen)-6)//".energy"
         open(11,file=filen(:len(filen)-6)//".energy")
         read(11,*) c_energy
         close(11)
      else
         rewind(1)
         do while (line .NE. "$cosmo_energy")
            read(1,*) line
         end do
         read (1,*)
         read (1,*) line,ld1,ld2,ld3,ld4,ld5,ld6,c_energy
      end if

      close(1)

   end subroutine read_cosmo

   subroutine initialize_param(filename,model,r_cav,disp_con, solvent, error)
      use element_dict
      use, intrinsic :: iso_fortran_env, only: output_unit
      use globals, only: param, cov_r, dG_shift

      type(DICT_STRUCT), pointer, intent(inout) :: r_cav, disp_con

      type(DICT_DATA) :: data1, r_c, d_c
      character(len=3) :: symbol
      character(len=*) :: filename, model, solvent
      logical :: g_exists
      integer :: i, io_error,dummy1
      character(len=100) :: home,param_path

      !> Error Handling
      type(error_type), allocatable :: error

      INQUIRE(file=filename, exist=g_exists)
      if (.NOT. g_exists) then
         Call fatal_error(error,"No Parameter File for CPCM-X found.") 
         return
      else
         open(1,file=filename)

         select case (trim(model))
            case("crs")

               ! Setting global COSMO-RS Parameters from parameter file

               read(1,*) param(1)
               param(2)=2.0_wp*param(1)
               do i=3,10
                  read(1,*) param(i)
               end do
               read(1,*) dG_shift

            case("sac")

               ! Setting global COSMO-SAC Parameters from parameter file

               read(1,*) param(1)
               io_error=0
               do i=2,8
                  read(1,*,iostat=io_error) param(i)
               end do
               read(1,*) dG_shift

            case("sac2010")

               do i=1,10
                  read(1,*) param(i)
               end do
               read(1,*) dG_shift

            case("sac2013")

               do i=1,10
                  read(1,*) param(i)
               end do
               read(1,*)
               io_error=0

               read(1,*) symbol, d_c%param
               Call dict_create(disp_con, trim(symbol), d_c)
               do while (io_error .GE. 0)
                  read(1,*,iostat=io_error) symbol, d_c%param
                  Call dict_add_key(disp_con,trim(symbol), d_c)
               end do

            case default
               write(*,*) model//" is not a defined model!"
               stop
         end select

      end if

      ! Negative makes no sense
      param(1)=abs(param(1))

         !Hard Coded Covalent Radii

         data1%param=0.31_wp
         call dict_create(cov_r, 'h', data1)
         data1%param=0.28_wp
         call dict_add_key(cov_r, 'he', data1)
         data1%param=1.28_wp
         call dict_add_key(cov_r, 'li', data1)
         data1%param=0.96_wp
         call dict_add_key(cov_r, 'be', data1)
         data1%param=0.84_wp
         call dict_add_key(cov_r, 'b', data1)
         data1%param=0.76_wp
         call dict_add_key(cov_r, 'c', data1)
         data1%param=0.71_wp
         call dict_add_key(cov_r, 'n', data1)
         data1%param=0.66_wp
         call dict_add_key(cov_r, 'o', data1)
         data1%param=0.57_wp
         call dict_add_key(cov_r, 'f', data1)
         data1%param=0.58_wp
         call dict_add_key(cov_r, 'ne', data1)
         data1%param=1.66_wp
         call dict_add_key(cov_r, 'na', data1)
         data1%param=1.41_wp
         call dict_add_key(cov_r, 'mg', data1)
         data1%param=1.21_wp
         call dict_add_key(cov_r, 'al', data1)
         data1%param=1.11_wp
         call dict_add_key(cov_r, 'si', data1)
         data1%param=1.07_wp
         call dict_add_key(cov_r, 'p', data1)
         data1%param=1.05_wp
         call dict_add_key(cov_r, 's', data1)
         data1%param=1.02_wp
         call dict_add_key(cov_r, 'cl', data1)
         data1%param=1.06_wp
         call dict_add_key(cov_r, 'ar', data1)
         data1%param=2.03_wp
         call dict_add_key(cov_r, 'k', data1)
         data1%param=1.76_wp
         call dict_add_key(cov_r, 'ca', data1)
         data1%param=1.70_wp
         call dict_add_key(cov_r, 'sc', data1)
         data1%param=1.60_wp
         call dict_add_key(cov_r, 'ti', data1)

   end subroutine initialize_param

   subroutine init_pr
      use globals
      use element_dict
      !Initialize PR Parameters

      character(2) :: symbol
      type(DICT_DATA) :: r_i,A,B
      integer:: i,io_error, dummy1
      character(len=100) :: param_path, home

      Call get_environment_variable("CSMHOME", home,dummy1,io_error,.TRUE.)
      if (io_error .EQ. 0) then
         param_path=trim(home)//"/pr.param"
      else if (io_error .EQ. 1) then
         param_path="pr.param"
      end if

      open(1,file=param_path)
      do i=1,6
         read(1,*) pr_param(i)
      !   write(*,*) pr_param(i)
      end do
      read(1,*)
      read(1,*) symbol, r_i%param, A%param, B%param

      !write(*,*) symbol, r_i%param, A%param, B%param
      Call dict_create(r_pr, trim(symbol), r_i)
      Call dict_create(A_dsp, trim(symbol), A)
      Call dict_create(B_dsp, trim(symbol), B)
      do while (io_error .GE. 0)
         read(1,*,iostat=io_error) symbol, r_i%param, A%param, B%param
         Call dict_add_key(r_pr, trim(symbol), r_i)
         Call dict_add_key(A_dsp, trim(symbol), A)
         Call dict_add_key(B_dsp, trim(symbol), B)
      end do
      close(1)
   end subroutine init_pr

end module initialize_cosmo
