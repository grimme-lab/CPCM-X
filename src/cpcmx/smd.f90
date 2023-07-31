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

module sdm
   implicit none
   private

   public :: calculate_cds

   interface calculate_cds
      module procedure :: calculate_cds_normal
      module procedure :: calculate_cds_isodens
      module procedure :: calculate_cds_internal
   end interface calculate_cds

   contains

      !> Calculate the CDS energy for a molecule input with internal SMD parameters
   subroutine calculate_cds_internal(species, symbols, coord, probe, solvent, dG_cds,internal_smd)
      use mctc_env, only : wp
      use numsa, only : surface_integrator, new_surface_integrator, get_vdw_rad_smd, grid_size &
         & , get_vdw_rad_cosmo
      use smd, only: init_smd, smd_param, calc_surft, smd_surft, calc_cds, ascii_cds
      use globals, only: dG_disp, BtoA, autokcal
      use iso_fortran_env, only: output_unit
      !> Unique chemical species in the input structure, shape: [nat]
      integer, intent(in) :: species(:)
      !> Element symbol for each chemical species, shape: [nsp]
      character(len=*), intent(in) :: symbols(:)
      !> Cartesian coordinates in Bohr, shape: [3, nat]
      real(wp), intent(in) :: coord(:, :)
      real(wp), allocatable :: coord_rev(:,:)
      !> Probe radius for surface area integration in Bohr
      real(wp), intent(in) :: probe
      !> CDS Energy
      real(wp), intent(out) :: dG_cds
      !> Accessible surface area in Bohr², shape: [nat]
      real(wp),allocatable :: surface(:)
      !> Derivative of surface area w.r.t. atomic displacements, shape: [3, nat, nat]
      real(wp),allocatable :: dsdr(:, :, :)
      !> Using default SMD Parameters?
      character(len=*), intent(in), dimension(:) :: internal_smd
      !>Laufen
      integer :: i, j
      !> Read Env
      integer :: dummy1,io_error
      !> Parameter Path and Solvent Name
      character(len=*) :: solvent
      character(len=:), allocatable :: path
      logical :: ex
      real(wp),allocatable :: cds(:)
      real(wp) :: cds_sm

      type(surface_integrator) :: sasa
      type(smd_param) :: param
      type(smd_surft) :: surft
      real(wp), allocatable :: rad(:)

      allocate (surface(size(species)))
      allocate (dsdr(3,size(species),size(species)))
      allocate(coord_rev(3,size(species)))

      select case (solvent)
      case ('water','Water','WATER','h2o','H2O','H2o')
         path = 'smd_h2o'
      case default
         path = 'smd_ot'
      end select

      !> Writing SMD Parameters to File
      open(output_unit, file=path,status='unknown')
      do i=1,size(internal_smd)
         write(output_unit,'(a)') trim(internal_smd(i))
      end do
      close(output_unit)

      do i=1,3
         do j=1,size(species)
            coord_rev(i,j)=coord(j,i)*(1/BtoA)
         end do
      end do

      rad = get_vdw_rad_smd(symbols)
      Call init_smd(param,solvent)

      if (param%alpha .lt. 0.43) then
        rad(8)=rad(8)+1.8*(0.43-param%alpha)
      end if

      call new_surface_integrator(sasa, species, rad, probe, grid_size(8))
      call sasa%get_surface(species, coord_rev, surface, dsdr)
      Call calc_surft(coord_rev,species,symbols,param,surft)
      Call calc_cds(surft,surface,cds,cds_sm)
      dG_disp= (sum(cds)+cds_sm)/1000
      dG_cds=dG_disp/autokcal

   end subroutine calculate_cds_internal

      !> Example implementation to calculate surface area for a molecule input
      subroutine calculate_cds_normal(species, symbols, coord, probe, solvent, path, dG_cds,default)
      use mctc_env, only : wp
      use numsa, only : surface_integrator, new_surface_integrator, get_vdw_rad_smd, grid_size &
         & , get_vdw_rad_cosmo
      use smd, only: init_smd, smd_param, calc_surft, smd_surft, calc_cds, ascii_cds
      use globals, only: dG_disp, BtoA, autokcal
      !> Unique chemical species in the input structure, shape: [nat]
      integer, intent(in) :: species(:)
      !> Element symbol for each chemical species, shape: [nsp]
      character(len=*), intent(in) :: symbols(:)
      !> Cartesian coordinates in Bohr, shape: [3, nat]
      real(wp), intent(in) :: coord(:, :)
      real(wp), allocatable :: coord_rev(:,:)
      !> Probe radius for surface area integration in Bohr
      real(wp), intent(in) :: probe
      !> CDS Energy
      real(wp), intent(out) :: dG_cds
      !> Accessible surface area in Bohr², shape: [nat]
      real(wp),allocatable :: surface(:)
      !> Derivative of surface area w.r.t. atomic displacements, shape: [3, nat, nat]
      real(wp),allocatable :: dsdr(:, :, :)
      !> Using default SMD Parameters?
      logical, intent(in) :: default
      !>Laufen
      integer :: i, j
      !> Read Env
      integer :: dummy1,io_error
      !> Parameter Path and Solvent Name
      character(len=*) :: path, solvent
      logical :: ex
      real(wp),allocatable :: cds(:)
      real(wp) :: cds_sm

      type(surface_integrator) :: sasa
      type(smd_param) :: param
      type(smd_surft) :: surft
      real(wp), allocatable :: rad(:)

      allocate (surface(size(species)))
      allocate (dsdr(3,size(species),size(species)))
      allocate(coord_rev(3,size(species)))


      if (.NOT. default) then
         select case (solvent)
            case ('h2o','water')
               INQUIRE(file=path//"smd_h2o",exist=ex)
               if (.not. ex) INQUIRE(file=path,exist=ex)
            case default
               INQUIRE(file=path//"smd_ot",exist=ex)
               if (.not. ex) INQUIRE(file=path,exist=ex)
         end select
         if (.NOT. ex) then
            write(*,*) "Parameter file for SMD model does not exists in the specified path."
            write(*,*) "Path: "//path
            write(*,*) "You can skip this check (using default parameters) by using the default_smd keyword!"
            stop
         end if
      end if
      do i=1,3
         do j=1,size(species)
            coord_rev(i,j)=coord(j,i)*(1/BtoA)
         end do
      end do

      rad = get_vdw_rad_smd(symbols)
      if (default) then
         Call init_smd(param,solvent)
      else
         Call init_smd(param,solvent,path)
      end if

      if (param%alpha .lt. 0.43) then
        rad(8)=rad(8)+1.8*(0.43-param%alpha)
      end if

      call new_surface_integrator(sasa, species, rad, probe, grid_size(8))
      call sasa%get_surface(species, coord_rev, surface, dsdr)
      Call calc_surft(coord_rev,species,symbols,param,surft)
      Call calc_cds(surft,surface,cds,cds_sm)
      dG_disp= (sum(cds)+cds_sm)/1000
      dG_cds=dG_disp/autokcal
end subroutine calculate_cds_normal

subroutine calculate_cds_isodens(species,symbols,coord,probe,solvent,path,dG_cds,rad)
  ! subroutine calculate_cds_normal(species, symbols, coord, probe, solvent, path, default)
      use mctc_env, only : wp
      use numsa, only : surface_integrator, new_surface_integrator, grid_size 
      use smd, only: init_smd, smd_param, calc_surft, smd_surft, calc_cds
      use globals, only: dG_disp, BtoA,autokcal
      !> Unique chemical species in the input structure, shape: [nat]
      integer, intent(in) :: species(:)
      !> Isodens Radii calculated by CPCM-X for by species, shape: [nat]
      real(wp), intent(in) :: rad(:)
      !> Element symbol for each chemical species, shape: [nsp]
      character(len=*), intent(in) :: symbols(:)
      !> Cartesian coordinates in Bohr, shape: [3, nat]
      real(wp), intent(in) :: coord(:, :)
      real(wp), allocatable :: coord_rev(:,:)
      !> Probe radius for surface area integration in Bohr
      real(wp), intent(in) :: probe
      !> CDS energy
      real(wp), intent(out) :: dG_cds
      !> Accessible surface area in Bohr², shape: [nat]
      real(wp),allocatable :: surface(:)
      !> Derivative of surface area w.r.t. atomic displacements, shape: [3, nat, nat]
      real(wp),allocatable :: dsdr(:, :, :)
      !>Laufen
      integer :: i, j
      !> Read Env
      integer :: dummy1,io_error
      !> Parameter Path and Solvent Name
      character(len=*) :: path, solvent
      logical :: ex
      real(wp),allocatable :: cds(:)
      real(wp) :: cds_sm

      type(surface_integrator) :: sasa
      type(smd_param) :: param
      type(smd_surft) :: surft

      allocate (surface(size(species)))
      allocate (dsdr(3,size(species),size(species)))
      allocate(coord_rev(3,size(species)))



      select case (solvent)
         case ('h2o','water')
            INQUIRE(file=path//"smd_h2o",exist=ex)
            if (.not. ex) INQUIRE(file=path,exist=ex)
         case default
            INQUIRE(file=path//"smd_ot",exist=ex)
            if (.not. ex) INQUIRE(file=path,exist=ex)
      end select
      if (.NOT. ex) then
         write(*,*) "Parameter file for SMD model does not exists in the specified path."
         write(*,*) "Path: "//path
         write(*,*) "You can skip this check (using default parameters) by using the default_smd keyword!"
         stop
      end if

      do i=1,3
         do j=1,size(species)
            coord_rev(i,j)=coord(j,i)*(1/BtoA)
         end do
      end do

      Call init_smd(param,solvent,path)

      ! if (param%alpha .lt. 0.43) then
      !   rad(8)=rad(8)+1.8*(0.43-param%alpha)
      ! end if

      call new_surface_integrator(sasa, species, rad, probe, grid_size(8))
      call sasa%get_surface(species, coord_rev, surface, dsdr)
      Call calc_surft(coord_rev,species,symbols,param,surft)
      Call calc_cds(surft,surface,cds,cds_sm)
      dG_disp= (sum(cds)+cds_sm)/1000
      dG_cds=dG_disp/autokcal
end subroutine calculate_cds_isodens




end module sdm
