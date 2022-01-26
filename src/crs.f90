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

module crs
   use mctc_env, only : wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use crs_broyden, only : broyden_mixer, new_broyden
   implicit none

contains


subroutine calcgas(E_cosmo,id_scr,gas_chem,area,sv,su,pot,element,ident,disp_con, T,r_cav)
   use globals
   use element_dict
   real(wp), intent(out) :: id_scr, gas_chem
   real(wp), intent(in) :: T, E_cosmo
   real(wp),dimension(:),allocatable, intent(in) :: area, sv, su, pot
   integer, allocatable, intent(in) :: ident(:)
   character(2), dimension(:), allocatable, intent(in) :: element
   type(DICT_STRUCT), pointer, intent(in) :: disp_con, r_cav
   !real(wp), dimension(10) :: param
   type(DICT_DATA) :: disp!, r_c
   real(wp) :: E_gas, dEreal, ediel, edielprime, vdW_gain, thermo, beta, avcorr
   integer :: dummy1, ioerror, i

   !> Check if gas phase energy exists.
   logical :: ex

   INQUIRE(file="gas.energy", exist=ex)
   if (.not. ex) error stop "No gas.energy file found. Use TM keyword or manually set up the gas phase energy."
   open(1,file="gas.energy")
   read(1,*,iostat=ioerror) E_gas
   if (ioerror .NE. 0) error stop "Problem while reading energies (check gas.energy file)."
   dEreal=(E_cosmo-E_gas)
   ediel=0
   edielprime=0
   do i=1,size(sv)
      ediel=ediel+(area(i)*pot(i)*su(i))
      edielprime=edielprime+(area(i)*pot(i)*sv(i))
   end do
   avcorr=(edielprime-ediel)/2.0_wp*0.8_wp
   write(*,*) "E_COSMO+dE: ", (E_cosmo+avcorr)*autokcal
   write(*,*) "E_gas: ", E_gas*autokcal
   dEreal=dEreal*autokcal
   id_scr=dEreal+avcorr*autokcal
   write(*,*) "E_COSMO-E_gas+dE: ", (E_cosmo-E_gas+avcorr)*autokcal
   write(*,*) "Ediel: ", ediel/2*autokcal
   write(*,*) "Averaging corr dE: ", avcorr*autokcal

   vdW_gain=0
   do i=1,size(area)
      disp=dict_get_key(disp_con, element(int(ident(i))))
      vdW_gain=vdW_gain+(area(i)*disp%param)
   end do
   write(*,*) "EvdW: ", vdW_gain
   write(*,*) "Area: ", sum(area)

   thermo=param(10)*R*jtokcal*T
   write(*,*) "thermostatic correction: ", thermo

   !!! RING CORRECTION IS MISSING ATM
   gas_chem=-id_scr+thermo-vdW_gain!-ring_corr
   !write(*,*) gas_chem


end subroutine calcgas


function E_dd(c_hb,alpha,f_corr,s_hb,sv1,svt1,sv2,svt2,ident,element,atom1,atom2,id2,ele2)
   real(wp), intent(in) :: c_hb, alpha, f_corr,s_hb
   real(wp), intent(in) :: sv1, svt1, sv2, svt2
   character(2), dimension(:), intent(in) :: element
   integer, dimension(:), intent(in) :: ident
   character(2), dimension(:), intent(in),optional :: ele2
   integer, dimension(:), intent(in),optional :: id2

   ! LOCAL

   character(2), dimension(:), allocatable :: element2
   real(wp), dimension(:), allocatable :: ident2
   integer, intent(in) :: atom1, atom2
   real(wp) :: E_dd
   real(wp) :: svdo, svac
   real(wp) :: E_misfit, E_hb

   !! Setting up optional parameters at the beginning

   if (present(id2)) then
      allocate(ident2(size(id2)))
      ident2(:)=id2(:)
   else
      allocate(ident2(size(ident)))
      ident2(:)=ident(:)
   end if

   if (present(ele2)) then
      allocate(element2(size(ele2)))
      element2(:)=ele2(:)
   else
      allocate(element2(size(element)))
      element2=element(:)
   end if

   !! Set Acceptor and Donor

   svac=0
   svdo=0

   if (sv1 .GT. sv2) then
      svac=sv1
      svdo=sv2
   else
      svac=sv2
      svdo=sv1
   end if

   !! Start E_dd Calculation

   E_hb=0.0_wp
   E_misfit=0.0_wp
   E_hb=c_hb*max(0.0_wp,svac-s_hb)*min(0.0_wp,svdo+s_hb)
   E_misfit=(alpha/2)*(sv1+sv2)&
      &*((sv1+sv2)+f_corr*(svt1+svt2))

   E_dd=E_hb+E_misfit

end function E_dd

subroutine iterate_solvent(pot_di,sv,svt,area,T,ident,element,edd)
   use globals
   !real(wp), dimension(10), intent(in) :: param
   real(wp), dimension(:), allocatable, intent(in) :: sv, svt, area
   integer, allocatable, intent(in) :: ident(:)
   character(2), dimension(:), allocatable, intent(in) :: element
   real(wp), dimension(:), allocatable, intent(inout) :: pot_di
   real(wp), intent(in) :: edd(:, :)
   real(wp), intent(in) :: T

   real(wp) :: temppot, area_total
   integer :: i, j

   area_total = sum(area)
   !! For mixed solvent, mole fraction needs to be introduced in the following loop
   do j=1,size(pot_di)
      temppot=0.0_wp
      do i=1,size(pot_di)
         temppot = temppot + area(i)*exp(-Edd(i,j)+pot_di(i))
      end do
      !exit
      pot_di(j)=-log(temppot/area_total)
   end do

end subroutine iterate_solvent


subroutine compute_solute(sol_pot,solv_pot,sv_sol,svt_sol,sv_solv,svt_solv,area_sol,area_solv,T,chem_pot_sol,&
      &ident_sol,ident_solv,elem_sol,elem_solv)
   use globals
   !real(wp), dimension(10), intent(in) :: param
   real(wp), intent(out) :: chem_pot_sol
   real(wp), dimension(:), allocatable, intent(in) :: sv_sol, svt_sol,sv_solv,svt_solv,area_solv,area_sol
   integer, allocatable, intent(in), dimension(:) :: ident_sol, ident_solv
   real(wp), dimension(:), allocatable, intent(inout) :: solv_pot,sol_pot
   character(2), dimension(:), allocatable, intent(in) :: elem_sol, elem_solv
   real(wp), dimension(:), allocatable :: W_v
   real(wp), intent(in) :: T

   real(wp) :: temppot, beta,temp2
   integer :: i, j

   allocate(W_v(size(sv_sol)))
   allocate(sol_pot(size(sv_sol)))
   write(output_unit,'(5x,a)') &
      "Calculate Solvent-Solute Interaction based on the converged Solvent Profile."
   W_v(:)=0.0_wp
   beta=(R*Jtokcal*T)/param(7)
   temppot=0.0_wp
   sol_pot(:)=0
   !! For mixed solvent, mole fraction needs to be introduced in the following loop
   do j=1,size(sol_pot)
      do i=1,size(solv_pot)
         !  if (i .NE. j) then
         temppot=temppot+(area_solv(i)*exp((-E_dd&
            &(param(5),param(3),param(4),param(6),sv_sol(j),svt_sol(j),sv_solv(i),svt_solv(i),ident_sol,elem_sol,&
            j,i,ident_solv,elem_solv)&
            &/beta)+solv_pot(i)))
         !  end if
         W_v(j)=W_v(j)+area_solv(i)
      end do
      sol_pot(j)=-log(sum(area_solv)**(-1)*temppot)
      temppot=0.0_wp
   end do
   temppot=0.0_wp
   temp2=0
   do i=1,size(sol_pot)
      temppot=temppot+(area_sol(i)*sol_pot(i))
      !    write(*,*) area_sol(i), sol_pot(i), area_sol(i)*sol_pot(i)
   end do

   chem_pot_sol=temppot*beta-(param(wp)*R*Jtokcal*T*log(sum(area_solv)))
   !write(*,*) chem_pot_sol
   temppot=0
   do i=1,size(solv_pot)
      temppot=temppot+(area_solv(i)*solv_pot(i))
   end do
   !write(*,*) beta*temppot
   write(output_unit,'(5x,a)') "Done!", &
   ""
end subroutine compute_solute

subroutine compute_solvent(pot_di,sv,svt,area,T,max_cycle,conv_crit,ident,element)
   use globals, only : param, autokcal
   !real(wp), dimension(10), intent(in) :: param
   real(wp), dimension(:), allocatable :: sv, svt, area
   integer, allocatable :: ident(:)
   real(wp), dimension(:), allocatable, intent(inout) :: pot_di
   character(2), dimension(:), allocatable, intent(in) :: element
   real(wp), intent(in) :: T
   real(4), intent(in) :: conv_crit
   integer, intent(in) :: max_cycle

   integer :: iter
   logical :: converged
   real(wp), allocatable :: edd(:, :)
   type(broyden_mixer) :: mixer
   real(wp), parameter :: mixer_damping = 0.4_wp

   write(output_unit,'(5x,a)') "", & 
      "Converging Solvent Sigma Profile."

   converged = .false.
   allocate(pot_di(size(sv)))
   allocate(edd(size(pot_di),size(pot_di)))
   pot_di(:) = 0.0_wp
   iter = 0
   call calculate_edd(edd,sv,svt,T,ident,element)
   call new_broyden(mixer,max_cycle,size(pot_di),mixer_damping)
   do while (.not.converged)
      if (iter > 0) then
         call mixer%next
         call mixer%get(pot_di)
      end if
      iter=iter+1
      call mixer%set(pot_di)
      Call iterate_solvent(pot_di,sv,svt,area,T,ident,element,edd)
      call mixer%diff(pot_di)
      converged = mixer%get_error() < conv_crit/autokcal
      if (iter >= max_cycle) exit
   end do
   if (.not.converged) then
      write(*,*) "Error, chemical potential could not be converged"
      error stop
   end if
   write(output_unit,'(5x,a,i3)') "Done! Chemical potential converged after Cycle ",iter

end subroutine compute_solvent

subroutine calculate_edd(edd, sv, svt, T, ident, element)
   use globals, only : R, Jtokcal, param
   real(wp), intent(inout) :: edd(:, :)
   real(wp), intent(in) :: sv(:), svt(:)
   integer, intent(in) :: ident(:)
   character(2), intent(in) :: element(:)
   real(wp), intent(in) :: T

   integer :: i, j
   real(wp) :: temppot, beta

   beta=(R*Jtokcal*T)/param(7)

   do j=1,size(sv)
      do i=1,size(sv)
         edd(i, j) = E_dd&
            &(param(5),param(3),param(4),param(6),sv(j),svt(j),sv(i),svt(i),ident,element,j,i)&
            &/beta
      end do
   end do
end subroutine calculate_edd

end module crs
