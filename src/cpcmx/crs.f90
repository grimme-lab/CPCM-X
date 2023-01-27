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

module crs
   use mctc_env, only : wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use crs_broyden, only : broyden_mixer, new_broyden
   implicit none

contains


subroutine calcgas(E_cosmo,E_gas,id_scr,area,sv,su,pot,element,ident, T,eta,dE_is,dE_cc)
   use globals
   use element_dict
   real(wp), intent(out) :: id_scr
   real(wp), intent(in) :: T, E_cosmo, eta, E_gas
   real(wp),dimension(:),allocatable, intent(in) :: area, sv, su, pot
   real(wp), intent(out) :: dE_is, dE_cc
   integer, allocatable, intent(in) :: ident(:)
   character(2), dimension(:), allocatable, intent(in) :: element

   real(wp) :: dEreal, ediel, edielprime, vdW_gain, thermo, beta, avcorr
   integer :: dummy1, ioerror, i

   !> Check if gas phase energy exists.
   logical :: ex

   dEreal=(E_cosmo-E_gas)
   ediel=0
   edielprime=0
   do i=1,size(sv)
      ediel=ediel+(area(i)*pot(i)*su(i))
      edielprime=edielprime+(area(i)*pot(i)*sv(i))
   end do
   avcorr=(edielprime-ediel)/2.0_wp*0.8_wp
   dEreal=dEreal*autokcal
   id_scr=dEreal+avcorr*autokcal
   thermo=eta*R*jtokcal*T

   ! vdW Correction is replaced by the SMD Model.
   ! vdW_gain=0
   ! do i=1,size(area)
   !    disp=dict_get_key(disp_con, element(int(ident(i))))
   !    vdW_gain=vdW_gain+(area(i)*disp%param)
   ! end do
   ! write(*,*) "EvdW: ", vdW_gain
   ! write(*,*) "Area: ", sum(area)
   
   dG_is=dEreal-thermo
   dG_cc=avcorr*autokcal
   dE_is=(dEreal-thermo)/autokcal
   dE_cc=avcorr


end subroutine calcgas


pure function E_dd(c_hb,alpha,f_corr,s_hb,sv1,svt1,sv2,svt2)
   real(wp), intent(in) :: c_hb
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: f_corr
   real(wp), intent(in) :: s_hb

   real(wp), intent(in) :: sv1, svt1, sv2, svt2


   real(wp) :: E_dd
   real(wp) :: svdo, svac
   real(wp) :: E_misfit, E_hb

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

subroutine iterate_solvent(param,pot_di,sv,svt,area,T,edd)
   use type, only: parameter_type

   type(parameter_type), intent(in) :: param
   real(wp), dimension(:), allocatable, intent(in) :: sv, svt, area
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


subroutine compute_solute(param,sol_pot,solv_pot,sv_sol,svt_sol,sv_solv,svt_solv,area_sol,area_solv,T,chem_pot_sol)
   use globals, only: Jtokcal, R
   use type, only: parameter_type

   type(parameter_type) :: param
   real(wp), intent(out) :: chem_pot_sol
   real(wp), dimension(:), allocatable, intent(in) :: sv_sol, svt_sol,sv_solv,svt_solv,area_solv,area_sol
   real(wp), dimension(:), allocatable, intent(inout) :: solv_pot,sol_pot
   !real(wp), dimension(:), allocatable :: W_v
   real(wp), intent(in) :: T

   real(wp) :: temppot, beta,temp2
   integer :: i, j

   !allocate(W_v(size(sv_sol)))
   allocate(sol_pot(size(sv_sol)))
   !write(output_unit,'(5x,a)') &
   !   "Calculate Solvent-Solute Interaction based on the converged Solvent Profile."
   !W_v(:)=0.0_wp
   beta=(R*Jtokcal*T)/param%aeff
   temppot=0.0_wp
   sol_pot(:)=0
   !! For mixed solvent, mole fraction needs to be introduced in the following loop
   do j=1,size(sol_pot)
      do i=1,size(solv_pot)
         !  if (i .NE. j) then
         temppot=temppot+(area_solv(i)*exp((-E_dd&
            &(param%chb,param%aprime,param%fcorr,param%shb,sv_sol(j),svt_sol(j),sv_solv(i),svt_solv(i))&
            &/beta)+solv_pot(i)))
         !  end if
         !W_v(j)=W_v(j)+area_solv(i)
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

   chem_pot_sol=temppot*beta-(param%lambda*R*Jtokcal*T*log(sum(area_solv)))
   !write(*,*) chem_pot_sol
   temppot=0
   do i=1,size(solv_pot)
      temppot=temppot+(area_solv(i)*solv_pot(i))
   end do
   !write(*,*) beta*temppot
   !write(output_unit,'(5x,a)') "Done!", &
   !""
end subroutine compute_solute

subroutine compute_solvent(param,pot_di,sv,svt,area,T,max_cycle,conv_crit,error)
   use mctc_env, only: error_type, fatal_error
   use globals, only : autokcal, to_string
   use type, only : parameter_type

   type(parameter_type), intent(in) :: param
   real(wp), dimension(:), allocatable :: sv, svt, area
   real(wp), dimension(:), allocatable, intent(inout) :: pot_di
   real(wp), intent(in) :: T
   real(wp), intent(in) :: conv_crit
   integer, intent(in) :: max_cycle

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iter
   logical :: converged
   real(wp), allocatable :: edd(:, :)
   type(broyden_mixer) :: mixer
   real(wp), parameter :: mixer_damping = 0.4_wp

   !write(output_unit,'(5x,a)') "", & 
   !   "Converging Solvent Sigma Profile."

   converged = .false.
   allocate(pot_di(size(sv)))
   allocate(edd(size(pot_di),size(pot_di)))
   pot_di(:) = 0.0_wp
   iter = 0
   call calculate_edd(param,edd,sv,svt,T)
   call new_broyden(mixer,max_cycle,size(pot_di),mixer_damping)
   do while (.not.converged)
      if (iter > 0) then
         call mixer%next
         call mixer%get(pot_di)
      end if
      iter=iter+1
      call mixer%set(pot_di)
      Call iterate_solvent(param,pot_di,sv,svt,area,T,edd)
      call mixer%diff(pot_di)
      converged = mixer%get_error() < conv_crit/autokcal
      if (iter >= max_cycle) exit
   end do

   if (.not.converged) then
      Call fatal_error(error,"Solvent Sigma Profile did not converge after "&
         &//to_string(max_cycle)//" cycles.")
         return
   end if
   !write(output_unit,'(5x,a,i3)') "Done! Chemical potential converged after Cycle ",iter

end subroutine compute_solvent

subroutine calculate_edd(param,edd, sv, svt, T)
   use globals, only : R, Jtokcal
   use type, only: parameter_type
   type(parameter_type), intent(in) :: param
   real(wp), intent(inout) :: edd(:, :)
   real(wp), intent(in) :: sv(:), svt(:)
   real(wp), intent(in) :: T

   integer :: i, j
   real(wp) :: temppot, beta

   beta=(R*Jtokcal*T)/param%aeff

   do j=1,size(sv)
      do i=1,size(sv)
         edd(i, j) = E_dd&
            &(param%chb,param%aprime,param%fcorr,param%shb,sv(j),svt(j),sv(i),svt(i))&
            &/beta
      end do
   end do
end subroutine calculate_edd

subroutine state_correction(density,mass,T,dG_state)
   use globals, only: R, autokcal, jtokcal
   !> Density of the Solvent (kg/m^3)
   real(wp), intent(in) :: density
   !> Atomic Mass of the Solvent (a.u.)
   real(wp), intent(in) :: mass
   !> Temperature of the Sytem
   real(wp), intent(in) :: T

   !> Correction energy bar-mol/mol -> mol/L-mol/L
   real(wp), intent(out) :: dG_state

   !> Molar Volume of the Gas in L/mol
   real(wp) :: V_m

   V_m=(R*T)/100

   dG_state=(-R*(Jtokcal)*T*log((density*V_m)/mass))/autokcal

end subroutine state_correction


end module crs
