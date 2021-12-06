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

module sac_mod
   use mctc_env, only : wp
   implicit none

   contains

   subroutine sac_gas(E_cosmo,id_scr,area,sv,su,pot)
         use globals
         use element_dict
         real(wp), intent(out) :: id_scr
         real(wp), intent(in) ::E_cosmo
         real(wp),dimension(:),allocatable, intent(in) :: area, sv, su, pot!,ident
         !character(2), dimension(:), allocatable, intent(in) :: element
       !  type(DICT_STRUCT), pointer :: dispa_con, dispb_con
         !real(wp), dimension(10) :: param
       !  type(DICT_DATA) :: a_disp,b_disp
         real(wp) :: E_gas, dEreal, ediel, edielprime, vdW_gain, thermo, beta, avcorr
         integer :: dummy1, ioerror, i

         open(1,file="energy")
         read(1,*,iostat=ioerror)
         if (ioerror .NE. 0) then
            write(*,*) "Problem while reading energies (check energy file)."
            error stop
         else
            read(1,*) dummy1,E_gas
         end if
         dEreal=(E_cosmo-E_gas)
         ediel=0
         edielprime=0
         do i=1,size(sv)
            ediel=ediel+(area(i)*pot(i)*su(i))
            edielprime=edielprime+(area(i)*pot(i)*sv(i))
         end do
         avcorr=(edielprime-ediel)/2.0_wp*0.8_wp
         write(*,*) "E_COSMO", E_cosmo*autokcal
         write(*,*) "E_COSMO+dE: ", (E_cosmo+avcorr)*autokcal
         write(*,*) "E_gas: ", E_gas*autokcal
         dEreal=dEreal*autokcal
         id_scr=dEreal+avcorr*autokcal
         write(*,*) "E_COSMO-E_gas", (E_cosmo-E_gas)*autokcal
         write(*,*) "E_COSMO-E_gas+dE: ", (E_cosmo-E_gas+avcorr)*autokcal
         write(*,*) "Ediel: ", ediel/2*autokcal
         write(*,*) "Averaging corr dE: ", avcorr*autokcal

         vdW_gain=0
!         do i=1,size(area)
 !           a_disp=dict_get_key(disp_cona, element(int(ident(i))))
 !           b_disp=dict_get_key(disp_conb, element(int(ident(i))))
  !          vdW_gain=vdW_gain+(area(i)*(a_disp%param*SysTemp+b_disp%param))
  !       end do
     !    write(*,*) "EvdW: ", vdW_gain
         write(*,*) "Area: ", sum(area)
        ! thermo=param(10)*R*jtokcal*SysTemp
    !     write(*,*) "thermostatic correction: ", thermo

         !!! RING CORRECTION IS MISSING ATM
      !   gas_chem=-id_scr+thermo!-vdW_gain!-ring_corr
       !  write(*,*) gas_chem
        dG_is=(E_cosmo-E_gas)*autokcal
        dG_cc=avcorr*autokcal

         if (ML) then
            open(5,file='ML.energy')
            write(5,'(F0.10A)',advance='no') dG_is,","
            write(5,'(F0.10A)',advance='no') dG_cc,","
            write(5,'(F0.10A)',advance='no') dG_disp,","
            close(5)
         end if
   end subroutine sac_gas

function E_dd1(sigma1,sigma2)
   use globals
   implicit none

   real(wp), intent(in) :: sigma1,sigma2
   real(wp), parameter :: EPS=3.667_wp, e0=2.395E-4
   real(wp) :: E_dd1,svdo,svac,E_misfit,E_hb, fpol, alpha, alphaprime,aef,s_hb,c_hb

   aef=param(5)
   c_hb=param(6)
   s_hb=abs(param(7)) !Sigma_Hb can't be negative

   fpol=(EPS-1.0_wp)/(EPS+0.5_wp)
   alpha=(0.3_wp*aef**(1.5))/e0
   alphaprime=fpol*alpha!param(8) !alphaprime is not really a parameter

   svac=0
   svdo=0

   if (sigma1 .GE. sigma2) then
      svac=sigma1
      svdo=sigma2
   else
      svac=sigma2
      svdo=sigma1
   end if
   E_hb=0.0_wp
   E_misfit=0.0_wp
   E_hb=c_hb*max(0.0_wp,svac-s_hb)*min(0.0_wp,svdo+s_hb)
  ! write(*,*) sigma1,sigma2
   E_misfit=(alphaprime/2.0_wp)*((sigma1+sigma2)**2.0_wp)
   E_dd1=E_misfit+E_hb
 !  Write(*,*) E_dd1

end function E_dd1


function E_dd3(sigma1,sigma2,t,s)
   use globals
   implicit none
   integer,intent(in) :: t,s
   real(wp), intent(in) :: sigma1,sigma2

   real(wp) :: E_dd3,E_es,E_hb,c_hb,c_es,A_es,B_es
   real(wp) :: c_ohoh,c_otot,c_ohot

   ! t and s are the sigma profile types (1=NH, 2=OH, 3=OT)


   c_ohoh=param(6)
   c_otot=param(7)
   c_ohot=param(8)
  c_hb=0.0_wp
   if ((sigma1*sigma2) .LT. 0.0_wp) then
      if ((s .EQ. 2) .AND. (t .EQ. 2)) then
         c_hb=c_ohoh
      else if ((s .EQ. 3) .AND. (t .EQ. 3)) then
         c_hb=c_otot
      else if (((s .EQ. 2) .AND. (t .EQ. 3)) .OR. ((s .EQ. 3) .AND. (t .EQ. 2))) then
         c_hb=c_ohot
      end if
   end if

   A_es=param(9)
   B_es=param(10)
   c_es=A_es+(B_es/(SysTemp**2)) !Temperaturabhängige elektrostatische interactiom

   E_hb=0.0_wp
   E_es=0.0_wp
   E_hb=c_hb*((sigma1-sigma2)**2)
   E_es=c_es*((sigma1+sigma2)**2)
   E_dd3=E_es-E_hb

end function E_dd3

subroutine sac_2005(profil,profil2,vcosmo1,vcosmo2,z1,z2)
   use globals
   use mctc_env, only: wp
   implicit none

   real(wp), dimension(0:50) :: profil,profil2
   !real(wp), dimension(1:9), intent(in) :: param
  ! real(wp), dimension(1:2,0:49) :: sigma_profiles
   real(wp) :: z1,z2
   real(wp) :: gam(0:50),maxsig,punit,profile(0:50), gam_saved(0:50),gam_sol(0:50)
   real(wp) :: gamma_solv, gamma_sol,gamma_test,gamma_test2, summ, mix_prof(0:50), mix_gam(0:50)
   real(wp) :: T, VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2
   real(wp) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation


  ! integer :: comp_num
   integer :: i,j, cycles
   logical :: not_conv

   z(1)=z1
   z(2)=z2

  ! comp_num=2

   VNORM=param(3)
   ANORM=param(2)
   not_conv=.TRUE.

   !! Calculate the SIGMA Profile of the mixture
  ! if (ML) open(5,file='ML.prof')
   do i=0,50
      mix_prof(i)=(z(1)*profil(i)+z(2)*profil2(i))/(z(1)*sum(profil)+z(2)*sum(profil2))
  !    if (ML) write(5,'(F0.10A)',advance='no') mix_prof(i),","
   end do
   if (ML) then
   !   close(5)
      open(5,file='ML.gamma')
   end if


      !! Mixed Activity Coefficient
      T=SysTemp
      maxsig=0.025_wp
      punit=0.001_wp
      mix_gam(:)=1.0
      gam_saved(:)=1.0
      summ=0.0_wp

      cycles=0
      do while (not_conv)
      cycles=cycles+1
         gam_saved(:)=mix_gam(:)
         do i=0,50
            do j=0,50
               summ=summ+mix_prof(j)*gam_saved(j)*exp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_wp*R*Jtokcal))
            end do
            mix_gam(i)=exp(-log(summ))
            mix_gam(i)=(mix_gam(i)+gam_saved(i))/2.0_wp
            summ=0.0_wp
         end do
         not_conv=.false.
         do i=0,50
            if (abs((mix_gam(i)-gam_saved(i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
        if (cycles .gt. 10000) then
          error stop "Error, mixed profile not converged"
        end if
      end do

      !! Pure Activity Coefficient of 1
      profile(:)=0.0_wp
      do i=0,50
         profile(i)=profil(i)/sum(profil)
      end do

      gam(:)=1.0
      gam_saved(:)=1.0
      summ=0.0_wp
      not_conv=.true.
      cycles=0
      do while (not_conv)
      cycles=cycles+1
         gam_saved(:)=gam(:)
         do i=0,50
            do j=0,50
               summ=summ+profile(j)*gam_saved(j)*exp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_wp*R*Jtokcal))
            end do
            gam(i)=exp(-log(summ))
            gam(i)=(gam(i)+gam_saved(i))/2.0_wp
            summ=0.0_wp
         end do
         not_conv=.false.
         do i=0,50
            if (abs((gam(i)-gam_saved(i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
        if (cycles .gt. 10000) then
          error stop "Error, pure profile 1 not converged"
        end if
      end do

      !! Pure Activity Coefficient of 2

      profile(:)=0.0_wp
      do i=0,50
         profile(i)=profil2(i)/sum(profil2)
      end do

      gam_sol(:)=1.0
      gam_saved(:)=1.0
      summ=0.0_wp
      not_conv=.true.
      cycles=0
      do while (not_conv)
      cycles=cycles+1
         gam_saved(:)=gam_sol(:)
         do i=0,50
            do j=0,50
               summ=summ+profile(j)*gam_saved(j)*exp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_wp*R*Jtokcal))
            end do
            gam_sol(i)=exp(-log(summ))
            gam_sol(i)=(gam_sol(i)+gam_saved(i))/2.0_wp
            summ=0.0_wp
         end do
         not_conv=.false.
         do i=0,50
            if (abs((gam_sol(i)-gam_saved(i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
        if (cycles .gt. 10000) then
          error stop "Error, pure profile 2 not converged"
        end if
      end do

      !! Staverman-Guggenheim equation
      coord=int(param(4))

      RNORM(1) = VCOSMO1/VNORM
      QNORM(1) = sum(profil)/ANORM
      RNORM(2) = VCOSMO2/VNORM
      QNORM(2) = sum(profil2)/ANORM

      bt=z(1)*QNORM(1)+z(2)*QNORM(2)
      bp=z(1)*RNORM(1)+z(2)*RNORM(2)


      do i=1,2
         Theta(i)=(z(i)*QNORM(i))/(bt)
         Phi(i)=(z(i)*RNORM(i))/(bp)
         L(i)=(coord/2.0_wp)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_wp)
      end do
      gammasg(1)=log(phi(1)/z(1))+(coord/2.0_wp)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
         &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))

      gammasg(2)=log(phi(2)/z(2))+(coord/2.0_wp)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
         &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))


      write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
      gamma_solv=0.0_wp
      gamma_sol=0.0_wp
      gamma_test=0.0_wp
      gamma_test2=0.0_wp
      do i=0,50
         if (ML) then
            write(5,'(F0.10A)',advance='no') mix_gam(i),","
            ! write(5,'(F0.10A)',advance='no') gam(i),","
            ! write(5,'(F0.10A)',advance='no') gam_sol(i),","
         end if
         !gamma_test=gamma_test+(profil(i)/param(5)*log(gam(i)))
         !gamma_test2=gamma_test2+(profil2(i)/param(5)*log(gam_sol(i)))
         gamma_solv=gamma_solv+(profil(i)/param(5)*(log(mix_gam(i)/gam(i))))
         gamma_sol=gamma_sol+(profil2(i)/param(5)*log(mix_gam(i)/gam_sol(i)))
      end do
      if (ML) close(5)
      !write(*,*) gamma_test, gamma_test2
      gamma_solv=exp(gammasg(1)+gamma_solv)
      gamma_sol=exp(gammasg(2)+gamma_sol)
      write(*,*) "Results for Mixture with Compound 1 x= ",z(1)," and Compound 2 x= ",z(2),"."
      write(*,*) "Gamma(1)= ",gamma_solv, "Gamma(2)= ", gamma_sol
      write(*,*) "lnGamma(1)= ", log(gamma_solv),"lnGamma(2)= ", log(gamma_sol)
      dG_res=log(gamma_sol)*SysTemp*R*Jtokcal


end subroutine sac_2005

subroutine sac_2010(profil,profil2,vcosmo1,vcosmo2)
   use globals
   implicit none

   real(wp), dimension(3,0:50) :: profil,profil2
   !real(wp), dimension(1:9), intent(in) :: param
  ! real(wp), dimension(1:2,0:49) :: sigma_profiles

   real(wp) :: gam(3,0:50),maxsig,punit,profile(3,0:50), gam_saved(3,0:50),gam_sol(3,0:50)
   real(wp) :: gamma_solv, gamma_sol,gamma_test, summ, mix_prof(3,0:50), mix_gam(3,0:50)
   real(wp) :: VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2
   real(wp) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation


  ! integer :: comp_num
   integer :: i,j,s,t, cycles
   logical :: not_conv

   z(1)=0.995_wp
   z(2)=0.005_wp

  ! comp_num=2

   VNORM=param(3)
   ANORM=param(2)
   not_conv=.TRUE.
   if (ML) open(5,file="ML.prof")
   !! Calculate the SIGMA Profile of the mixture
   do s=1,3
      do i=0,50
         mix_prof(s,i)=(z(1)*profil(s,i)+z(2)*profil2(s,i))/(z(1)*sum(profil)+z(2)*sum(profil2))
         if (ML) write(5,'(F0.10A)',advance='no') mix_prof(s,i),","
      end do
   end do

   !! Mixed Activity Coefficient
   ! T=SysTemp
   maxsig=0.025_wp
   punit=0.001_wp
   mix_gam=1.0
   gam_saved=1.0
   summ=0.0_wp
  cycles=0
   do while (not_conv)
   cycles=cycles+1
      gam_saved(:,:)=mix_gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+mix_prof(s,j)*gam_saved(t,j)*exp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            mix_gam(t,i)=exp(-log(summ))
            mix_gam(t,i)=(mix_gam(t,i)+gam_saved(t,i))/2.0_wp
            summ=0.0_wp
         end do
      end do
      not_conv=.false.
 !     write(*,*) mix_gam

      do t=1,3
         do i=0,50
            if (abs((mix_gam(t,i)-gam_saved(t,i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
      end do
      if (cycles .gt. 10000) then
        error stop "Error, not converged"
      end if
   end do

   if (ML) then
   !   do t=1,3
   !      do i=0,50
   !         write(5,'(F12.10A)',advance='no') mix_gam(t,i),","
   !      end do
   !   end do
      close(5)
   end if
   !! Pure Activity Coefficient of 1
   profile=0.0_wp
   do t=1,3
      do i=0,50
         profile(t,i)=profil(t,i)/sum(profil)
      end do
   end do

   gam=1.0
   gam_saved=1.0
   summ=0.0_wp
   not_conv=.true.
   cycles=0
   do while (not_conv)
   cycles=cycles+1
      gam_saved(:,:)=gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*exp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam(t,i)=exp(-log(summ))
            gam(t,i)=(gam(t,i)+gam_saved(t,i))/2.0_wp
            summ=0.0_wp
         end do
      end do
      not_conv=.false.
      do t=1,3
         do i=0,50
            if (abs((gam(t,i)-gam_saved(t,i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
      end do
      if (cycles .gt. 10000) then
        error stop "Error, not converged"
      end if
   end do


   !! Pure Activity Coefficient of 2

   profile=0.0_wp
   cycles=0
   do t=1,3
      do i=0,50
         profile(t,i)=profil2(t,i)/sum(profil2)
      end do
   end do

   gam_sol=1.0
   gam_saved=1.0
   summ=0.0_wp
   not_conv=.true.

   do while (not_conv)
   cycles=cycles+1
      gam_saved(:,:)=gam_sol(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*exp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam_sol(t,i)=exp(-log(summ))
            gam_sol(t,i)=(gam_sol(t,i)+gam_saved(t,i))/2.0_wp
            summ=0.0_wp
         end do
      end do
      not_conv=.false.
      do t=1,3
         do i=0,50
            if (abs((gam_sol(t,i)-gam_saved(t,i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
      end do
      if (cycles .gt. 10000) then
        error stop "Error, not converged"
      end if
   end do

   !! Staverman-Guggenheim equation
   coord=int(param(4))

   RNORM(1) = VCOSMO1/VNORM
   QNORM(1) = sum(profil)/ANORM
   RNORM(2) = VCOSMO2/VNORM
   QNORM(2) = sum(profil2)/ANORM

   bt=z(1)*QNORM(1)+z(2)*QNORM(2)
   bp=z(1)*RNORM(1)+z(2)*RNORM(2)


   do i=1,2
      Theta(i)=(z(i)*QNORM(i))/(bt)
      Phi(i)=(z(i)*RNORM(i))/(bp)
      L(i)=(coord/2.0_wp)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_wp)
   end do
   gammasg(1)=log(phi(1)/z(1))+(coord/2.0_wp)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
      &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))

   gammasg(2)=log(phi(2)/z(2))+(coord/2.0_wp)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
      &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))


   write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
   gamma_solv=0.0_wp
   gamma_sol=0.0_wp
   gamma_test=0.0_wp
   do t=1,3
      do i=0,50
         gamma_test=gamma_test+(profil2(t,i)/param(5)*log(gam_sol(t,i)))
         gamma_solv=gamma_solv+(profil(t,i)/param(5)*(log(mix_gam(t,i)/gam(t,i))))
         gamma_sol=gamma_sol+(profil2(t,i)/param(5)*log(mix_gam(t,i)/gam_sol(t,i)))
      end do
   end do
   write(*,*) gamma_test, gamma_sol, gammasg(2)
   gamma_solv=exp(gammasg(1)+gamma_solv)
   gamma_sol=exp(gammasg(2)+gamma_sol)
   write(*,*) "Results for Mixture with Compound 1 x= ",z(1)," and Compound 2 x= ",z(2),"."
   write(*,*) "Gamma(1)= ",gamma_solv, "Gamma(2)= ", gamma_sol
   write(*,*) "lnGamma(1)= ", log(gamma_solv),"lnGamma(2)= ", log(gamma_sol)
  ! write(*,*) param(10)
   !write(*,*) gammasg(1), gammasg(2)

   dG_res=log(gamma_sol)*SysTemp*R*Jtokcal

end subroutine sac_2010


subroutine sac_2013(profil,profil2,vcosmo1,vcosmo2,sac_disp)
   use globals
   implicit none

   real(wp), dimension(3,0:50) :: profil,profil2
   !real(wp), dimension(1:9), intent(in) :: param
  ! real(wp), dimension(1:2,0:49) :: sigma_profiles
   real(wp), dimension(:) :: sac_disp

   real(wp) :: gam(3,0:50),maxsig,punit,profile(3,0:50), gam_saved(3,0:50),gam_sol(3,0:50)
   real(wp) :: gamma_solv, gamma_sol,gamma_test, summ, mix_prof(3,0:50), mix_gam(3,0:50)
   real(wp) :: VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2, A, omega
   real(wp) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation
   real(wp) :: gammadisp(2)

  ! integer :: comp_num
   integer :: i,j,s,t
   logical :: not_conv

   omega=0.27027 !Needs some dispersion flags, should be negative for water + hb_only etc.
   z(1)=0.995_wp
   z(2)=0.005_wp

  ! comp_num=2

   VNORM=param(3)
   ANORM=param(2)
   not_conv=.TRUE.

   !! Calculate the SIGMA Profile of the mixture
   do s=1,3
      do i=0,50
         mix_prof(s,i)=(z(1)*profil(s,i)+z(2)*profil2(s,i))/(z(1)*sum(profil)+z(2)*sum(profil2))
      end do
   end do

   !! Mixed Activity Coefficient
   ! T=SysTemp
   maxsig=0.025_wp
   punit=0.001_wp
   mix_gam=1.0
   gam_saved=1.0
   summ=0.0_wp


   do while (not_conv)
      gam_saved(:,:)=mix_gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+mix_prof(s,j)*gam_saved(t,j)*exp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            mix_gam(t,i)=exp(-log(summ))
            mix_gam(t,i)=(mix_gam(t,i)+gam_saved(t,i))/2.0_wp
            summ=0.0_wp
         end do
      end do
      not_conv=.false.
 !     write(*,*) mix_gam

      do t=1,3
         do i=0,50
            if (abs((mix_gam(t,i)-gam_saved(t,i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
      end do
   end do
!stop
   !! Pure Activity Coefficient of 1
   profile=0.0_wp
   do t=1,3
      do i=0,50
         profile(t,i)=profil(t,i)/sum(profil)
      end do
   end do

   gam=1.0
   gam_saved=1.0
   summ=0.0_wp
   not_conv=.true.

   do while (not_conv)
      gam_saved(:,:)=gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*exp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam(t,i)=exp(-log(summ))
            gam(t,i)=(gam(t,i)+gam_saved(t,i))/2.0_wp
            summ=0.0_wp
         end do
      end do
      not_conv=.false.
      do t=1,3
         do i=0,50
            if (abs((gam(t,i)-gam_saved(t,i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
      end do
   end do


   !! Pure Activity Coefficient of 2

   profile=0.0_wp
   do t=1,3
      do i=0,50
         profile(t,i)=profil2(t,i)/sum(profil2)
      end do
   end do

   gam_sol=1.0
   gam_saved=1.0
   summ=0.0_wp
   not_conv=.true.

   do while (not_conv)
      gam_saved(:,:)=gam_sol(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*exp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam_sol(t,i)=exp(-log(summ))
            gam_sol(t,i)=(gam_sol(t,i)+gam_saved(t,i))/2.0_wp
            summ=0.0_wp
         end do
      end do
      not_conv=.false.
      do t=1,3
         do i=0,50
            if (abs((gam_sol(t,i)-gam_saved(t,i))) .LT. 0.000001) then
               cycle
            else
               not_conv=.true.
               exit
            end if
         end do
      end do
   end do

   !! Staverman-Guggenheim equation
   coord=int(param(4))

   RNORM(1) = VCOSMO1/VNORM
   QNORM(1) = sum(profil)/ANORM
   RNORM(2) = VCOSMO2/VNORM
   QNORM(2) = sum(profil2)/ANORM

   bt=z(1)*QNORM(1)+z(2)*QNORM(2)
   bp=z(1)*RNORM(1)+z(2)*RNORM(2)


   do i=1,2
      Theta(i)=(z(i)*QNORM(i))/(bt)
      Phi(i)=(z(i)*RNORM(i))/(bp)
      L(i)=(coord/2.0_wp)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_wp)
   end do
   gammasg(1)=log(phi(1)/z(1))+(coord/2.0_wp)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
      &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))

   gammasg(2)=log(phi(2)/z(2))+(coord/2.0_wp)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
      &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))


   !! Additional dispersion correction
   A=omega*(0.5*(sac_disp(1)+sac_disp(2))-sqrt(sac_disp(1)*sac_disp(2)))
   gammadisp(1)=A*z(2)*z(2)
   gammadisp(2)=A*z(1)*z(1)

   write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
   gamma_solv=0.0_wp
   gamma_sol=0.0_wp
   do t=1,3
      do i=0,50
         gamma_solv=gamma_solv+(profil(t,i)/param(5)*(log(mix_gam(t,i)/gam(t,i))))
         gamma_sol=gamma_sol+(profil2(t,i)/param(5)*log(mix_gam(t,i)/gam_sol(t,i)))
      end do
   end do
   !write(*,*) gamma_solv, gamma_sol
   gamma_solv=exp(gammasg(1)+gamma_solv+gammadisp(1))
   gamma_sol=exp(gammasg(2)+gamma_sol+gammadisp(2))
   write(*,*) "Results for Mixture with Compound 1 x= ",z(1)," and Compound 2 x= ",z(2),"."
   write(*,*) "Gamma(1)= ",gamma_solv, "Gamma(2)= ", gamma_sol
   write(*,*) "lnGamma(1)= ", log(gamma_solv),"lnGamma(2)= ", log(gamma_sol)
  ! write(*,*) param(10)
   !write(*,*) gammasg(1), gammasg(2)

   dG_res=log(gamma_sol)*SysTemp*R*Jtokcal
end subroutine sac_2013


subroutine sac2013_disp(nam,is_bonded,ident,elements,disp_con,sac_disp)
   use element_dict
   implicit none

   logical, dimension(:,:), allocatable, intent(in) :: is_bonded

   real(wp), dimension(:), allocatable, intent(in) :: ident
   character(2), dimension(:), allocatable, intent(in) :: elements
   character(len=*) :: nam
   type(DICT_STRUCT), pointer :: disp_con

   type(DICT_DATA) :: disp

   character(2) :: symbol
   character(2) :: bond_string
   character(3) :: data_string

   real(wp), intent(out) :: sac_disp

   real(wp) :: atom_disp

   integer :: i,j, bond_count, disp_atoms

   sac_disp=0.0_wp
   disp_atoms=0

   do i=1,int(maxval(ident))
   bond_count=0
   symbol=elements(i)

   select case (trim(symbol))
      case default
         do j=1,int(maxval(ident))
            if (is_bonded(i,j)) then
               bond_count=bond_count+1
            end if
         end do

         write(bond_string,"(I1)") bond_count
         data_string=trim(symbol)//trim(bond_string)
      case ('h')
         bond_string='t'
         if (nam .EQ. 'h2o') then
            bond_string='2'
         else
            do j=1,int(maxval(ident))
               if (is_bonded(i,j)) then
                  if (elements(j) .EQ. 'o') then
                     bond_string='o'
                     exit
                  else if (elements(j) .EQ. 'n') then
                     bond_string='n'
                     exit
                  end if
               end if
            end do
         end if
         data_string=trim(symbol)//trim(bond_string)
   end select

   disp=dict_get_key(disp_con, trim(data_string))
   atom_disp=disp%param
   sac_disp=sac_disp+atom_disp
   if (atom_disp .NE. 0) disp_atoms=disp_atoms+1

   !write(*,*) data_string, atom_disp
   end do

   sac_disp=sac_disp/disp_atoms
   !write(*,*) sac_disp




end subroutine sac2013_disp

end module sac_mod

