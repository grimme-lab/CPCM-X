module sac_mod

   contains
function E_dd1(sigma1,sigma2)
   use globals
   implicit none
   real(8), intent(in) :: sigma1,sigma2
   !real(8), parameter :: EPS=3.667_8, e0=2.395E-4
   real(8) :: E_dd1,svdo,svac,E_misfit,E_hb, fpol, alpha, alphaprime,aef,s_hb,c_hb

   aef=param(5)
   c_hb=param(6)
   s_hb=param(7)

  ! fpol=(EPS-1.0_8)/(EPS+0.5_8)
  ! alpha=(0.3_8*aef**(1.5))/e0
   alphaprime=param(8)
   
   svac=0
   svdo=0

   if (sigma1 .GE. sigma2) then
      svac=sigma1
      svdo=sigma2
   else
      svac=sigma2
      svdo=sigma1
   end if
   E_hb=0.0_8
   E_misfit=0.0_8
   E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
  ! write(*,*) sigma1,sigma2
   E_misfit=alphaprime/2*(sigma1+sigma2)**2.0_8
   E_dd1=E_misfit+E_hb
 !  Write(*,*) E_dd1

end function E_dd1


function E_dd3(sigma1,sigma2,t,s)
   use globals
   implicit none
   integer,intent(in) :: t,s
   real(8), intent(in) :: sigma1,sigma2

   real(8) :: E_dd3,E_es,E_hb,c_hb,c_es,A_es,B_es
   real(8) :: c_ohoh,c_otot,c_ohot

   ! t and s are the sigma profile types (1=NH, 2=OH, 3=OT)


   c_ohoh=param(6)
   c_otot=param(7)
   c_ohot=param(8)
  c_hb=0.0_8
   if ((sigma1*sigma2) .LT. 0.0_8) then
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

   E_hb=0.0_8
   E_es=0.0_8
   E_hb=c_hb*((sigma1-sigma2)**2)
   E_es=c_es*((sigma1+sigma2)**2)
   E_dd3=E_es-E_hb

end function E_dd3

subroutine sac_2005(profil,profil2,vcosmo1,vcosmo2)
   use globals
   implicit none

   real(8), dimension(0:50) :: profil,profil2
   !real(8), dimension(1:9), intent(in) :: param
  ! real(8), dimension(1:2,0:49) :: sigma_profiles

   real(8) :: gam(0:50),maxsig,punit,profile(0:50), gam_saved(0:50),gam_sol(0:50)
   real(8) :: gamma_solv, gamma_sol,gamma_test, summ, mix_prof(0:50), mix_gam(0:50)
   real(8) :: T, VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2
   real(8) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation

   
  ! integer :: comp_num
   integer :: i,j
   logical :: not_conv

   z(1)=0.995_8
   z(2)=0.005_8
   
  ! comp_num=2

   VNORM=param(3)
   ANORM=param(2)
   not_conv=.TRUE.

   !! Calculate the SIGMA Profile of the mixture

   do i=0,50
      mix_prof(i)=(z(1)*profil(i)+z(2)*profil2(i))/(z(1)*sum(profil)+z(2)*sum(profil2))
   end do

   !! Mixed Activity Coefficient
   T=SysTemp
   maxsig=0.025_8
   punit=0.001_8
   mix_gam(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8


   do while (not_conv)
      gam_saved(:)=mix_gam(:)
      do i=0,50
         do j=0,50
            summ=summ+mix_prof(j)*gam_saved(j)*dexp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_8*R*Jtokcal))
         end do
         mix_gam(i)=exp(-log(summ))
         mix_gam(i)=(mix_gam(i)+gam_saved(i))/2.0_8
         summ=0.0_8
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
   end do

   !! Pure Activity Coefficient of 1
   profile(:)=0.0_8
   do i=0,50
      profile(i)=profil(i)/sum(profil)
   end do

   gam(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8
   not_conv=.true. 
   do while (not_conv)
      gam_saved(:)=gam(:)
      do i=0,50
         do j=0,50
            summ=summ+profile(j)*gam_saved(j)*dexp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_8*R*Jtokcal))
         end do
         gam(i)=exp(-log(summ))
         gam(i)=(gam(i)+gam_saved(i))/2.0_8
         summ=0.0_8
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
   end do

   !! Pure Activity Coefficient of 2

   profile(:)=0.0_8
   do i=0,50
      profile(i)=profil2(i)/sum(profil2)
   end do
   
   gam_sol(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8
   not_conv=.true.
   do while (not_conv)
      gam_saved(:)=gam_sol(:)
      do i=0,50
         do j=0,50
            summ=summ+profile(j)*gam_saved(j)*dexp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_8*R*Jtokcal))
         end do
         gam_sol(i)=exp(-log(summ))
         gam_sol(i)=(gam_sol(i)+gam_saved(i))/2.0_8
         summ=0.0_8
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
      L(i)=(coord/2.0_8)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_8)
   end do
   gammasg(1)=log(phi(1)/z(1))+(coord/2.0_8)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
      &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))
   
   gammasg(2)=log(phi(2)/z(2))+(coord/2.0_8)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
      &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))


   write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
   gamma_solv=0.0_8
   gamma_sol=0.0_8
   do i=0,50
      gamma_solv=gamma_solv+(profil(i)/param(5)*(log(mix_gam(i)/gam(i))))
      gamma_sol=gamma_sol+(profil2(i)/param(5)*log(mix_gam(i)/gam_sol(i)))
   end do
   gamma_solv=exp(gammasg(1)+gamma_solv)
   gamma_sol=exp(gammasg(2)+gamma_sol)
   write(*,*) "Results for Mixture with Compound 1 x= ",z(1)," and Compound 2 x= ",z(2),"."
   write(*,*) "Gamma(1)= ",gamma_solv, "Gamma(2)= ", gamma_sol
   write(*,*) "lnGamma(1)= ", log(gamma_solv),"lnGamma(2)= ", log(gamma_sol)



end subroutine sac_2005

subroutine sac_2010(profil,profil2,vcosmo1,vcosmo2)
   use globals
   implicit none

   real(8), dimension(3,0:50) :: profil,profil2
   !real(8), dimension(1:9), intent(in) :: param
  ! real(8), dimension(1:2,0:49) :: sigma_profiles

   real(8) :: gam(3,0:50),maxsig,punit,profile(3,0:50), gam_saved(3,0:50),gam_sol(3,0:50)
   real(8) :: gamma_solv, gamma_sol,gamma_test, summ, mix_prof(3,0:50), mix_gam(3,0:50)
   real(8) :: VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2
   real(8) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation

   
  ! integer :: comp_num
   integer :: i,j,s,t
   logical :: not_conv

   z(1)=0.995_8
   z(2)=0.005_8
   
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
   maxsig=0.025_8
   punit=0.001_8
   mix_gam=1.0
   gam_saved=1.0
   summ=0.0_8


   do while (not_conv)
      gam_saved(:,:)=mix_gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+mix_prof(s,j)*gam_saved(t,j)*dexp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            mix_gam(t,i)=exp(-log(summ))
            mix_gam(t,i)=(mix_gam(t,i)+gam_saved(t,i))/2.0_8
            summ=0.0_8
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
   profile=0.0_8
   do t=1,3
      do i=0,50
         profile(t,i)=profil(t,i)/sum(profil)
      end do
   end do

   gam=1.0
   gam_saved=1.0
   summ=0.0_8
   not_conv=.true.
   
   do while (not_conv)
      gam_saved(:,:)=gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*dexp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam(t,i)=exp(-log(summ))
            gam(t,i)=(gam(t,i)+gam_saved(t,i))/2.0_8
            summ=0.0_8
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

   profile=0.0_8
   do t=1,3
      do i=0,50
         profile(t,i)=profil2(t,i)/sum(profil2)
      end do
   end do
   
   gam_sol=1.0
   gam_saved=1.0
   summ=0.0_8
   not_conv=.true.

   do while (not_conv)
      gam_saved(:,:)=gam_sol(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*dexp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam_sol(t,i)=exp(-log(summ))
            gam_sol(t,i)=(gam_sol(t,i)+gam_saved(t,i))/2.0_8
            summ=0.0_8
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
      L(i)=(coord/2.0_8)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_8)
   end do
   gammasg(1)=log(phi(1)/z(1))+(coord/2.0_8)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
      &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))
   
   gammasg(2)=log(phi(2)/z(2))+(coord/2.0_8)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
      &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))


   write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
   gamma_solv=0.0_8
   gamma_sol=0.0_8
   do t=1,3
      do i=0,50
         gamma_solv=gamma_solv+(profil(t,i)/param(5)*(log(mix_gam(t,i)/gam(t,i))))
         gamma_sol=gamma_sol+(profil2(t,i)/param(5)*log(mix_gam(t,i)/gam_sol(t,i)))
      end do
   end do
   !write(*,*) gamma_solv, gamma_sol
   gamma_solv=exp(gammasg(1)+gamma_solv)
   gamma_sol=exp(gammasg(2)+gamma_sol)
   write(*,*) "Results for Mixture with Compound 1 x= ",z(1)," and Compound 2 x= ",z(2),"."
   write(*,*) "Gamma(1)= ",gamma_solv, "Gamma(2)= ", gamma_sol
   write(*,*) "lnGamma(1)= ", log(gamma_solv),"lnGamma(2)= ", log(gamma_sol)
  ! write(*,*) param(10) 
   !write(*,*) gammasg(1), gammasg(2)

end subroutine sac_2010


subroutine sac_2013(profil,profil2,vcosmo1,vcosmo2,sac_disp)
   use globals
   implicit none

   real(8), dimension(3,0:50) :: profil,profil2
   !real(8), dimension(1:9), intent(in) :: param
  ! real(8), dimension(1:2,0:49) :: sigma_profiles
   real(8), dimension(:) :: sac_disp

   real(8) :: gam(3,0:50),maxsig,punit,profile(3,0:50), gam_saved(3,0:50),gam_sol(3,0:50)
   real(8) :: gamma_solv, gamma_sol,gamma_test, summ, mix_prof(3,0:50), mix_gam(3,0:50)
   real(8) :: VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2, A, omega
   real(8) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation
   real(8) :: gammadisp(2)
   
  ! integer :: comp_num
   integer :: i,j,s,t
   logical :: not_conv

   omega=0.27027 !Needs some dispersion flags, should be negative for water + hb_only etc.
   z(1)=0.995_8
   z(2)=0.005_8
   
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
   maxsig=0.025_8
   punit=0.001_8
   mix_gam=1.0
   gam_saved=1.0
   summ=0.0_8


   do while (not_conv)
      gam_saved(:,:)=mix_gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+mix_prof(s,j)*gam_saved(t,j)*dexp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            mix_gam(t,i)=exp(-log(summ))
            mix_gam(t,i)=(mix_gam(t,i)+gam_saved(t,i))/2.0_8
            summ=0.0_8
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
   profile=0.0_8
   do t=1,3
      do i=0,50
         profile(t,i)=profil(t,i)/sum(profil)
      end do
   end do

   gam=1.0
   gam_saved=1.0
   summ=0.0_8
   not_conv=.true.
   
   do while (not_conv)
      gam_saved(:,:)=gam(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*dexp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam(t,i)=exp(-log(summ))
            gam(t,i)=(gam(t,i)+gam_saved(t,i))/2.0_8
            summ=0.0_8
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

   profile=0.0_8
   do t=1,3
      do i=0,50
         profile(t,i)=profil2(t,i)/sum(profil2)
      end do
   end do
   
   gam_sol=1.0
   gam_saved=1.0
   summ=0.0_8
   not_conv=.true.

   do while (not_conv)
      gam_saved(:,:)=gam_sol(:,:)
      do t=1,3
         do i=0,50
            do s=1,3
               do j=0,50
                  summ=summ+profile(s,j)*gam_saved(t,j)*dexp((-E_dd3((i*punit)-maxsig,(j*punit)-maxsig,t,s))/(R*SysTemp*Jtokcal))
               end do
            end do
            gam_sol(t,i)=exp(-log(summ))
            gam_sol(t,i)=(gam_sol(t,i)+gam_saved(t,i))/2.0_8
            summ=0.0_8
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
      L(i)=(coord/2.0_8)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_8)
   end do
   gammasg(1)=log(phi(1)/z(1))+(coord/2.0_8)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
      &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))
   
   gammasg(2)=log(phi(2)/z(2))+(coord/2.0_8)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
      &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))


   !! Additional dispersion correction
   A=omega*(0.5*(sac_disp(1)+sac_disp(2))-sqrt(sac_disp(1)*sac_disp(2)))
   gammadisp(1)=A*z(2)*z(2)
   gammadisp(2)=A*z(1)*z(1)

   write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
   gamma_solv=0.0_8
   gamma_sol=0.0_8
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

end subroutine sac_2013


subroutine sac2013_disp(nam,is_bonded,ident,elements,disp_con,sac_disp)
   use element_dict
   implicit none

   logical, dimension(:,:), allocatable, intent(in) :: is_bonded

   real(8), dimension(:), allocatable, intent(in) :: ident
   character(2), dimension(:), allocatable, intent(in) :: elements
   character(len=*) :: nam
   type(DICT_STRUCT), pointer :: disp_con

   type(DICT_DATA) :: disp

   character(2) :: symbol
   character(2) :: bond_string
   character(3) :: data_string

   real(8), intent(out) :: sac_disp

   real(8) :: atom_disp

   integer :: i,j, bond_count, disp_atoms

   sac_disp=0.0_8
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
      
