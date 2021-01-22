module sac_mod

   contains
function E_dd1(sigma1,sigma2)

   real(8), intent(in) :: sigma1,sigma2
   real(8), parameter :: EPS=3.667_8, aef= 7.5, e0=2.395E-4, c_hb=85580.0, s_hb=0.0084
   real(8) :: E_dd1,svdo,svac,E_misfit,E_hb, fpol, alpha, alphaprime

   fpol=(EPS-1.0_8)/(EPS+0.5_8)
   alpha=(0.3_8*aef**(1.5))/e0
   alphaprime=fpol*alpha

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



subroutine sac_2005(profil,profil2,vcosmo1,vcosmo2)
   use globals
   implicit none

   real(8), dimension(0:50) :: profil,profil2

  ! real(8), dimension(1:2,0:49) :: sigma_profiles

   real(8) :: gam(0:50),maxsig,punit,profile(0:50), gam_saved(0:50),gam_sol(0:50)
   real(8) :: gamma_solv, gamma_sol,gamma_test, summ, mix_prof(0:50), mix_gam(0:50)
   real(8) :: T, VNORM, ANORM, RNORM(2), QNORM(2), vcosmo1, z(2),vcosmo2
   real(8) :: Theta(2), Phi(2), L(2), coord, gammasg(2), bt, bp !SG Equation

   
  ! integer :: comp_num
   integer :: i,j
   logical :: not_conv

   z(1)=0.3_8
   z(2)=0.7_8
   
  ! comp_num=2

   VNORM=66.69_8
   ANORM=79.53_8
 !  vcosmo1=vcosmo1*BtoA**3
 !  vcosmo2=vcosmo2*BtoA**3
  ! write(*,*) vcosmo
   not_conv=.TRUE.

   !! Calculate the SIGMA Profile of the mixture

   do i=0,50
      mix_prof(i)=(z(1)*profil(i)+z(2)*profil2(i))/(z(1)*sum(profil)+z(2)*sum(profil2))
   end do
  ! write(*,*) E_dd1(-0.025_8,-0.016_8)
  ! write(*,*) mix_prof
  ! stop
   !! Mixed Activity Coefficient
   T=298.15_8
   maxsig=0.025_8
   punit=0.001_8
   mix_gam(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8

   !write(*,*) (1*punit)-maxsig
   
   do while (not_conv)
      gam_saved(:)=mix_gam(:)
      do i=0,50!size(mix_prof)-1
         do j=0,50!size(mix_prof)-1
        ! write(*,*) j*punit-maxsig, summ, mix_prof(j), gam_saved(j), E_dd1((i*punit)-maxsig,(j*punit)-maxsig)
            summ=summ+mix_prof(j)*gam_saved(j)*dexp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_8*R*Jtokcal))
         end do
        ! write(*,*) summ
        ! stop
         mix_gam(i)=exp(-log(summ))
         mix_gam(i)=(mix_gam(i)+gam_saved(i))/2.0_8
         summ=0.0_8
      end do
      not_conv=.false.
      do i=0,50!size(mix_gam)-1
         if (abs((mix_gam(i)-gam_saved(i))) .LT. 0.000001) then
            cycle
         else
            not_conv=.true.
            exit
         end if
      end do
   end do

   !write(*,*) mix_gam
   !stop
   !! Pure Activity Coefficient of 1
   profile(:)=0.0_8
  ! write(*,*) profil
   do i=0,50
      profile(i)=profil(i)/sum(profil)
   end do
   !write(*,*) profile

   gam(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8
   not_conv=.true. 
   do while (not_conv)
      gam_saved(:)=gam(:)
      do i=0,50
         do j=0,50
    !     write(*,*) summ
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

  ! write(*,*) gam

   !! Pure Activity Coefficient of 2
   !write(*,*) profil2
   profile(:)=0.0_8
  ! write(*,*) sum(profil2)
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

 !  write(*,*) gam_sol
   !! Staverman-Guggenheim equation
   coord=10

   RNORM(1) = VCOSMO1/VNORM 
   QNORM(1) = sum(profil)/ANORM
   RNORM(2) = VCOSMO2/VNORM
   QNORM(2) = sum(profil2)/ANORM

 !  write(*,*) RNORM(1), QNORM(1)
 !  write(*,*) RNORM(2), QNORM(2)
   bt=z(1)*QNORM(1)+z(2)*QNORM(2)
   bp=z(1)*RNORM(1)+z(2)*RNORM(2)
  
   
   do i=1,2
      Theta(i)=(z(i)*QNORM(i))/(bt)
      Phi(i)=(z(i)*RNORM(i))/(bp)
      L(i)=(coord/2.0_8)*(RNORM(i)-QNORM(i))-(RNORM(i)-1.0_8)
   !   write(*,*) Theta(i), Phi(i), L(i)
   end do
   gammasg(1)=log(phi(1)/z(1))+(coord/2.0_8)*QNORM(1)*log(Theta(1)/Phi(1))+L(1)-&
      &(Phi(1)/z(1))*(z(1)*L(1)+z(2)*L(2))
   
   gammasg(2)=log(phi(2)/z(2))+(coord/2.0_8)*QNORM(2)*log(Theta(2)/Phi(2))+L(2)-&
      &(Phi(2)/z(2))*(z(1)*L(1)+z(2)*L(2))
!     write(*,*) gammasg(1), gammasg(2)


   write(*,*) "COSMO-SAC Acitivity Coefficient Prediction:"
  ! write(*,*) mue
   gamma_solv=0.0_8
   gamma_sol=0.0_8
  ! gamma_test=0.0_8
 !  write(*,*) gam_sol
   do i=0,50
  !    write(*,*) mue(i)
      gamma_solv=gamma_solv+(profil(i)/7.5_8*(log(mix_gam(i)/gam(i))))
      gamma_sol=gamma_sol+(profil2(i)/7.5_8*log(mix_gam(i)/gam_sol(i)))
     ! gamma_test=gamma_test+(profil2(i)/7.5_8*log(gam_sol(i)))
   end do
 !  write(*,*) gamma_test*R*T*jtokcal
 !  write(*,*) gamma_sol*R*T*jtokcal
   gamma_solv=exp(gammasg(1)+gamma_solv)
   gamma_sol=exp(gammasg(2)+gamma_sol)
   write(*,*) "Results for Mixture with Compound 1 x= ",z(1)," and Compound 2 x= ",z(2),"."
   write(*,*) "Gamma(1)= ",gamma_solv, "Gamma(2)= ", gamma_sol
   write(*,*) "lnGamma(1)= ", log(gamma_solv),"lnGamma(2)= ", log(gamma_sol)



end subroutine sac_2005

end module sac_mod
      
