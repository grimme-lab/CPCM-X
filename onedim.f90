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
  ! Write(*,*) E_dd1

end function E_dd1



subroutine onedim(profil,profil2,vcosmo)
   use globals
   implicit none

   real(8), dimension(:) :: profil,profil2

   real(8) :: gam(0:49),maxsig,punit,profile(0:49), gam_saved(0:49),gam_sol(0:49)
   real(8) :: gamma_solv, gamma_sol, summ
   real(8) :: T, VNORM, ANORM, RNORM, QNORM, vcosmo

   integer :: i,j,z
   logical :: not_conv

   VNORM=66.69_8
   ANORM=79.53_8
   vcosmo=vcosmo*BtoA**3
  ! write(*,*) vcosmo
   not_conv=.TRUE.

   !! Pure Activity Coefficient of Solvent
   profile(:)=0.0_8
  ! write(*,*) profil
   do i=0,49
      profile(i)=profil(i)/sum(profil)
   end do
   T=298.15_8
   maxsig=0.025_8
   punit=0.001_8
   gam(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8
   do while (not_conv)
      gam_saved(:)=gam(:)
      do i=0,size(profile)-1
         do j=0,size(profile)-1
            summ=summ+profile(j)*gam_saved(j)*dexp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_8*R*Jtokcal))
         end do
         gam(i)=exp(-log(summ))
         gam(i)=(gam(i)+gam_saved(i))/2.0_8
         summ=0.0_8
      end do
      not_conv=.false.
      do i=0,size(gam)-1
         if (abs((gam(i)-gam_saved(i))) .LT. 0.000001) then
            cycle
         else
            not_conv=.true.
            exit
         end if
      end do
   end do
  
   !! Pure Activity Coefficient of Solute
  ! write(*,*) profil2
   profile(:)=0.0_8
   do i=0,49
      profile(i)=profil2(i)/sum(profil2)
   end do
 !  write(*,*) profile
   T=298.15_8
   maxsig=0.025_8
   punit=0.001_8
   gam_sol(:)=1.0
   gam_saved(:)=1.0
   summ=0.0_8
   not_conv=.true.
   do while (not_conv)
      gam_saved(:)=gam_sol(:)
      do i=0,size(profile)-1
         do j=0,size(profile)-1
            summ=summ+profile(j)*gam_saved(j)*dexp((-E_dd1((i*punit)-maxsig,(j*punit)-maxsig))/(298.15_8*R*Jtokcal))
         end do
         gam_sol(i)=exp(-log(summ))
         gam_sol(i)=(gam_sol(i)+gam_saved(i))/2.0_8
         summ=0.0_8
      end do
      not_conv=.false.
      do i=0,size(gam_sol)-1
         if (abs((gam_sol(i)-gam_saved(i))) .LT. 0.000001) then
            cycle
         else
            not_conv=.true.
            exit
         end if
      end do
   end do

   !! Staverman-Guggenheim equation

     RNORM = VCOSMO/VNORM 
     QNORM = sum(profil2)/ANORM



   write(*,*) "One dimensional calculation done!"
  ! write(*,*) mue
   gamma_solv=0.0_8
   gamma_sol=0.0_8

 !  write(*,*) gam_sol
   do i=0,size(profil)-1
  !    write(*,*) mue(i)
      gamma_solv=gamma_solv+profil2(i)/7.5*(log(gam(i))-log(gam_sol(i)))
      gamma_sol=gamma_sol+profil2(i)/7.5*log(gam(i))
   end do
   write(*,*) gamma_solv, gamma_sol*R*T*Jtokcal



end subroutine onedim
      
