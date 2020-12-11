function E_dd1(sigma1,sigma2)

   real(8), intent(in) :: sigma1,sigma2
   real(8), parameter :: alpha = 11609, c_hb=7400, s_hb=0.0082
   real(8) :: E_dd1,svdo,svac,E_misfit,E_hb

   svac=0
   svdo=0

   if (sigma1 .GT. sigma2) then
      svac=sigma1
      svdo=sigma2
   else
      svac=sigma2
      svdo=sigma1
   end if
   E_hb=0.0_8
   E_misfit=0.0_8
  ! E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
  ! write(*,*) sigma1,sigma2
   E_misfit=alpha/2*(sigma1+sigma2)**2.0_8
   E_dd1=E_misfit+E_hb


end function E_dd1



subroutine onedim(profil,profil2)
   implicit none

   real(8), dimension(:) :: profil,profil2

   real(8) :: mue(0:49),maxsig,punit,profile(0:49), mue_saved(0:49)
   real(8) :: chem_pot_solv, chem_pot_sol

   integer :: i,j,z
   profile(:)=0.0_8
   do i=0,49
      profile(i)=profil(i)/sum(profil)
   end do
   maxsig=0.025_8
   punit=0.001_8

   mue(:)=0
   mue_saved(:)=0
   do z=1,105
   do i=0,size(profile)
      do j=0,size(profile)
       mue_saved(i)=mue_saved(i)+profile(j)*dexp(-E_dd1((i*punit)-maxsig,(j*punit)-maxsig)+mue(j))
      end do
       mue_saved(i)=-log(mue_saved(i))
   end do

      mue(:)=mue_saved(:)
      mue_saved(:)=0.0_8
 !  do i=1,size(profile)
 !     write(*,*)  mue(i)
 !  end do
 !  write(*,*)
   end do
  
   write(*,*) "One dimensional calculation done!"
  ! write(*,*) mue
   chem_pot_solv=0.0_8
   chem_pot_sol=0.0_8
   do i=0,size(profil)
  !    write(*,*) mue(i)
      chem_pot_solv=chem_pot_solv+profil(i)*mue(i)
      chem_pot_sol=chem_pot_sol+profil2(i)*mue(i)
   end do
   write(*,*) chem_pot_solv*0.0832-0.52, chem_pot_sol*0.0832-0.52



end subroutine onedim
      
