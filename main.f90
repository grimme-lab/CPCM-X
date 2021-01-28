program COSMO
   use element_dict
   use globals
   use sort
   use initialize_cosmo
   use sigma_av
   use sac_mod
   use bonding
   use profile
   implicit none
   integer :: i,j,z
   real(8), dimension(:), allocatable :: solute_su, solute_area, solute_sv, solute_sv0,solvent_pot,solute_pot
   real(8), dimension(:), allocatable :: solvent_su, solvent_area, solvent_sv, solvent_sv0, solute_svt, solvent_svt
   real(8), dimension(:), allocatable :: sol_pot, solv_pot, solvent_ident, solute_ident
   real(8), dimension(:,:), allocatable :: solvent_xyz, solute_xyz, solvat_xyz, solat_xyz
   character(2), dimension(:), allocatable :: solute_elements, solvent_elements, solute_hb, solvent_hb
   logical, dimension(:,:), allocatable :: solute_bonds, solvent_bonds
   real(8), dimension(3,0:50) :: solvent_sigma3, solute_sigma3
   character(20) :: solvent, solute
   !real(8), dimension(10) :: param
   real(8), dimension(5) :: T_a
   real(8) :: id_scr,gas_chem,chem_pot_sol, temp, temp2, T, solute_volume, solvent_volume,&
      &solute_energy, solvent_energy, solvent_sigma(0:50), solute_sigma(0:50),sac_disp(2)
   logical :: gas,sig_in,sac
  

   type(DICT_STRUCT), pointer :: r_cav, disp_con
  
   gas=.TRUE.
   
   !! Read Command Line Arguments and set Parameters accordingly

   Call getargs(solvent,solute,T,sig_in)
   Call initialize_param(r_cav,disp_con)
   !! Read Sigma Profiles (--sigma) or create Sigma Profiles from .cosmo files (default)
   T=SysTemp
   if (sig_in) then

      write(*,*) "Reading Sigma Profile"
      Call read_singlesig(solvent_sigma,trim(solvent)//".sigma",solvent_volume)
      Call read_singlesig(solute_sigma,trim(solute)//".sigma",solute_volume)
      
      Call read_triplesig(solvent_sigma3,trim(solvent)//".sigma",solvent_volume)
      Call read_triplesig(solute_sigma3,trim(solute)//".sigma",solute_volume)
   else

      write(*,*) "Creating Sigma Profile from COSMO data"

      Call read_cosmo(trim(solvent)//".cosmo",solvent_elements,solvent_ident,solvent_xyz,solvent_su,&
          &solvent_area,solvent_pot,solvent_volume,solvent_energy,solvat_xyz)
      Call read_cosmo(trim(solute)//".cosmo",solute_elements,solute_ident, solute_xyz, solute_su,&
         &solute_area,solute_pot,solute_volume,solute_energy,solat_xyz)
  
      Call average_charge(param(1), solvent_xyz,solvent_su,solvent_area,solvent_sv)
      Call average_charge(param(1), solute_xyz, solute_su, solute_area, solute_sv)


      Call det_bonds(solute_ident,solat_xyz,solute_elements,solute_bonds)
      Call hb_grouping(solute_ident,solute_elements,solute_bonds,solute_hb)
      Call det_bonds(solvent_ident,solvat_xyz,solvent_elements,solvent_bonds)
      Call hb_grouping(solvent_ident,solvent_elements,solvent_bonds,solvent_hb)

      Call single_sigma(solvent_sv,solvent_area,solvent_sigma,trim(solvent))
      Call single_sigma(solute_sv,solute_area,solute_sigma,trim(solute))

      Call split_sigma(solvent_sv,solvent_area,solvent_hb,solvent_ident,solvent_elements,&
         &solvent_sigma3,trim(solvent))
      Call split_sigma(solute_sv,solute_area,solute_hb,solute_ident,solute_elements,&
         &solute_sigma3,trim(solute))
   end if
   
   select case (trim(model))
      case ("sac")! Do a COSMO-SAC calculation instead COSMO-RS
         Call sac_2005(solvent_sigma,solute_sigma,solvent_volume,solute_volume)
         stop
      case("sac2010")
         Call sac_2010(solvent_sigma3,solute_sigma3,solvent_volume,solute_volume)

      case("sac2013")
         Call sac2013_disp(trim(solvent),solvent_bonds,solvent_ident,solvent_elements,disp_con,sac_disp(1))
         Call sac2013_disp(trim(solute),solute_bonds,solute_ident,solute_elements,disp_con,sac_disp(2))
         Call sac_2013(solvent_sigma3,solute_sigma3,solvent_volume,solute_volume,sac_disp)
      case ("crs")
   

         !! COSMO-RS calculation starts here !!

         ! Calculate sv0,svt for COSMO-RS

         Call average_charge(param(2), solvent_xyz,solvent_su,solvent_area,solvent_sv0)
         Call ortho_charge(solvent_sv,solvent_sv0,solvent_svt)
      
         Call average_charge(param(2), solute_xyz, solute_su, solute_area, solute_sv0)
         Call ortho_charge(solute_sv,solute_sv0,solute_svt)

         ! Calcualtion of Gas Phase energies

         if (gas) then
            Call calcgas(solute_energy,id_scr,gas_chem,solute_area,solute_sv,solute_su,&
               &solute_pot,solute_elements,solute_ident,disp_con, T,r_cav)
         end if

         ! Computation of COSMO-RS equations (here may be something wrong atm)

         Call compute_solvent(solv_pot,solvent_sv,solvent_svt,solvent_area,T,500,0.0001,solvent_ident,solvent_hb)
         Call compute_solute(sol_pot,solv_pot,solute_sv,solute_svt,solvent_sv,&
         &solvent_svt,solute_area,solvent_area,T,chem_pot_sol,solute_ident,solvent_ident,solute_elements,solvent_hb)
  
      write(*,*) "calc_gas_chem: ", gas_chem
      write(*,*) "calc_sol_chem: ", chem_pot_sol
      write(*,*) "G_solvshift: ", chem_pot_sol-gas_chem-4.28!-R*T*Jtokcal*log((solute_volume*(BtoA**3.0_8)*N_a*1000_8*10E-30)/22.414)

      deallocate(solute_su,solute_sv,solute_svt,solute_sv0,solvent_su,solvent_sv,&
         &solvent_svt,solvent_sv0,solvent_area,solute_area,solvent_xyz,solute_xyz,&
         &solv_pot,sol_pot)
      end select
   
   
   contains


      subroutine calcgas(E_cosmo,id_scr,gas_chem,area,sv,su,pot,element,ident,disp_con, T,r_cav)
         use globals
         use element_dict
         real(8), intent(out) :: id_scr, gas_chem
         real(8), intent(in) :: T, E_cosmo
         real(8),dimension(:),allocatable, intent(in) :: area, sv, su, pot,ident
         character(2), dimension(:), allocatable, intent(in) :: element
         type(DICT_STRUCT), pointer, intent(in) :: disp_con, r_cav
         !real(8), dimension(10) :: param
         type(DICT_DATA) :: disp, r_c
         real(8) :: E_gas, dEreal, ediel, edielprime, vdW_gain, thermo, beta, avcorr
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
         avcorr=(edielprime-ediel)/2.0_8*0.8_8
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
         write(*,*) gas_chem
         

      end subroutine calcgas


      function E_dd(c_hb,alpha,f_corr,s_hb,sv1,svt1,sv2,svt2,ident,element,atom1,atom2,id2,ele2)
         implicit none
         real(8), intent(in) :: c_hb, alpha, f_corr,s_hb
         real(8), intent(in) :: sv1, svt1, sv2, svt2
         character(2), dimension(:), allocatable, intent(in) :: element
         real(8), dimension(:), allocatable, intent(in) :: ident
         character(2), dimension(:), allocatable, intent(in),optional :: ele2
         real(8), dimension(:), allocatable, intent(in),optional :: id2

         ! LOCAL

         character(2), dimension(:), allocatable :: element2
         real(8), dimension(:), allocatable :: ident2
         integer, intent(in) :: atom1, atom2
         integer :: nl1,nl2,nr1,nr2
         real(8) :: E_dd
         real(8) :: svdo, svac
         real(8) :: E_misfit, E_hb

         logical :: hbnof

       !  hbnof=.true.
         
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

         !! Set Neighbors, not sufficient information in COSMO files

      !   nl1=atom1-1
      !   if (nl1 .LT. 1) nl1=size(ident)
      !   nl2=atom2-1
      !   if (nl2 .LT. 1) nl2=size(ident2)
      !   nr1=atom1+1
      !   if (nr1 .GT. size(ident)) nr1=1
      !   nr2=atom2+1
      !   if (nr2 .GT. size(ident2)) nr2=1


         !! Start E_dd Calculation
         
         E_hb=0.0_8
         E_misfit=0.0_8
        ! if ((element(int(ident(atom1))) .EQ. 'OH') .AND. (element(int(ident(atom2))) .EQ. 'OH')) then
            E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
        ! end if
         E_misfit=(alpha/2)*(sv1+sv2)&
                 &*((sv1+sv2)+f_corr*(svt1+svt2))


 !        if (element(int(ident(atom1))) .EQ. 'h') then
 !           if (element2(int(ident2(atom2))) .EQ. 'o') then
 !              E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
 !           end if
 !        else if (element(int(ident(atom1))) .EQ. 'o') then
 !           if (element2(int(ident2(atom2))) .EQ. 'h') then
 !              E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
 !           end if
 !        end if

         !! HBNOF does not work, because provided information is not sufficient

 !        if (hbnof) then
 !           if (E_hb .NE. 0.0_8) then
 !              select case (element(int(ident(atom1))))
 !                 case default
 !                    E_hb=0.0_8
 !                 case ('h')
 !                    if (((element(int(ident(nr1))) .EQ. 'o') .OR. ((element(int(ident(nl1))) .EQ. 'o')))) then
 !                    else
 !                       E_hb=0.0_8
 !                    end if
 !                 case ('o')
 !                    if (element2(int(ident2(atom2))) .EQ. 'h') then
 !                       if ((element2(int(ident2(nr2))) .EQ. 'o') .OR. ((element2(int(ident2(nl2)))) .EQ. 'o')) then
 !                       !   write(*,*) "yes"
 !                       else
 !                          E_hb=0.0_8
 !                       end if
 !                    else
 !                       E_hb=0.0_8
 !                    end if  
 !              end select
 !
 !           end if
 !       end if
            
            

               


         E_dd=E_hb+E_misfit

      end function E_dd

      subroutine iterate_solvent(pot_di,sv,svt,area,T,ident,element)
         use globals
         implicit none
         !real(8), dimension(10), intent(in) :: param
         real(8), dimension(:), allocatable, intent(in) :: sv, svt, area,ident
         character(2), dimension(:), allocatable, intent(in) :: element
         real(8), dimension(:), allocatable, intent(inout) :: pot_di
         real(8), dimension(:), allocatable :: W_v 
         real(8), intent(in) :: T

         real(8), dimension(:), allocatable :: test_di
         real(8) :: temppot, beta
         integer :: i, j

         allocate(W_v(size(sv)))
         W_v(:)=0.0_8
    !     allocate(test_di(size(pot_di)))
      !   test_di(:)=pot_di(:)
         beta=(R*Jtokcal*T)/param(7)
         temppot=0.0_8
         W_v=0.0_8
         !! For mixed solvent, mole fraction needs to be introduced in the following loop
         do j=1,size(pot_di)
            do i=1,size(pot_di)
       !        if (i .NE. j) then
               temppot=temppot+(area(i)*dexp((-E_dd&
                      &(param(5),param(3),param(4),param(6),sv(j),svt(j),sv(i),svt(i),ident,element,j,i)&
                      &/beta+pot_di(i))))
               W_v(j)=W_v(j)+area(i)
        !       end if
            end do
            !exit
            pot_di(j)=-log(temppot/W_v(j))
            temppot=0.0_8
         end do
   !      deallocate(test_di)

      end subroutine iterate_solvent


      subroutine compute_solute(sol_pot,solv_pot,sv_sol,svt_sol,sv_solv,svt_solv,area_sol,area_solv,T,chem_pot_sol,&
                      &ident_sol,ident_solv,elem_sol,elem_solv)
         use globals
         implicit none
         !real(8), dimension(10), intent(in) :: param
         real(8), intent(out) :: chem_pot_sol
         real(8), dimension(:), allocatable, intent(in) :: sv_sol, svt_sol,sv_solv,svt_solv,area_solv,area_sol,ident_sol,ident_solv
         real(8), dimension(:), allocatable, intent(inout) :: solv_pot,sol_pot
         character(2), dimension(:), allocatable, intent(in) :: elem_sol, elem_solv
         real(8), dimension(:), allocatable :: W_v 
         real(8), intent(in) :: T

         real(8) :: temppot, beta,temp2
         integer :: i, j

         allocate(W_v(size(sv_sol)))
         allocate(sol_pot(size(sv_sol)))
         W_v(:)=0.0_8
         beta=(R*Jtokcal*T)/param(7)
         temppot=0.0_8
         sol_pot(:)=0
         !! For mixed solvent, mole fraction needs to be introduced in the following loop
         do j=1,size(sol_pot)
            do i=1,size(solv_pot)
             !  if (i .NE. j) then
               temppot=temppot+(area_solv(i)*dexp((-E_dd&
                      &(param(5),param(3),param(4),param(6),sv_sol(j),svt_sol(j),sv_solv(i),svt_solv(i),ident_sol,elem_sol,&
                      j,i,ident_solv,elem_solv)&
                      &/beta)+solv_pot(i)))
             !  end if
               W_v(j)=W_v(j)+area_solv(i)
            end do
            sol_pot(j)=-log(sum(area_solv)**(-1)*temppot)
            temppot=0.0_8
         end do
         temppot=0.0_8
         temp2=0
         do i=1,size(sol_pot)
            temppot=temppot+(area_sol(i)*sol_pot(i))
        !    write(*,*) area_sol(i), sol_pot(i), area_sol(i)*sol_pot(i)
         end do

         chem_pot_sol=temppot*beta-(param(8)*R*Jtokcal*T*log(sum(area_solv)))
         write(*,*) chem_pot_sol
         temppot=0
         do i=1,size(solv_pot)
            temppot=temppot+(area_solv(i)*solv_pot(i))
         end do
         write(*,*) beta*temppot
      end subroutine compute_solute

      subroutine compute_solvent(pot_di,sv,svt,area,T,max_cycle,conv_crit,ident,element)
         use globals
         implicit none
         !real(8), dimension(10), intent(in) :: param
         real(8), dimension(:), allocatable :: sv, svt, area, ident
         real(8), dimension(:), allocatable, intent(inout) :: pot_di
         character(2), dimension(:), allocatable, intent(in) :: element
         real(8), intent(in) :: T
         real(4), intent(in) :: conv_crit
         integer, intent(in) :: max_cycle
         
         real(8) :: dummy
         integer :: i, j, z
         real(8), dimension(:), allocatable :: saved_potdi,av
         logical :: not_conv
         
         !sv=anint(sv*1000)/1000
         !svt=anint(svt*1000)/1000
   
         not_conv=.TRUE.
         allocate(pot_di(size(sv)))
         allocate(saved_potdi(size(sv)))
       !  allocate(av(size(sv)))
         pot_di(:)=0.0_8
         saved_potdi(:)=0.0_8
       !  av(:)=0.0_8
         i=0
         do while (not_conv) 
            i=i+1
        !    write(*,*) pot_di
      !      do z=1,size(pot_di)
      !         av(z)=(saved_potdi(z)+pot_di(z))/2.0_8
      !      end do
            saved_potdi(:)=pot_di(:)
!            write(*,*) sum(pot_di)
            Call iterate_solvent(pot_di,sv,svt,area,T,ident,element)
       !     pot_di(:)=av(:)
            not_conv=.FALSE.
            do j=1,size(pot_di)
               if (abs(saved_potdi(j)-pot_di(j)) .LE. (conv_crit/autokcal)) then
                  cycle
               else
                  not_conv=.TRUE.
                  exit
               end if
            end do
            if (i .GE. max_cycle) then
               exit
            end if
         end do
       !  write(*,*) pot_di(:)
       !  Call reorder(sv,pot_di,area)
         if (not_conv) then
            write(*,*) "Error, chemical potential could not be converged"
            error stop
         else
            write(*,*) "Done! Chemical potential converged after Cycle ",i
         end if


      end subroutine compute_solvent

      subroutine sigma_profile(sv,area,sigma,nam)

         real(8), dimension(:), allocatable,intent(in) :: sv,area

         real(8), dimension(:),intent(out) :: sigma(0:50)

         character(len=*), intent(in) :: nam

         integer :: sigma_min, sigma_max, i, j,tmp

         real(4) :: punit

         real(4), parameter :: sig_width=0.025_8
         integer, parameter :: n_sig=51
         real(4) :: counter(0:n_sig-1)

         real(8) :: profile(0:n_sig-1), chdval(0:n_sig-1), temp

         punit=0.001


         profile(:)=0.0_8
         counter(:)=0.0
         
         do i=0,n_sig-1
            profile(i) = 0.0_8
            counter(i) = -sig_width+punit*i
         end do
         do i= 1, size(sv)
            temp = sv(i)
           
            tmp = int((temp-counter(0))/punit)
 
            if (tmp<0) tmp=0
            if (tmp>n_sig-1) tmp=n_sig-1
            profile(tmp) = profile(tmp)+area(i)*(counter(tmp+1)-temp)/punit
            profile(tmp+1) = profile(tmp+1)+area(i)*(temp-counter(tmp))/punit
         end do
         open(unit=2,file=nam//"_sigma.txt",action="write",status="replace")
         
        ! write(*,*) profile
         do i=0,size(profile)-1
            write(2,*) counter(i),";", profile(i)!/sum(area)
         end do
         close(2)

         
         sigma(:)=profile(:)

      end subroutine

      subroutine read_sigma(sigma,nam,volume)

         character(len=*), intent(in) :: nam
         real(8), intent(out) :: volume
         real(8), dimension(:), intent(out) :: sigma (0:50)

         character(len=20) :: dummy1, dummy2
         integer :: io_error,i
         real(8) :: dummy3

         open(unit=2,file=nam)

         read(2,*)
         read(2,*)
         read(2,*) 
         read(2,*) dummy1,dummy2,volume
     	 io_error=0
         sigma(:) = 0.0_8
         i=0
     	 do i=0,50
     	    read(2,*,iostat=io_error) dummy3, sigma(i)
          !  write(*,*) i, dummy3, sigma(i)
     	 end do
     !	 write(*,*) sigma
         close(2)
     	
      end subroutine read_sigma


end program COSMO

