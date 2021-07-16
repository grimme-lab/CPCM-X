module crs

      contains


      subroutine calcgas(E_cosmo,id_scr,gas_chem,area,sv,su,pot,element,ident,disp_con, T,r_cav)
         use globals
         use element_dict
         real(8), intent(out) :: id_scr, gas_chem
         real(8), intent(in) :: T, E_cosmo
         real(8),dimension(:),allocatable, intent(in) :: area, sv, su, pot
         integer, allocatable, intent(in) :: ident(:)
         character(2), dimension(:), allocatable, intent(in) :: element
         type(DICT_STRUCT), pointer, intent(in) :: disp_con, r_cav
         !real(8), dimension(10) :: param
         type(DICT_DATA) :: disp!, r_c
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
         integer, dimension(:), allocatable, intent(in) :: ident
         character(2), dimension(:), allocatable, intent(in),optional :: ele2
         integer, dimension(:), allocatable, intent(in),optional :: id2

         ! LOCAL

         character(2), dimension(:), allocatable :: element2
         real(8), dimension(:), allocatable :: ident2
         integer, intent(in) :: atom1, atom2
         real(8) :: E_dd
         real(8) :: svdo, svac
         real(8) :: E_misfit, E_hb


         
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
         
         E_hb=0.0_8
         E_misfit=0.0_8
            E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
         E_misfit=(alpha/2)*(sv1+sv2)&
                 &*((sv1+sv2)+f_corr*(svt1+svt2))


            
            

               


         E_dd=E_hb+E_misfit

      end function E_dd

      subroutine iterate_solvent(pot_di,sv,svt,area,T,ident,element)
         use globals
         implicit none
         !real(8), dimension(10), intent(in) :: param
         real(8), dimension(:), allocatable, intent(in) :: sv, svt, area
         integer, allocatable, intent(in) :: ident(:)
         character(2), dimension(:), allocatable, intent(in) :: element
         real(8), dimension(:), allocatable, intent(inout) :: pot_di
         real(8), dimension(:), allocatable :: W_v 
         real(8), intent(in) :: T

         real(8) :: temppot, beta
         integer :: i, j

         allocate(W_v(size(sv)))
         W_v(:)=0.0_8
         beta=(R*Jtokcal*T)/param(7)
         temppot=0.0_8
         W_v=0.0_8
         !! For mixed solvent, mole fraction needs to be introduced in the following loop
         do j=1,size(pot_di)
            do i=1,size(pot_di)
               temppot=temppot+(area(i)*dexp((-E_dd&
                      &(param(5),param(3),param(4),param(6),sv(j),svt(j),sv(i),svt(i),ident,element,j,i)&
                      &/beta+pot_di(i))))
               W_v(j)=W_v(j)+area(i)
            end do
            !exit
            pot_di(j)=-log(temppot/W_v(j))
            temppot=0.0_8
         end do

      end subroutine iterate_solvent


      subroutine compute_solute(sol_pot,solv_pot,sv_sol,svt_sol,sv_solv,svt_solv,area_sol,area_solv,T,chem_pot_sol,&
                      &ident_sol,ident_solv,elem_sol,elem_solv)
         use globals
         implicit none
         !real(8), dimension(10), intent(in) :: param
         real(8), intent(out) :: chem_pot_sol
         real(8), dimension(:), allocatable, intent(in) :: sv_sol, svt_sol,sv_solv,svt_solv,area_solv,area_sol
         integer, allocatable, intent(in), dimension(:) :: ident_sol, ident_solv
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
         real(8), dimension(:), allocatable :: sv, svt, area
         integer, allocatable :: ident(:)
         real(8), dimension(:), allocatable, intent(inout) :: pot_di
         character(2), dimension(:), allocatable, intent(in) :: element
         real(8), intent(in) :: T
         real(4), intent(in) :: conv_crit
         integer, intent(in) :: max_cycle
         
         integer :: i, j
         real(8), dimension(:), allocatable :: saved_potdi
         logical :: not_conv
         
   
         not_conv=.TRUE.
         allocate(pot_di(size(sv)))
         allocate(saved_potdi(size(sv)))
         pot_di(:)=0.0_8
         saved_potdi(:)=0.0_8
         i=0
         do while (not_conv) 
            i=i+1
            saved_potdi(:)=pot_di(:)
            Call iterate_solvent(pot_di,sv,svt,area,T,ident,element)
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
         if (not_conv) then
            write(*,*) "Error, chemical potential could not be converged"
            error stop
         else
            write(*,*) "Done! Chemical potential converged after Cycle ",i
         end if


      end subroutine compute_solvent

end module crs
