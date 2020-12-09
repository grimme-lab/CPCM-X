program COSMO
   use element_dict
   use globals
   implicit none
   integer :: i,j
   real(8), dimension(:), allocatable :: solute_su, solute_area, solute_sv, solute_sv0,solvent_pot,solute_pot
   real(8), dimension(:), allocatable :: solvent_su, solvent_area, solvent_sv, solvent_sv0, solute_svt, solvent_svt
   real(8), dimension(:), allocatable :: sol_pot, solv_pot
   real(8), dimension(:,:), allocatable :: solvent_xyz, solute_xyz
   character(2), dimension(:), allocatable :: solute_ident, solvent_ident
   real(8), dimension(10) :: param
   real(8) :: id_scr,gas_chem,chem_pot_sol, T, temp
   logical :: gas
   

   type(DICT_STRUCT), pointer :: r_cav, disp_con
   
   gas=.TRUE.

   T=298.15_8

   Call initialize_param(param,r_cav,disp_con)


   Call read_cosmo("h2o.cosmo",solvent_ident,solvent_xyz,solvent_su,solvent_area,solvent_pot)
   Call average_charge(param(1), solvent_xyz,solvent_su,solvent_area,solvent_sv)
   Call average_charge(param(2), solvent_xyz,solvent_su,solvent_area,solvent_sv0)
   Call ortho_charge(solvent_sv,solvent_sv0,solvent_svt)
   Call read_cosmo("meo.cosmo",solute_ident, solute_xyz, solute_su, solute_area,solute_pot)
   Call average_charge(param(1), solute_xyz, solute_su, solute_area, solute_sv)
   Call average_charge(param(2), solute_xyz, solute_su, solute_area, solute_sv0)
   Call ortho_charge(solute_sv,solute_sv0,solute_svt)
  ! write(*,*) solute_su
   if (gas) then
      Call calcgas(id_scr,gas_chem,solute_area,solute_sv,solute_su,solute_pot,solute_ident,disp_con,param, T,r_cav)
   end if

   Call compute_solvent(solv_pot,param,solvent_sv,solvent_svt,solvent_area,T,1,0.001)
   Call compute_solute(sol_pot,solv_pot,param,solute_sv,solute_svt,solvent_sv,&
         &solvent_svt,solute_area,solvent_area,T,chem_pot_sol)
  
 !     temp=0
 !     do j=1,size(solvent_sv)
 !        temp=temp+E_dd(param(5),param(3),param(4),param(6),solvent_sv(1),solvent_svt(1),solute_sv(j),solute_svt(j))
 !        write(*,*) temp
 !     end do
  
   
   
   contains

      subroutine read_cosmo(compound,ident,xyz,charges,area,pot)
         use globals
         implicit none
         character(9), intent(in) :: compound
         character(30) :: line 
         real(8), dimension(:), allocatable,intent(out) :: charges,area,pot
         real(8), dimension(:,:), allocatable, intent(out) :: xyz
         character(2), allocatable, dimension(:),intent(out) :: ident
         integer :: i, io_error, dummy1, dummy2, num
         real(8) :: dummy3, dummy4, dummy5
         real(8), dimension(:), allocatable :: dummy_ident
         character(2) :: element
         
         

         open(1,file=compound)
         io_error=0 
         do while (io_error .GE. 0)
            read(1,*,iostat=io_error) line
         end do
         read(line,*) num
         allocate(charges(num))
         allocate(ident(num))
         allocate(dummy_ident(num))
         allocate(area(num))
         allocate(xyz(num,3))
         allocate(pot(num))
         rewind(1)
         ident(:)="XX"
         do while (line .NE. "$segment_information")
            read(1,*) line
         end do
         
         num=1
         dummy4=0
         do while (.TRUE.)
            read(1,'(A1)',advance='no',iostat=io_error) line
            if (line=="#") then
               read(1,*)
               cycle
            else if (io_error .LT. 0) then
               exit
            else
               read(1,*) dummy1,dummy_ident(num),xyz(num,1),xyz(num,2),xyz(num,3),dummy3,area(num),charges(num),pot(num)
               pot(num)=pot(num)*BtoA 
               charges(num)=charges(num)*autokcal
               num=num+1
            end if
         end do
         rewind(1)
         do while (line .NE. "#atom")
            read(1,*) line
         end do

         do while (.TRUE.)
            read(1,'(A1)',advance='no',iostat=io_error) line
            if (line .NE. "$") then
               read(1,*) dummy1, dummy3, dummy4, dummy5, element
               do i=1,size(dummy_ident)
                  if (dummy1==dummy_ident(i)) then
                     ident(i)=element
                  end if
               end do
            else
               exit
            end if
         end do
         
         close(1)
         deallocate(dummy_ident)
      end subroutine read_cosmo

      function distance(xyz1,xyz2)
         use globals
         implicit none
         real(8), dimension(3), intent(in) :: xyz1, xyz2
         real(8) :: distance

         distance=BtoA*(sqrt((xyz1(1)-xyz2(1))**2.0_8+(xyz1(2)-xyz2(2))**2+(xyz1(3)-xyz2(3))**2.0_8))
         
      end function distance

      subroutine average_charge (r_av, xyz, charges, area,av_charge)
         implicit none
         real(8), dimension(:), allocatable, intent(in) :: charges, area
         real(8), dimension(:,:), allocatable, intent(in) :: xyz
         real(8), intent(in) :: r_av
         real(8), dimension(:), allocatable, intent(out) :: av_charge
         real(8) :: tmpcharge, tmpcounter, tmpdenominator, r_u2, r_av2
         real(8), parameter :: pi = 4*atan(1.0_8)
         integer :: num, i, j

         num = size(charges)
         allocate(av_charge(num))
         r_av2=r_av**2.0_8
         r_u2=0.0_8
         tmpcharge=0.0_8
         tmpcounter=0.0_8
         tmpdenominator=0.0_8
         do i=1,num
            do j=1,num
               !tmpcounter=tmpcounter+(charges(j)/(area(j)*2.0_8*pi*r_av)*&
                !         &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_8)/(r_av2))))
               r_u2=(area(j)/pi)
              ! tmpdenominator=tmpdenominator+(area(j)/(area(j)*2.0_8*pi*r_av)*&
               !          &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_8)/(r_av2))))
               tmpcounter=tmpcounter+(charges(j)*((r_u2*r_av2)/(r_u2+r_av2))*&
                         &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_8)/(r_u2+r_av2))))
               tmpdenominator=tmpdenominator+(((r_u2*r_av2)/(r_u2+r_av2))*&
                             &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_8)/(r_u2+r_av2))))
               r_u2=0.0_8
            end do
            tmpcharge=tmpcounter/tmpdenominator
            tmpcounter=0.0_8
            tmpdenominator=0.0_8
            av_charge(i)=tmpcharge
         end do
      end subroutine average_charge

      subroutine ortho_charge (v,v0,vt)
         implicit none
         real(8), dimension(:), allocatable, intent(in) :: v,v0
         real(8), dimension(:), allocatable, intent(out) :: vt

         integer :: i

         allocate(vt(size(v)))

         do i=1,size(v)
            vt(i)=v0(i)-0.816*v(i)
         end do

      end subroutine

      subroutine initialize_param(param,r_cav,disp_con)
         implicit none

         real(8), dimension(10) :: param
         type(DICT_STRUCT), pointer, intent(inout) :: r_cav, disp_con

         type(DICT_DATA) :: data1
         character(len=2) :: symbol
         
         !Setting global COSMO Parameters - TO DO: read from file

         param(1)=0.50_8
         param(2)=1.00_8
         param(3)=1288.0_8
         param(4)=2.40_8
         param(5)=7400.0_8
         param(6)=0.00820_8
         param(7)=7.10_8
         param(8)=0.140_8
         param(9)=-0.210_8
         param(10)=-9.150_8

         !Setting cavity radii - to do: read from file

         data1%param = 1.30_8
         call dict_create(r_cav, 'h', data1)
         data1%param = 2.00_8
         call dict_add_key(r_cav, 'c', data1)
         data1%param = 1.72_8
         call dict_add_key(r_cav, 'o', data1)

         !Setting dispersion contant -- need to be read from file

         data1%param = -0.041_8
         call dict_create(disp_con, 'h', data1)
         data1%param = -0.037_8
         call dict_add_key(disp_con, 'c', data1)
         data1%param = -0.042_8
         call dict_add_key(disp_con, 'o', data1)



      end subroutine initialize_param

      subroutine calcgas(id_scr,gas_chem,area,sv,su,pot,ident,disp_con,param, T,r_cav)
         use globals
         use element_dict
         real(8), intent(out) :: id_scr, gas_chem
         real(8), intent(in) :: T
         real(8),dimension(:),allocatable, intent(in) :: area, sv, su, pot
         character(2), dimension(:), allocatable, intent(in) :: ident
         type(DICT_STRUCT), pointer, intent(in) :: disp_con, r_cav
         real(8), dimension(10) :: param
         type(DICT_DATA) :: disp, r_c
         real(8) :: E_gas, E_solv,  dEreal, ediel, edielprime, vdW_gain, thermo, beta
         integer :: dummy1, ioerror, i 

         open(1,file="energy")
         read(1,*,iostat=ioerror)
         if (ioerror .NE. 0) then
            write(*,*) "Problem while reading energies (check energy file)."
            error stop
         else
            read(1,*) dummy1,E_gas
            read(1,*) dummy1,E_solv
         end if
         dEreal=(E_gas-E_solv)
         ediel=0
         edielprime=0
         do i=1,size(sv)
            ediel=ediel+(area(i)*pot(i)*su(i))
            edielprime=edielprime+(area(i)*pot(i)*sv(i))
         end do
         write(*,*) "Ediel: ", ediel/2 
         id_scr=(dEreal+(0.4*edielprime)-(0.4*ediel))
         dEreal=dEreal*autokcal
         write(*,*) id_scr, dEreal
         vdW_gain=0

         do i=1,size(area)
            disp=dict_get_key(disp_con, ident(i))
            vdW_gain=vdW_gain+(area(i)*(BtoA**2)*disp%param)
         end do
         thermo=param(10)*R*jtokcal*T
         write(*,*) vdW_gain, thermo, id_scr 

         beta=1!(k_b*Jtokcal*T*N_A)/param(7)
         !!! RING CORRECTION IS MISSING ATM
         gas_chem=beta*(-id_scr-vdW_gain-thermo)!-ring_corr
         write(*,*) gas_chem 
         

      end subroutine calcgas

      function E_dd(c_hb,alpha,f_corr,s_hb,sv1,svt1,sv2,svt2)
         implicit none
         real(8), intent(in) :: c_hb, alpha, f_corr,s_hb
         real(8), intent(in) :: sv1, svt1, sv2, svt2
        
         real(8) :: E_dd
         real(8) :: svdo, svac
         real(8) :: E_misfit, E_hb
         
         svac=0
         svdo=0

         if (sv1 .GE. sv2) then
            svac=sv1
            svdo=sv2
         else
            svac=sv2
            svdo=sv1
         end if
         E_hb=0.0_8
         E_misfit=0.0_8
      !   E_hb=c_hb*max(0.0_8,svac-s_hb)*min(0.0_8,svdo+s_hb)
      
         E_misfit=(alpha/2)*(sv1+svt1)&
                 &*((sv1+sv2)+f_corr*(svt1+svt2))
         write(*,*) alpha
         E_dd=E_hb+E_misfit
         !write(*,*) E_dd

      end function E_dd

      subroutine iterate_solvent(pot_di,param,sv,svt,area,T)
         use globals
         implicit none
         real(8), dimension(10), intent(in) :: param
         real(8), dimension(:), allocatable, intent(in) :: sv, svt, area
         real(8), dimension(:), allocatable, intent(inout) :: pot_di
         real(8), dimension(:), allocatable :: W_v 
         real(8), intent(in) :: T

         real(8) :: temppot, beta
         integer :: i, j

         allocate(W_v(size(sv)))
         W_v(:)=0.0_8
         beta=(k_b*Jtokcal*T*N_A)/param(7)
         temppot=0.0_8
         W_v=0.0_8
         !! For mixed solvent, mole fraction needs to be introduced in the following loop
         do j=1,size(pot_di)
            do i=1,size(sv)
               temppot=temppot+(area(i)*exp((-beta**(-1.0_8))*&
                      &E_dd(param(5),param(3),param(4),param(6),sv(j),svt(j),sv(i),svt(i))&
                      &+pot_di(i)))
              ! write(*,*) temppot
               W_v(j)=W_v(j)+area(i)
            end do
            pot_di(j)=-log(W_v(j)**(-1)*temppot)
            temppot=0.0_8
         end do

      end subroutine iterate_solvent


      subroutine compute_solute(sol_pot,solv_pot,param,sv_sol,svt_sol,sv_solv,svt_solv,area_sol,area_solv,T,chem_pot_sol)
         use globals
         implicit none
         real(8), dimension(10), intent(in) :: param
         real(8), intent(out) :: chem_pot_sol
         real(8), dimension(:), allocatable, intent(in) :: sv_sol, svt_sol,sv_solv,svt_solv,area_solv,area_sol
         real(8), dimension(:), allocatable, intent(inout) :: solv_pot,sol_pot
         real(8), dimension(:), allocatable :: W_v 
         real(8), intent(in) :: T

         real(8) :: temppot, beta
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
               temppot=temppot+(area_solv(i)*exp((-beta**(-1.0_8))*&
                      &E_dd(param(5),param(3),param(4),param(6),sv_sol(j),svt_sol(j),sv_solv(i),svt_solv(i))&
                      &+solv_pot(i)))
               W_v(j)=W_v(j)+area_solv(i)
            end do
            sol_pot(j)=-log(W_v(j)**(-1)*temppot)
            temppot=0.0_8
         end do
         temppot=0.0_8
         !write(*,*) sol_pot
         do i=1,size(sol_pot)
            temppot=temppot+(area_sol(i)*sol_pot(i))
         end do
         !write(*,*) temppot
         !chem_pot_sol=sum(sol_pot)*beta
         chem_pot_sol=beta*temppot-(param(8)*R*Jtokcal*T*log(sum(area_solv)))
         write(*,*) chem_pot_sol
      end subroutine compute_solute

      subroutine compute_solvent(pot_di,param,sv,svt,area,T,max_cycle,conv_crit)
         use globals
         implicit none
         real(8), dimension(10), intent(in) :: param
         real(8), dimension(:), allocatable, intent(in) :: sv, svt, area
         real(8), dimension(:), allocatable, intent(inout) :: pot_di
         real(8), intent(in) :: T
         real(4), intent(in) :: conv_crit
         integer, intent(in) :: max_cycle
         
         real(8) :: dummy
         integer :: i, j
         real(8), dimension(:), allocatable :: saved_potdi
         logical :: not_conv
         
         not_conv=.TRUE.
         allocate(pot_di(size(sv)))
         allocate(saved_potdi(size(sv)))
         pot_di(:)=0.0_8
         i=0
         do while (not_conv) 
            i=i+1
            write(*,*) "Cycle ",i
            saved_potdi(:)=pot_di(:)
            write(*,*) sum(pot_di)
            Call iterate_solvent(pot_di,param,sv,svt,area,T)
            not_conv=.FALSE.
            do j=1,size(pot_di)
               if (abs(saved_potdi(j)-pot_di(j)) .LE. (conv_crit)) then
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



end program COSMO

