program COSMO
   use element_dict
   use globals
   use sort
   use initialize_cosmo
   use sigma_av
   use sac_mod
   use bonding
   use profile
   use pr
   use crs
   use mctc_env, only : wp, get_argument, fatal_error, error_type
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use sdm
   implicit none
   integer :: oh_sol, nh_sol, near_sol
   real(8), dimension(:), allocatable :: solute_su, solute_area, solute_sv, solute_sv0,solvent_pot,solute_pot
   real(8), dimension(:), allocatable :: solvent_su, solvent_area, solvent_sv, solvent_sv0, solute_svt, solvent_svt
   real(8), dimension(:), allocatable :: sol_pot, solv_pot
   real(8), dimension(:,:), allocatable :: solvent_xyz, solute_xyz, solvat_xyz, solat_xyz, solat2
   character(2), dimension(:), allocatable :: solute_elements, solvent_elements, solute_hb, solvent_hb
   logical, dimension(:,:), allocatable :: solute_bonds, solvent_bonds
   logical, dimension(:), allocatable :: solute_rings
   real(8), dimension(3,0:50) :: solvent_sigma3, solute_sigma3
   character(20) :: solvent, solute
   !real(8), dimension(10) :: param
   real(8) :: id_scr,gas_chem,chem_pot_sol, T, solute_volume, solvent_volume,&
      &solute_energy, solvent_energy, solvent_sigma(0:50), solute_sigma(0:50),sac_disp(2)
   logical :: gas,sig_in
   integer :: sol_nat, i
   integer, allocatable :: int_ident(:),solute_ident(:),solvent_ident(:)
   real(wp), allocatable :: surface(:), dsdr(:,:,:)

   type :: configuration
      character(len=:), allocatable :: input
      character(len=:), allocatable :: smd_solvent
      character(len=:), allocatable :: csm_solvent
      character(len=:), allocatable :: csm_solute
      real(wp) :: T
      real(wp) :: probe
      real(wp) :: z1,z2
      character(len=:), allocatable :: sac_param_path
      character(len=:), allocatable :: smd_param_path
      logical :: ML, sig_in, prof, smd_default
      character(len=:), allocatable :: model
   end type configuration

   type(configuration) :: config
   type(error_type), allocatable :: error

  

   type(DICT_STRUCT), pointer :: r_cav, disp_con
  
   gas=.TRUE.
   !! ------------------------------------------------------------ 
   !! Read Command Line Arguments and set Parameters accordingly
   !! ------------------------------------------------------------

   Call get_arguments(config,error)
   Call initialize_param(config%sac_param_path,config%model,r_cav,disp_con,config%csm_solvent)
   
   if (config%ML) then
      Call init_pr
      write(*,*) "Machine Learning Mode selected. Will Only Write an ML.data file." !! ML Mode deprecated
      ML=.TRUE.
   end if

   !! ----------------------------------------------------------------------------------
   !! Read Sigma Profiles (--sigma) - Not the default case
   !! ----------------------------------------------------------------------------------
   T=SysTemp
   if (config%sig_in) then

      write(*,*) "Reading Sigma Profile"
      Call read_singlesig(solvent_sigma,config%csm_solvent,solvent_volume)
      Call read_singlesig(solute_sigma,config%csm_solute,solute_volume)
      
      Call read_triplesig(solvent_sigma3,config%csm_solvent,solvent_volume)
      Call read_triplesig(solute_sigma3,config%csm_solute,solute_volume)
   else
   !! ----------------------------------------------------------------------------------
   !! Create the Sigma Profile from COSMO files
   !! ----------------------------------------------------------------------------------

      write(*,*) "Creating Sigma Profile from COSMO data"

   !! ------------------------------------------------------------------------------------
   !! Read necessary COSMO Data
   !! ------------------------------------------------------------------------------------
      Call read_cosmo(config%csm_solvent,solvent_elements,solvent_ident,solvent_xyz,solvent_su,&
          &solvent_area,solvent_pot,solvent_volume,solvent_energy,solvat_xyz)
      Call read_cosmo(config%csm_solute,solute_elements,solute_ident, solute_xyz, solute_su,&
         &solute_area,solute_pot,solute_volume,solute_energy,solat_xyz)
 
   !! ------------------------------------------------------------------------------------
   !! Sigma Charge Averaging and creating of a single Sigma Profile for Solute and Solvent
   !! ------------------------------------------------------------------------------------

      Call average_charge(param(1), solvent_xyz,solvent_su,solvent_area,solvent_sv)
      Call average_charge(param(1), solute_xyz, solute_su, solute_area, solute_sv)
      Call single_sigma(solvent_sv,solvent_area,solvent_sigma,"solvent")
      Call single_sigma(solute_sv,solute_area,solute_sigma,"solute")

   !! ------------------------------------------------------------------------------------
   !! Determination of HB Grouping and marking of Atom that are able to form HBs.
   !! Determination of Atoms in Rings, necessary for the PR2018 EOS (only ML Model)
   !! ------------------------------------------------------------------------------------
   if ((config%ML) .OR. (.NOT. (config%model .EQ. "sac"))) then
      Call det_bonds(solute_ident,solat_xyz,solute_elements,solute_bonds,oh_sol,nh_sol)
      Call hb_grouping(solute_ident,solute_elements,solute_bonds,solute_hb)
      Call det_bonds(solvent_ident,solvat_xyz,solvent_elements,solvent_bonds)
      Call hb_grouping(solvent_ident,solvent_elements,solvent_bonds,solvent_hb)
      
      if (config%ML) Call det_rings(solute_ident,solute_bonds,solute_rings,near_sol)
   end if

   !! ------------------------------------------------------------------------------------
   !! Creation of a splitted Sigma Profile, necessary for sac2010/sac2013
   !! ------------------------------------------------------------------------------------

      if (.NOT. (config%model .EQ. "sac")) then
      Call split_sigma(solvent_sv,solvent_area,solvent_hb,solvent_ident,solvent_elements,&
            &solvent_sigma3,"solvent")
      Call split_sigma(solute_sv,solute_area,solute_hb,solute_ident,solute_elements,&
            &solute_sigma3,"solute")
      end if

   !! ------------------------------------------------------------------------------------
   !! Exit here if you only want Sigma Profiles to be created 
   !! ------------------------------------------------------------------------------------
      if (onlyprof) then;
         write(*,*) "Only Profile mode choosen, exiting."
         stop
      end if
   end if
   

   !! ------------------------------------------------------------------------------------
   !! Choice of the different post COSMO Models (sac,sac2010,sac2013,COSMO-RS)
   !! ------------------------------------------------------------------------------------

   select case (trim(config%model))
      case ("sac")
         !Calculation of the Gas Phase (ideal gas --> ideal conductor)
         Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)
         !Calculation of the Solvent Phase (ideal conductor --> real solution)
         Call sac_2005(solvent_sigma,solute_sigma,solvent_volume,solute_volume,config%z1,config%z2)
         !Calculation of NES contributions (real gas --> ideal gas?)
         if (config%ML) Call pr2018(solute_area,solute_elements,solute_ident,oh_sol,nh_sol,near_sol)

         allocate (int_ident(maxval(solute_ident)))
         do i=1,maxval(solute_ident)
            int_ident(i)=i
         end do
   
         Call calculate_cds(int_ident,solute_elements,solat_xyz,config%probe,&
         &config%smd_solvent,config%smd_param_path,config%smd_default)

      case("sac2010")
         
         Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)
         Call sac_2010(solvent_sigma3,solute_sigma3,solvent_volume,solute_volume)
         if (config%ML) Call pr2018(solute_area,solute_elements,solute_ident,oh_sol,nh_sol,near_sol)

         allocate (int_ident(maxval(solute_ident)))
         do i=1,maxval(solute_ident)
            int_ident(i)=i
         end do
   
         Call calculate_cds(int_ident,solute_elements,solat_xyz,config%probe,&
         &config%smd_solvent,config%smd_param_path,config%smd_default)
   !! ------------------------------------------------------------------------------------
   !! The SAC 2013 Routine is not fully implemented and not supported anymore
   !! ------------------------------------------------------------------------------------
    !  case("sac2013")

     !    Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)
     !    Call sac2013_disp(trim(solvent),solvent_bonds,solvent_ident,solvent_elements,disp_con,sac_disp(1))
     !    Call sac2013_disp(trim(solute),solute_bonds,solute_ident,solute_elements,disp_con,sac_disp(2))
     !    Call sac_2013(solvent_sigma3,solute_sigma3,solvent_volume,solute_volume,sac_disp)
     !    Call pr2018(solute_area,solute_elements,solute_ident,oh_sol,nh_sol,near_sol)
      case ("crs")
   

         !! COSMO-RS calculation starts here !!

         ! Calculate sv0,svt for COSMO-RS

         Call average_charge(param(1)*2.0_wp, solvent_xyz,solvent_su,solvent_area,solvent_sv0)
         Call ortho_charge(solvent_sv,solvent_sv0,solvent_svt)
      
         Call average_charge(param(1)*2.0_wp, solute_xyz, solute_su, solute_area, solute_sv0)
         Call ortho_charge(solute_sv,solute_sv0,solute_svt)

         ! Calculation of Gas Phase energies

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
         deallocate(solute_su,solute_sv,solute_sv0,solvent_su,solvent_sv,&
         &solvent_sv0,solvent_area,solute_area,solvent_xyz,solute_xyz,&
         &solv_pot,sol_pot)
         stop
      end select
      if (config%ML) then
         write(*,*) "Writing ML data in ML.data"
         Call System("paste --delimiters='' ML.energy ML.gamma ML.pr > ML.data")
         Call System ("rm ML.energy ML.pr")
      else if (config%model .NE. "crs") then
         write(*,*) "Free Energy contributions:"
         write(*,*) "Ideal State (dG_is):", dG_is
         write(*,*) "Averaging correction (dG_cc):", dG_cc
         write(*,*) "restoring free energy (dG_res):", dG_res
         write(*,*) "Resulting chemical potential in mixture:", dG_is+dG_cc+dG_res
         write(*,*) "SMD Contribution (dG_CDS):", dG_disp
         write(*,*) "Systematic empirical shift (dG_shift)", dG_shift
         write(*,*) "-------------------------------------------------"
         write(*,*) "solvation free energy: ", dG_is+dG_cc+dG_res+dG_disp+dG_shift
      end if


     ! deallocate(solute_su,solute_sv,solute_svt,solute_sv0,solvent_su,solvent_sv,&
     !    &solvent_svt,solvent_sv0,solvent_area,solute_area,solvent_xyz,solute_xyz,&
     !    &solv_pot,sol_pot)
contains

subroutine get_arguments(config, error)
   type(configuration), intent(out) :: config
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   real(wp) :: val
   character(len=:), allocatable :: arg

   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      ! case("--help")
      !    call help(output_unit)
      !    stop
      ! case("--version")
      !    call version(output_unit)
      !   stop
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end select
   end do

    if (.not.(allocated(config%input))) then
       if (.not.allocated(error)) then
   !       call help(output_unit)
          error stop
       end if
    end if
   Call read_input(config,error)

end subroutine get_arguments 

   !> Subroutine to Read the COSMO-SACMD Input File
subroutine read_input(config,error)
   type(configuration) :: config
   type(error_type), allocatable :: error

   character(len=100) :: sac_param_path, smd_param_path, line

   integer :: io_error, i, n,j
   logical :: ex, started

   !> Check if the COSMO-SACMD Input File Exists.
   ex=.false.
   INQUIRE(file=config%input,exist=ex)
   IF (.NOT. ex) error stop "No Input File."

   !> Set Defaults
   config%T=298.15_wp
   config%ML=.FALSE.
   config%sig_in=.FALSE.
   config%prof=.FALSE.
   config%smd_default=.FALSE.

   Open(input_unit,file=config%input)
   Read(input_unit,'(A)',iostat=io_error,err=255) line
   Call move_line(line,config%sac_param_path)
   Read(input_unit,'(A)',iostat=io_error,err=255) line
   Call move_line(line,config%smd_param_path)
   Read(input_unit,'(A)',iostat=io_error,err=255) line

   j=1
   do i=1,len(trim(line))+1
      if (line(i:i) .EQ. " ") then
         select case(line(j:i-1))
            case('ML')
               config%ML=.true.
            case('sac','sac2010','sac2013')
               config%model=line(j:i-1)
            case('onlyprof')
               config%prof=.true.
            case('sigma_in','sig_in')
               config%sig_in=.true.
            case('smd_default', 'default_smd')
               config%smd_default=.true.
         end select
         j=i+1
      end if
   end do 

   Read(input_unit,'(A)',iostat=io_error,err=255) line !Comment Line
   Read(input_unit,'(A)',iostat=io_error,err=255) line
   Call move_line(line,config%csm_solvent)
   Read(input_unit,'(A)',iostat=io_error,err=255) line
   Call move_line(line,config%csm_solute)
   Read(input_unit,*,iostat=io_error,err=255) line, config%probe
   Call move_line(line,config%smd_solvent)
   Read(input_unit,*,iostat=io_error,err=255) config%T
   SysTemp=config%T 
   Read(input_unit,*,iostat=io_error,err=255) config%z1, config%z2
   
255 if (io_error .NE. 0) error stop "Check Input File."
end subroutine read_input

subroutine move_line(line,aline)
   character(*), intent(in) :: line
   character(:), allocatable, intent(inout) :: aline

   allocate(character(len(trim(line))) :: aline)
   aline=trim(line)
end subroutine move_line

end program COSMO

