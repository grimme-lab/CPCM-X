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

program COSMOX
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
   use qc_calc, only: qc_cal
   use mctc_env, only : wp, get_argument, fatal_error, error_type
   use crs_timer, only: timer_type, format_time
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use sdm
   implicit none
   character(len=*), parameter :: prog_name = "csx"
   integer :: oh_sol, nh_sol, near_sol
   real(wp), dimension(:), allocatable :: solute_su, solute_area, solute_sv, solute_sv0,solvent_pot,solute_pot
   real(wp), dimension(:), allocatable :: solvent_su, solvent_area, solvent_sv, solvent_sv0, solute_svt, solvent_svt
   real(wp), dimension(:), allocatable :: sol_pot, solv_pot
   real(wp), dimension(:,:), allocatable :: solvent_xyz, solute_xyz, solvat_xyz, solat_xyz, solat2
   character(2), dimension(:), allocatable :: solute_elements, solvent_elements, solute_hb, solvent_hb
   logical, dimension(:,:), allocatable :: solute_bonds, solvent_bonds
   logical, dimension(:), allocatable :: solute_rings
   real(wp), dimension(3,0:50) :: solvent_sigma3, solute_sigma3
   character(20) :: solvent, solute
   !real(wp), dimension(10) :: param
   real(wp) :: id_scr,gas_chem,chem_pot_sol, T, solute_volume, solvent_volume,&
      &solute_energy, solvent_energy, solvent_sigma(0:50), solute_sigma(0:50),sac_disp(2)
   logical :: gas,sig_in
   integer :: sol_nat, i
   integer, allocatable :: int_ident(:),solute_ident(:),solvent_ident(:)
   real(wp), allocatable :: surface(:), dsdr(:,:,:)

   type(timer_type) :: timer

   type :: configuration
      character(len=:), allocatable :: input
      character(len=:), allocatable :: smd_solvent
      character(len=:), allocatable :: csm_solvent
      character(len=:), allocatable :: csm_solute
      real(wp) :: T
      real(wp) :: probe
      real(wp) :: z1,z2
      real(wp) :: qc_eps
      character(len=:), allocatable :: sac_param_path
      character(len=:), allocatable :: smd_param_path
      character(len=:), allocatable :: database
      character(len=:), allocatable :: qc_calc
      character(len=:), allocatable :: xyz_input
      logical :: ML, sig_in, prof, smd_default, time
      character(len=:), allocatable :: model
   end type configuration

   type(configuration) :: config
   type(error_type), allocatable :: error

  

   type(DICT_STRUCT), pointer :: r_cav, disp_con
  
   gas=.TRUE.
   !! ------------------------------------------------------------ 
   !! Read Command Line Arguments and set Parameters accordingly
   !! ------------------------------------------------------------
   Call timer%push("total")
   Call get_arguments(config,error)
   if (allocated(error)) then
      write(error_unit,'(a)') error%message
      error stop
   end if
   Call echo_init(config)
   Call initialize_param(config%sac_param_path,config%model,r_cav,disp_con,config%csm_solvent) 
   if (config%ML) then
      Call init_pr
      write(*,*) "Machine Learning Mode selected. Will Only Write an ML.data file." !! ML Mode deprecated
      ML=.TRUE.
   end if
   !! ----------------------------------------------------------------------------------
   !! Read Sigma Profiles (--sigma) - Not the default case
   !! ----------------------------------------------------------------------------------
   T=config%T
   SysTemp=T
   if (config%sig_in) then

      write(*,*) "Reading Sigma Profile"
      Call read_singlesig(solvent_sigma,config%csm_solvent,solvent_volume)
      Call read_singlesig(solute_sigma,config%csm_solute,solute_volume)
      
      Call read_triplesig(solvent_sigma3,config%csm_solvent,solvent_volume)
      Call read_triplesig(solute_sigma3,config%csm_solute,solute_volume)
   else   
   !! ----------------------------------------------------------------------------------
   !! Creating COSMO Files with QC packages
   !! ----------------------------------------------------------------------------------
      if (allocated(config%qc_calc)) then
         Call timer%push("qc_calc")
         select case(config%qc_calc)
            case('TM')
               Call qc_cal(config%qc_eps,config%csm_solute,config%smd_solvent)
            case('ORCA')
               Call qc_cal(config%xyz_input,error)
            case default
               error stop "Chosen program "//config%qc_calc//"not supported"
         end select
         Call check_error(error)
         Call timer%pop() 
      end if 
   !! ----------------------------------------------------------------------------------
   !! Create the Sigma Profile from COSMO files
   !! ----------------------------------------------------------------------------------
      
      Call timer%push("sigma_av")
   !! ------------------------------------------------------------------------------------
   !! Read necessary COSMO Data
   !! ------------------------------------------------------------------------------------
      write(output_unit,'(10x,a)') &
         " ------------------------------------------------- ",&
         "|                   Calculation                   |",&
         " ------------------------------------------------- ", &
         ""
      write(output_unit,'(5x,a)') "Reading COSMO data."
      Call read_cosmo(config%csm_solvent,solvent_elements,solvent_ident,solvent_xyz,solvent_su,&
          &solvent_area,solvent_pot,solvent_volume,solvent_energy,solvat_xyz,config%database)
      Call read_cosmo(config%csm_solute,solute_elements,solute_ident, solute_xyz, solute_su,&
         &solute_area,solute_pot,solute_volume,solute_energy,solat_xyz,config%database)
 
   !! ------------------------------------------------------------------------------------
   !! Sigma Charge Averaging and creating of a single Sigma Profile for Solute and Solvent
   !! ------------------------------------------------------------------------------------
      write(output_unit,'(5x,a)') "Creating Sigma Profile from COSMO data."
      Call average_charge(param(1), solvent_xyz,solvent_su,solvent_area,solvent_sv)
      Call average_charge(param(1), solute_xyz, solute_su, solute_area, solute_sv)
      Call single_sigma(solvent_sv,solvent_area,solvent_sigma,"solvent")
      Call single_sigma(solute_sv,solute_area,solute_sigma,"solute")
      Call timer%pop()
   !! ------------------------------------------------------------------------------------
   !! Determination of HB Grouping and marking of Atom that are able to form HBs.
   !! Determination of Atoms in Rings, necessary for the PR2018 EOS (only ML Model)
   !! ------------------------------------------------------------------------------------
   if ((config%ML) .OR. (.NOT. config%model .EQ. "sac")) then
      Call timer%push("bondings") 
      write(output_unit,'(5x,a)') "Determine Ring atoms and HB groups."
      Call det_bonds(solute_ident,solat_xyz,solute_elements,solute_bonds,oh_sol,nh_sol)
      Call hb_grouping(solute_ident,solute_elements,solute_bonds,solute_hb)
      Call det_bonds(solvent_ident,solvat_xyz,solvent_elements,solvent_bonds)
      Call hb_grouping(solvent_ident,solvent_elements,solvent_bonds,solvent_hb)
      
      Call det_rings(solute_ident,solute_bonds,solute_rings,near_sol)
      Call timer%pop()
   end if

   !! ------------------------------------------------------------------------------------
   !! Creation of a splitted Sigma Profile, necessary for sac2010/sac2013
   !! ------------------------------------------------------------------------------------

      if (.NOT. (config%model .EQ. "sac")) then
         Call timer%push("sigma_split")
         Call split_sigma(solvent_sv,solvent_area,solvent_hb,solvent_ident,solvent_elements,&
            &solvent_sigma3,"solvent")
         Call split_sigma(solute_sv,solute_area,solute_hb,solute_ident,solute_elements,&
            &solute_sigma3,"solute")
         Call timer%pop()
      end if

   !! ------------------------------------------------------------------------------------
   !! Exit here if you only want Sigma Profiles to be created 
   !! ------------------------------------------------------------------------------------
      if (config%prof) then;
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
         write(output_unit,'(5x,a)') "Calculating Gas Phase Energies.", &
            ""
         Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)         
         !if (gas) then
         !   Call calcgas(solute_energy,id_scr,gas_chem,solute_area,solute_sv,solute_su,&
         !      &solute_pot,solute_elements,solute_ident,disp_con, T,r_cav)
         !end if

         ! Computation of COSMO-RS equations (here may be something wrong atm)
         Call timer%push("solv")
         Call compute_solvent(solv_pot,solvent_sv,solvent_svt,solvent_area,T,500,0.0001,solvent_ident,solvent_hb)
         Call timer%pop()
         Call timer%push("solu")
         Call compute_solute(sol_pot,solv_pot,solute_sv,solute_svt,solvent_sv,&
         &solvent_svt,solute_area,solvent_area,T,chem_pot_sol,solute_ident,solvent_ident,solute_elements,solvent_hb)
         Call timer%pop()
         allocate (int_ident(maxval(solute_ident)))
         do i=1,maxval(solute_ident)
            int_ident(i)=i
         end do
         
         Call timer%push("cds")
         Call calculate_cds(int_ident,solute_elements,solat_xyz,config%probe,&
         &config%smd_solvent,config%smd_param_path,config%smd_default)
         Call timer%pop()

         dG_res=chem_pot_sol+param(9)*near_sol
         deallocate(solute_su,solute_sv,solute_sv0,solvent_su,solvent_sv,&
         &solvent_sv0,solvent_area,solute_area,solvent_xyz,solute_xyz,&
         &solv_pot,sol_pot)
      end select
      
      write(output_unit,'(10x,a)') &
         " ------------------------------------------------- ",&
         "|                     Results                     |",&
         " ------------------------------------------------- ", &
         ""
      
      if (config%ML) then
         write(*,*) "Writing ML data in ML.data"
         Call System("paste --delimiters='' ML.energy ML.gamma ML.pr > ML.data")
         Call System ("rm ML.energy ML.pr")
      else
         write(output_unit,'(4x,a)') repeat('-',73)
         write(output_unit,'(5x,a,t55,a,t66,a)') &
            "Free Energy contributions:", "[Eh]", " [kcal/mol]"
         write(output_unit,'(5x,a,t50,E13.5,t65,F10.5)') &
         "Ideal State (dG_is):", dG_is/autokcal, dG_is, &
         "Averaging correction (dG_cc):", dG_cc/autokcal, dG_cc, &
         "restoring free energy (dG_res):", dG_res/autokcal, dG_res, &
         "Resulting chemical potential in mixture:", (dG_is+dG_cc+dG_res)/autokcal,&
         & dG_is+dG_cc+dG_res, &
         "SMD Contribution (dG_CDS):", dG_disp/autokcal, dG_disp, &
         "Systematic empirical shift (dG_shift)", dG_shift/autokcal, dG_shift
         write(output_unit,'(4x,a)') repeat('-',73)
         write(output_unit,'(5x,a,t50,E13.5,t65,F10.5)') &
         "solvation free energy: ", (dG_is+dG_cc+dG_res+dG_disp+dG_shift)/autokcal,&
         & dG_is+dG_cc+dG_res+dG_disp+dG_shift
         write(output_unit,*) ""
      end if

      Call timer%pop()

      if (config%time) then
      block
         integer :: i
         real(wp) :: ttime, stime
         character(len=*), parameter :: label(*) = [character(len=20) :: &
         & "qc_calc","sigma_av","bondings","sigma_split","solv","solu","cds"]
         ttime=timer%get("total")
         write(*,*)
         write(*,*) "Timings:"
         write(*,*) "total"//repeat(" ", 18)//format_time(ttime)
         do i = 1,size(label)
            stime = timer%get(label(i))
            write(*,*) label(i)//repeat(" ",3)//format_time(stime) 
         end do
      end block
      end if


     ! deallocate(solute_su,solute_sv,solute_svt,solute_sv0,solvent_su,solvent_sv,&
     !    &solvent_svt,solvent_sv0,solvent_area,solute_area,solvent_xyz,solute_xyz,&
     !    &solv_pot,sol_pot)
contains

subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options]/<inputfile>"

   write(unit, '(a)') &
      "Calculates the solvation free energy of a compound in a solvent.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "    --solvent", "Specify a solvent and uses configuration given in config.toml.", &
      "    --newrc", "Creates a sample configuration file (config.toml).", &
      "    --newinput", "Creates a sample input file csx.input in the currect working directory.", &
      "    --keyword", "Shows a list of possible Keywords for the csx.input file.", &
      "    --inp", "Allows to specify an input file for advanced configuration.", &
      "    --help", "Show this help message"
   write(unit, '(a)')
   
end subroutine help

subroutine sample(filename,rc)
   !> Path, where the Sample should be created.
   character(len=*), intent(in) :: filename
   !> Is TRUE for a sample configuration file and FALSE for a sample input.
   logical, intent(in) :: rc

   integer :: unit

   open(newunit=unit,file=filename)
   if (rc) then
      write(unit,'(a)') &
         '# This is a Sample CSX Configuration File', &
         '# To work with this file, an environmental variable CSXHOME has to be set.', &
         '# This File has to be placed in the CSXHOME path.', &
         '# All Parameters and the Database need to be placed in CSXHOME.', &
         '', &
         '# Path or Filename for the parameters for Water. (Path in respective to CSXHOME)', &
         'smd_h2o="smd_h2o"', &
         'crs_h2o="crs.param_h2o"', &
         '', &
         '# Path or Filename for the parameters for other solvents.', &
         'smd_ot="smd_ot"', &
         'crs_ot="crs.param_ot"', &
         '', &
         '# Path to the Solvent Database', &
         '# Note, that the .cosmo files in the database need to be named according to SMD solvent names.', &
         'DB="DATABASE-COSMO"', &
         '', &
         '# Temperature in Kelvin', &
         'Temperature=298.15', &
         '', &
         '#Probe radius for the CDS term in Å (default=0.4)', &
         'r_probe=0.4'
   else
      write(unit,'(a)') &
         "/path/to/crs/parameter/file.param   &
         &#This needs to directly point to the respective parameter file (crs.param_h2o or crs.param_ot)", &
         "/path/to/smd/parameters   #This needs to point to the folder, where the smd_h2o and smd_ot files are included.", &
         "KEYWORDS", &
         "#Comment line.", &
         "/path/to/solvent.cosmo   #Needs to point at the correct Solvent from the Database", &
         "solute.cosmo   #Needs to point to the solute. Will be automatically created with the TM keyword.", &
         "solvent 0.4   #Will set the solvent used for the SMD(CDS) part and the probe radius (default=0.4).", &
         "298.15   #Sets the temperature."
   end if
   close(unit)

end subroutine sample

subroutine print_keywords(unit)
   integer, intent(in) :: unit

   write(unit, '(a)') &
      "Keywords are used in the keyword line in the .input file.", &
      "", &
      "Available Keywords:", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "    crs", "Invokes the standard COSMO-X model.", &
      "    sac, sac2010", "Invokes an SAC based model with or without HB splitting (needs different parameters).", &
      "    TM/TM=epsilon", "Starts with single point calculation for the solute. Needs control file. (default: epsilon=infinity)", &
      "    time", "Shows additional Information about the time needed for various steps of the algorithm.", &
      "    onlyprof", "Only calculates a Sigma Profile and prints it in a .sigma file.", &
      "    sigma_in", "Expects Sigma Profiles instead of .cosmo files (only for SAC based models).", &
      "    smd_default", "Uses SMD default Parameters instead of fitted Parameters (use only if you know what you are doing).", &
      "    DB=path", "Optionally defines the Path to a COSMO file database (e.g. DATABASE-COSMO)"
   write(unit, '(a)')
end subroutine print_keywords


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
      case("--help", "-h")
         call help(output_unit)
         stop
      case("--newinput")
         call sample("csx.input",.FALSE.)
         write(output_unit, '(a)') "[Info] Sample input file 'csx.input' created."
         stop
      case("--newrc")
         call sample("config.toml",.TRUE.)
         write(output_unit, '(a)') "[Info] Sample config file 'config.toml' created."
         stop
      case ("--keyword", "--keywords")
         call print_keywords(output_unit)
         stop
      case ("--solvent", "-s", "--solv")
         iarg=iarg+1
         call get_argument(iarg,arg)
         if ((.not.allocated(config%smd_solvent)) .AND. (.not.allocated(config%input))) then
            call use_default(config,arg,error)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      ! case("--version")
      !    call version(output_unit)
      !   stop
      case ("--inp", "--input")
         iarg=iarg+1
         call get_argument(iarg,arg)
         if ((.not.allocated(config%input)) .AND. (.not.allocated(config%smd_solvent))) then
            call move_alloc(arg, config%input)
            Call read_input(config,error)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case default
         if ((.not.allocated(config%xyz_input))) then
            call move_alloc(arg, config%xyz_input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end select
   end do

    if ((.not.(allocated(config%input))) .AND. (.not. (allocated(config%smd_solvent)))) then
       if (.not.allocated(error)) then
          call help(output_unit)
          stop
       end if
    end if

end subroutine get_arguments 

   !> Subroutine to Read the COSMO-SACMD Input File
subroutine read_input(config,error)
   type(configuration) :: config
   type(error_type), allocatable :: error

   character(len=100) :: sac_param_path, smd_param_path, line

   integer :: io_error, i, n,j, equal
   logical :: ex, started

   character(len=:), allocatable :: keyword, substring

   !> Set Defaults
   config%T=298.15_wp
   config%ML=.FALSE.
   config%sig_in=.FALSE.
   config%prof=.FALSE.
   config%smd_default=.FALSE.
   config%time=.FALSE.
   config%qc_eps=0
   config%probe=0.4
   !> Check if the COSMO-SACMD Input File Exists.
   ex=.false.
   INQUIRE(file=config%input,exist=ex)
   IF (.NOT. ex) then
      write(output_unit,'(a)') "Error: No Input File defined.", &
      ""
      Call help(output_unit)
      stop
   end if

   Call move_line("NONE",config%database)

   Open(input_unit,file=config%input)
   Read(input_unit,'(A)',iostat=io_error,err=255) line
   Call move_line(line,config%sac_param_path)
   Read(input_unit,'(A)',iostat=io_error,err=255) line
   Call move_line(line,config%smd_param_path)
   Read(input_unit,'(A)',iostat=io_error,err=255) line

   j=1
   do i=1,len(trim(line))+1
      if (line(i:i) .EQ. " ") then
         equal=index(line(j:i-1),'=')
         if (equal .ne. 0) then
            Call move_line(line(j:(j+equal-2)),keyword)
            Call move_line(line((j+equal):i-1),substring)
         else
            Call move_line(line(j:i-1),keyword)
         end if
         
         select case(keyword)
            case ('TM','tm')
               if (allocated(config%qc_calc)) then
                  Call fatal_error(error,"Too many Arguments for QC calculation.")
                  return
               end if
               Call move_line("TM",config%qc_calc)
               if (equal .ne. 0) then
                  select case(substring)
                     case ('default','minnesota')
                        config%qc_eps=-1
                     case ('infinity')
                        config%qc_eps=0
                     case default
                        read(substring,*) config%qc_eps
                  end select
               end if
            case ('ORCA','Orca','orca')
               if (allocated(config%qc_calc)) then
                  Call fatal_error(error,"Too many Arguments for QC calculation.")
                  return
               end if
               Call move_line("ORCA",config%qc_calc)
               if (equal .ne. 0) then
                  select case(substring)
                     case ('default','minnesota')
                        config%qc_eps=-1
                     case ('infinity')
                        config%qc_eps=0
                     case default
                        read(substring,*) config%qc_eps
                  end select
               end if
            case ('DB')
               if (equal .ne. 0) Call move_line(substring,config%database)
            case('time')
               config%time=.true.
            case('ML')
               config%ML=.true.
            case('sac','sac2010','sac2013','crs')
               config%model=line(j:i-1)
            case('onlyprof')
               config%prof=.true.
            case('sigma_in','sig_in')
               config%sig_in=.true.
            case('smd_default', 'default_smd')
               config%smd_default=.true.
         end select
         deallocate(keyword)
         if (allocated(substring)) deallocate(substring)
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
   config%z1=0.995
   config%z2=0.005
   
255 if (io_error .NE. 0) error stop "Check Input File."
end subroutine read_input

subroutine use_default(config, solv, error)
   use mctc_env_system, only : get_variable
   use tomlf, only : toml_table, toml_parse, toml_error, toml_key, get_value
   !> Solvent used for default configuration
   character(:), allocatable, intent(inout) :: solv
   !> Configuration Type
   type(configuration), intent(inout) :: config
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Toml unit
   type(toml_table), allocatable :: config_table
   !> Toml error handling
   type(toml_error), allocatable :: config_error
   
   
   type(toml_key), allocatable, dimension(:) :: list
   integer :: stat
   character(len=255) :: line
   character(:), allocatable :: home, line2
   logical :: ex
   character(len=10) :: control, command
   integer :: nconf
   
   !> Set Defaults
   config%T=298.15_wp
   config%ML=.FALSE.
   config%sig_in=.FALSE.
   config%prof=.FALSE.
   config%smd_default=.FALSE.
   config%time=.FALSE.
   config%qc_eps=0
   config%probe=0.4
   if (solv .eq. "h2o") call move_line("water",solv)
   call move_line(solv,config%smd_solvent)
   call move_line(solv//".cosmo",config%csm_solvent)
   call move_line("solute.cosmo",config%csm_solute)
   call move_line("crs",config%model)
   Call get_variable("CSXHOME",home)
   
   if (.not.allocated(home)) then
      call fatal_error(error, "CSXHOME Variable ist not set.")
      RETURN
   end if

   ex=.false.
  
   if (home(len(home):len(home)) .ne. "/") call move_line(home//"/",home)
  
   INQUIRE(file=home//"config.toml",exist=ex)

   if (.not. ex) then 
      call fatal_error(error, "No config.toml found in "//home)
      return
   end if

   open(input_unit,file=home//"config.toml")
   call toml_parse(config_table,input_unit,config_error)
   close(input_unit)
   if (allocated(config_table)) then
      call config_table%get_keys(list)
      do nconf=1,size(list)
         select case(list(nconf)%key)
         case("prog") 
               call get_value(config_table,list(nconf),line2)
               call move_line(line2,config%qc_calc)
         case("smd_h2o")
            if (solv .eq. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(home//line2,config%smd_param_path)
            end if
         case("smd_ot")
            if (solv .ne. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(home//line2,config%smd_param_path)
            end if
         case("crs_h2o")
            if (solv .eq. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(home//line2,config%sac_param_path)
            end if
         case("crs_ot")
            if (solv .ne. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(home//line2,config%sac_param_path)
            end if
         case("DB")
            call get_value(config_table,list(nconf),line2)
            call move_line(home//line2,config%database)
         case("Temperature")
            call get_value(config_table,list(nconf),config%T)
         case("r_probe")
            call get_value(config_table,list(nconf),config%probe)
         case default
            call fatal_error(error, "Unrecognized key in your config file: "//list(nconf)%key)
            return
         end select
      end do
   else
      call fatal_error(error, config_error%message)
   end if

   ! do while (.TRUE.)
   !    if (index(line,"$") .eq. 0) read(11,'(a)',iostat=stat) line
   !    write(*,*) line
   !    if (stat .lt. 0) exit
   !    if (index(line,"$") .eq. 0) cycle
   !    line=line(2:len(line))
   !    select case(solv)
   !    case("H2O", "Water", "water")
   !       if (trim(line) .eq. "H2O") then
   !          do while (.TRUE.)
   !             read(11,'(a)',iostat=stat) line
   !             if (index(line, "$") .ne. 0) exit
   !             control=line(1:index(line,"=")-1)
   !             command=trim(line(index(line,"=")+1:len(line)))
   !             select case(trim(control))
   !             case("smd")
   !                call move_line(trim(command),config%smd_param_path)
   !             case("crs")
   !                call move_line(trim(command),config%sac_param_path)
   !             end select
   !          end do
   !       end if
   !    case default
   !       if (trim(line) .eq. "OT") then
   !          do while (.TRUE.)
   !             read(11,'(a)',iostat=stat) line
   !             if (index(line, "$") .ne. 0) exit
   !             control=line(1:index(line,"=")-1)
   !             command=trim(line(index(line,"=")+1:len(line)))
   !             select case(trim(control))
   !             case("smd")
   !                call move_line(trim(command),config%smd_param_path)
   !             case("crs")
   !                call move_line(trim(command),config%sac_param_path)
   !             end select
   !          end do
   !       end if
   !    end select

   !    if (trim(line) .eq. "DB") then
   !       read(11,'(a)',iostat=stat) line
   !       control=line(1:index(line,"=")-1)
   !       command=trim(line(index(line,"=")+1:len(line)))
   !       call move_line(trim(command),config%database)
   !    end if

   ! end do


end subroutine use_default

subroutine move_line(line,aline,hignore)
   !> Line to write into the allocatable unit
   character(*), intent(in) :: line
   !> Ignores everything after an hashtag (default=true)
   logical, intent(in), optional :: hignore
   !> Allocatable character array to be set to line
   character(:), allocatable, intent(inout) :: aline

   integer :: i
   logical :: ignore

   ignore=.true.

   if (present(hignore)) then
      if (.not. hignore) ignore=.false.
   end if 

   if (allocated(aline)) deallocate(aline)

   if (ignore) then
      do i= 1,len(trim(line))
         if (line(i:i) .EQ. "#") then 
            allocate(character(len(trim(line(1:i-1)))) :: aline)
            aline=trim(line(1:i-1))
            exit
         end if 
         if (i .EQ. len(trim(line))) then
            allocate(character(len(trim(line(1:i)))) :: aline)
            aline=trim(line(1:i))
            exit
         end if 
      end do
   else
      allocate(character(len(trim(line))) :: aline)
      aline=trim(line)
   end if
end subroutine move_line

subroutine echo_init(config)
   type(configuration) :: config


   write(output_unit,'(10x,a)') &
      " ------------------------------------------------- ",&
      "|                 Initialization                  |",&
      " ------------------------------------------------- ", &
      ""

   write(output_unit,'(5x,a,t35,a)') &
      "SMD Parameter Path:", config%smd_param_path, &
      "CRS Parameter Path:", config%sac_param_path, &
      "Solvent:", config%smd_solvent, &
      "Corresponding COSMO File:", config%csm_solvent

   if (.NOT. allocated(config%qc_calc)) write(output_unit,'(5x,a,t35,a)') &
                           "Solute COSMO File:", config%csm_solute
      
end subroutine echo_init

subroutine check_error(error)
   type(error_type), intent(in), allocatable :: error
   
   if (allocated(error)) then
      write(error_unit,'(a)') error%message
      error stop
   end if

end subroutine check_error

end program COSMOX

