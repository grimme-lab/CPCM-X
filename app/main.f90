! This file is part of CPCM-X.
! SPDX-Identifier: LGPL-3.0-or-later
!
! CPCM-X is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CPCM-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with CPCM-X.  If not, see <https://www.gnu.org/licenses/>.

program CPCMX
   use cpx, only: calculation_type, parameter_type, initialize_param, load_solvent, read_cosmo, init_pr,&
      &AtomicMass, Density
   use globals
   use isodens, only: get_isodens_radii
   use qc_calc, only: qc_cal, orcatocosmo
   use mctc_env, only : wp, get_argument, fatal_error, error_type
   use crs_timer, only: timer_type, format_time
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   implicit none
   character(len=*), parameter :: prog_name = "cpx"
   logical :: ex
   integer :: ioerror

   real(wp), allocatable :: isodens_rad(:)
   integer :: i


   type(timer_type) :: timer

   type(calculation_type) :: calc
   type(configuration_type) :: config
   type(parameter_type) :: parameter
   type(error_type), allocatable :: error

  

  
   !! ------------------------------------------------------------ 
   !! Read Command Line Arguments and set Parameters accordingly
   !! ------------------------------------------------------------
   Call timer%push("total")
   Call get_arguments(config,error)
   Call check_error(error)
   Call echo_init(config)
   if (config%internal) then
      Call initialize_param(config%qc_calc,config%smd_solvent,calc,error)
   else
      Call initialize_param(config%sac_param_path,calc%param,error)
   end if
   Call check_error(error) 
   if (config%ML) then
      Call init_pr
      write(*,*) "Machine Learning Mode selected. Will Only Write an ML.data file." !! ML Mode deprecated
      ML=.TRUE.
   end if
   !! ----------------------------------------------------------------------------------
   !! Read Sigma Profiles (--sigma) - Not the default case
   !! ----------------------------------------------------------------------------------
   SysTemp=config%T
   !if (config%sig_in) then
   !
   !   write(*,*) "Reading Sigma Profile"
   !   Call read_singlesig(solvent_sigma,config%csm_solvent,solvent_volume)
   !   Call read_singlesig(solute_sigma,config%csm_solute,solute_volume)
   !   
   !   Call read_triplesig(solvent_sigma3,config%csm_solvent,solvent_volume)
   !   Call read_triplesig(solute_sigma3,config%csm_solute,solute_volume)
   !else   
   !! ----------------------------------------------------------------------------------
   !! Creating COSMO Files with QC packages
   !! ----------------------------------------------------------------------------------
      if (allocated(config%qc_calc)) then
         Call timer%push("qc_calc")
         select case(config%qc_calc)
            case('tm')
               Call qc_cal(config%qc_eps,config%csm_solute, error, config%isodens, config%smd_solvent)
            case('orca')
               Call qc_cal(config%xyz_input,error)
            case('P-gTB', 'gtb')
               Call qc_cal(0.6_wp,0.20_wp,error)
            case('xtb')
               Call qc_cal('inf','2',error)
            case default
               write(error_unit,'(a,a,a)') "Chosen program "//config%qc_calc//" not supported"
               error stop
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
      if (config%internal) then
         Call load_solvent(config%smd_solvent,calc%solvent,error)
      else
         Call read_cosmo(config%csm_solvent,calc%solvent,config%database,error)
      end if
      call check_error(error)
      Call read_cosmo(config%csm_solute,calc%solute,config%database,error)
      call check_error(error)
 
   !! ------------------------------------------------------------------------------------
   !! Sigma Charge Averaging and creating of a single Sigma Profile for Solute and Solvent
   !! ------------------------------------------------------------------------------------
      write(output_unit,'(5x,a)') "Creating Sigma Profile from COSMO data."
      Call calc%average_charge(error)
      call check_error(error)
      Call calc%sigma(.true.)
      Call timer%pop()
   !! ------------------------------------------------------------------------------------
   !! Determination of HB Grouping and marking of Atom that are able to form HBs.
   !! Determination of Atoms in Rings, necessary for the PR2018 EOS and ER Correction
   !! ------------------------------------------------------------------------------------
   if ((config%ML) .OR. (.NOT. config%model .EQ. "sac")) then
      Call timer%push("bondings") 
      write(output_unit,'(5x,a)') "Determine Ring atoms and HB groups."
      Call calc%init_bonding()
      Call timer%pop()
   end if

   !! ------------------------------------------------------------------------------------
   !! Creation of a splitted Sigma Profile, necessary for sac2010/sac2013
   !! ------------------------------------------------------------------------------------

   if (.NOT. (config%model .EQ. "sac")) then
      Call timer%push("sigma_split")
      Call calc%split_sigma(.true.)
      Call timer%pop()
   end if

   !! ------------------------------------------------------------------------------------
   !! Exit here if you only want Sigma Profiles to be created 
   !! ------------------------------------------------------------------------------------
   if (config%prof) then;
      write(*,*) "Only Profile mode choosen, exiting."
      stop
   end if
   !end if
   

   !! ------------------------------------------------------------------------------------
   !! Choice of the different post COSMO Models (sac,sac2010,sac2013,CPCM-X)
   !! CPCM-X is currently the only recommended Model (and since the refactor the only working one).
   !! ------------------------------------------------------------------------------------

   !select case (trim(config%model))
   !   case ("sac")
   !      !Calculation of the Gas Phase (ideal gas --> ideal conductor)
   !      Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)
   !      !Calculation of the Solvent Phase (ideal conductor --> real solution)
   !      Call sac_2005(solvent_sigma,solute_sigma,solvent_volume,solute_volume,config%z1,config%z2)
   !      !Calculation of NES contributions (real gas --> ideal gas?)
   !      if (config%ML) Call pr2018(solute_area,solute_elements,solute_ident,oh_sol,nh_sol,near_sol)
   !
   !      allocate (int_ident(maxval(solute_ident)))
   !      do i=1,maxval(solute_ident)
   !         int_ident(i)=i
   !      end do
   !
   !      Call calculate_cds(int_ident,solute_elements,solat_xyz,config%probe,&
   !      &config%smd_solvent,config%smd_param_path,config%smd_default)
   !
   !   case("sac2010")
   !      
   !      Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)
   !      Call sac_2010(solvent_sigma3,solute_sigma3,solvent_volume,solute_volume)
   !      if (config%ML) Call pr2018(solute_area,solute_elements,solute_ident,oh_sol,nh_sol,near_sol)
   !
   !      allocate (int_ident(maxval(solute_ident)))
   !      do i=1,maxval(solute_ident)
   !         int_ident(i)=i
   !      end do
   !
   !      Call calculate_cds(int_ident,solute_elements,solat_xyz,config%probe,&
   !      &config%smd_solvent,config%smd_param_path,config%smd_default)
   !! ------------------------------------------------------------------------------------
   !! The SAC 2013 Routine is not fully implemented and not supported anymore
   !! ------------------------------------------------------------------------------------
   !  case("sac2013")
   !
   !    Call sac_gas(solute_energy,id_scr,solute_area,solute_sv,solute_su,solute_pot)
   !    Call sac2013_disp(trim(solvent),solvent_bonds,solvent_ident,solvent_elements,disp_con,sac_disp(1))
   !    Call sac2013_disp(trim(solute),solute_bonds,solute_ident,solute_elements,disp_con,sac_disp(2))
   !    Call sac_2013(solvent_sigma3,solute_sigma3,solvent_volume,solute_volume,sac_disp)
   !    Call pr2018(solute_area,solute_elements,solute_ident,oh_sol,nh_sol,near_sol)
   !  case ("crs")
   
         !! CPCM-RS calculation starts here !!

         ! Calculation of Solvent Phase energies
         write(output_unit,'(5x,a)') "Calculating Solvent Phase Energies.", &
            ""

         !> Read in gas phase energy
         open(1,file="gas.energy",iostat=ioerror)
         read(1,*,iostat=ioerror) calc%solute%energy_gas
         if (ioerror .NE. 0) Call fatal_error(error,"Problem while reading energies (check gas.energy file).")
         Call check_error(error)

         Call timer%push("calc")
         Call calc%solv("crs",error,T=298.15_wp,max_cycle=500,conv_crit=0.0001_wp)
         Call check_error(error)
         call calc%state_correction(density(config%smd_solvent),AtomicMass(calc%solvent%element),config%T)
         Call timer%pop()

         write(output_unit,'(10x,a)') &
         " ------------------------------------------------- ",&
         "|                Gas Phase Results                |",&
         " ------------------------------------------------- ", &
         ""
         
         write(output_unit,'(5x,a,t30,F15.8,2x,a)') &
         "E_COSMO:", calc%solute%energy, "Eh", &
         "E_COSMO+dE:", (calc%solute%energy+calc%dG_cc),"Eh",   &
         "E_gas:", calc%solute%energy,"Eh", &
         "E_COSMO-E_gas", (calc%id_scr-calc%solute%energy),"Eh", &
         "E_COSMO-E_gas+dE:", (calc%id_scr-calc%solute%energy+calc%dG_cc),"Eh", &
         "Averaging corr dE:", calc%dG_cc,"Eh", &
         "thermostatic correction: ", calc%param%eta*R*jtokcal*config%T, "Eh", &
         "Area:", sum(calc%solute%area), "Å"
         write(output_unit,'(5x,a,t30,F15.8,2x,a)') &
         "Solvent density:", density(config%smd_solvent), "g/l", &
         "Solvent atomic mass:", AtomicMass(calc%solvent%element), "a.u.", &
         "State correction", calc%dG_ss, "Eh", &
         ""
         
         Call timer%push("cds")
         if (config%isodens) then
            Call get_isodens_radii(calc%solute%xyz,calc%solute%id,calc%solute%atom_xyz,isodens_rad)
            write(output_unit,'(10x,a)') &
            " ------------------------------------------------- ",&
            "|                 Isodensity Radii                 |",&
            " ------------------------------------------------- ", &
            ""
            write(output_unit,'(5x,a)'), &
            "Isodensity Flag used, calculated isodensity radii:",&
            ""
            write(output_unit,'(10x,a,t30,a)'), &
               "Atom Number:", "[A]"
            do i=1,maxval(calc%solute%id)
               write(output_unit,'(10x,I0,t30,F4.2)'),&
                  i, isodens_rad(i)
            end do
            write(output_unit,'(a)') ""
            Call calc%cds(config%probe,&
            &config%smd_solvent,config%smd_param_path,isodens_rad)
         else
            if (allocated(calc%smd_param)) then
               Call calc%cds(config%probe,config%smd_solvent)
            else
               Call calc%cds(config%probe,&
               &config%smd_solvent,config%smd_param_path,config%smd_default)
            end if
         end if
         Call timer%pop()

      !end select
      
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
         "Ideal State (dG_is):", calc%dG_is, calc%dG_is*autokcal, &
         "Averaging correction (dG_cc):", calc%dG_cc, calc%dG_cc*autokcal, &
         "restoring free energy (dG_res):", calc%dG_res, calc%dG_res*autokcal, &
         "Resulting chemical potential in mixture:", &
         &calc%dG_is+calc%dG_cc+calc%dG_res, &
         &(calc%dG_is+calc%dG_cc+calc%dG_res)*autokcal, &
         "SMD Contribution (dG_CDS):", calc%dG_smd, calc%dG_smd*autokcal, &
         "Standard state correction:", calc%dG_ss, calc%dG_ss*autokcal, &
         "Systematic empirical shift (dG_shift)", calc%dG_shift, calc%dG_shift*autokcal
         write(output_unit,'(4x,a)') repeat('-',73)
         write(output_unit,'(5x,a,t50,E13.5,t65,F10.5)') &
         "solvation free energy: ", calc%dG(), calc%dG()*autokcal
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
contains

subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options]/<inputfile>"

   write(unit, '(a)') &
      "Calculates the solvation free energy of a compound in a solvent.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "    --solvent", "Specify a solvent and uses configuration given in cpcmx.toml.", &
      "","For Orca this needs an .xyz file as input (e.g. csx inp.xyz --solvent water)",&
      "    --prog", "Overwrites the qc program chosen in the configuration file.", &
      "    --newrc", "Creates a sample configuration file (cpcmx.toml).", &
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
         '# This is a Sample CPX Configuration File', &
         '# To work with this file, an environmental variable CPXHOME has to be set.', &
         '# This File has to be placed in the CPXHOME or in your home path.', &
         '# Database need to be placed in CPXHOME.', &
         '', &
         '# Default QC Program for the Single Point Calculations', &
         'prog="TM"', &
         '', &
         '# Path or Filename for the parameters for Water. (Path in respective to DB)', &
         '# Uses CPXHOME, if Database is not specified', &
         'smd_h2o="smd_h2o"', &
         'crs_h2o="crs.param_h2o"', &
         '', &
         '# Path or Filename for the parameters for other solvents.', &
         'smd_ot="smd_ot"', &
         'crs_ot="crs.param_ot"', &
         '', &
         '# Path to the Solvent Database in respective to CPXHOME', &
         '# Note, that the .cosmo files in the database need to be named according to SMD solvent names.', &
         '# DB path needs to contain subfolders for the QC programs, if used with qc mode.', &
         'DB="DB"', &
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
      "    crs", "Invokes the standard CPCM-X model.", &
      "    sac, sac2010", "Invokes an SAC based model with or without HB splitting (needs different parameters).", &
      "    TM/TM=epsilon", "Starts with single point calculation for the solute. Needs control file. (default: epsilon=infinity)", &
      "    ORCA", "Starts with single point calculation for the solute with epsilon=infinity. Needs .xyz file.", &
      "    time", "Shows additional Information about the time needed for various steps of the algorithm.", &
      "    onlyprof", "Only calculates a Sigma Profile and prints it in a .sigma file.", &
      "    sigma_in", "Expects Sigma Profiles instead of .cosmo files (only for SAC based models).", &
      "    smd_default", "Uses SMD default Parameters instead of fitted Parameters (use only if you know what you are doing).", &
      "    DB=path", "Optionally defines the Path to a COSMO file database (e.g. DATABASE-COSMO)"
   write(unit, '(a)')
end subroutine print_keywords


subroutine get_arguments(config, error)
   use mctc_env_system, only : get_variable
   use data, only: solvent_name
   type(configuration_type), intent(out) :: config
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   real(wp) :: val
   character(len=:), allocatable :: arg, home
   logical :: ex


   Call get_variable("CPXHOME",home)
   
   if (.not.allocated(home)) then
      Call move_line("/",home)
   else
      if (home(len(home):len(home)) .ne. "/") call move_line(home//"/",home)
   end if
   ex=.false.
  

   config%isodens=.false.
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
      case("--isodens")
         config%isodens=.true.
      case("--newrc")
         call sample("cpcmx.toml",.TRUE.)
         write(output_unit, '(a)') "[Info] Sample config file 'cpcmx.toml' created."
         stop
      case ("--keyword", "--keywords")
         call print_keywords(output_unit)
         stop
      case ("--solvent", "-s", "--solv")
         iarg=iarg+1
         call get_argument(iarg,arg)
         if ((.not.allocated(config%smd_solvent)) .AND. (.not.allocated(config%input))) then
            call use_default(config,solvent_name(arg),home,error)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case ("--prog") 
         iarg=iarg+1
         call get_argument(iarg,arg)
         call move_alloc(arg, config%qc_calc)
      ! case("--version")
      !    call version(output_unit)
      !   stop
      case ("--internal")
         config%internal=.true.
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

   if (config%internal) then
      select case (config%qc_calc)
      case default
         Call fatal_error(error,"Internal mode only available for xTB")
         return
      case ('xtb','xTB')
         return
      end select
   end if

   if (config%isodens) then
      config%qc_calc="tm"
      inquire(file=config%database//"/isodens/"//config%sac_param_path, exist=ex) 
      if (ex) then
         config%database=config%database//"/isodens"
      else
         call fatal_error(error,"Database for isodens mode does not seem to be set up properly.")
         return
      end if
   end if

   if ((allocated(config%qc_calc)) .AND. (.NOT. config%isodens)) then
      inquire(file=config%database//"/"//config%qc_calc//"/"//config%sac_param_path, exist=ex)
      if (ex) then
         config%database=config%database//"/"//config%qc_calc
      else
         call fatal_error(error,"Database for "//config%qc_calc//" does not seem to be set up properly.")
         return
      end if
   end if

   inquire(file=config%sac_param_path, exist=ex)
   if (.not. ex) then
      inquire(file=config%database//"/"//config%sac_param_path, exist=ex)
      if (ex) then
         call move_line(config%database//"/"//config%sac_param_path, config%sac_param_path)
      else
         inquire(file=home//config%sac_param_path, exist=ex)
         if (ex) then
            call move_line(home//config%sac_param_path, config%sac_param_path)
         else
            Call fatal_error(error, "No crs Parameter File for CPCM-X found.")
         end if
      end if
   end if

   inquire(file=config%smd_param_path, exist=ex)
   if (.not. ex) then
      inquire(file=config%database//"/"//config%smd_param_path, exist=ex)
      if ((ex) .and. (.not. config%smd_default)) then
         call move_line(config%database//"/"//config%smd_param_path, config%smd_param_path)
      else
         inquire(file=home//config%smd_param_path, exist=ex)
         if (ex) then
            call move_line(home//config%smd_param_path, config%smd_param_path)
         else
            Call fatal_error(error, "No smd Parameter File for CPCM-X found.&
            &You can skip this check and use default parameters with the default flag.")
         end if
      end if
   end if
 
    if ((.not.(allocated(config%input))) .AND. (.not. (allocated(config%smd_solvent)))) then
       if (.not.allocated(error)) then
          call help(output_unit)
          stop
       end if
    end if

    if ((.not.(allocated(config%xyz_input))) .AND. (config%qc_calc .eq. "orca")) then
       if (.not.allocated(error)) then
          call help(output_unit)
          stop
       end if
    end if

end subroutine get_arguments 

   !> Subroutine to Read the CPCM-SACMD Input File
subroutine read_input(config,error)
   type(configuration_type) :: config
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
   !> Check if the CPCM-SACMD Input File Exists.
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

subroutine use_default(config, solv, home, error)
   use mctc_env_system, only : get_variable
   use tomlf, only : toml_table, toml_parse, toml_error, toml_key, get_value
   !> Solvent used for default configuration
   character(len=*), intent(in) :: solv
   !> Configuration Type
   type(configuration_type), intent(inout) :: config
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Toml unit
   type(toml_table), allocatable :: config_table
   !> Toml error handling
   type(toml_error), allocatable :: config_error

   character(len=:), allocatable, intent(in) :: home
   
   
   type(toml_key), allocatable, dimension(:) :: list
   integer :: stat
   character(len=255) :: line
   character(:), allocatable :: line2
   logical :: ex
   character(len=10) :: control, command
   integer :: nconf, io

   character(len=:), allocatable :: user
  
   if (solv .eq. '') then
      Call fatal_error(error,'The Solvent you specified is not available.')
      return
   end if
   !> Set Defaults
   config%T=298.15_wp
   config%ML=.FALSE.
   config%sig_in=.FALSE.
   config%prof=.FALSE.
   config%smd_default=.FALSE.
   config%time=.FALSE.
   config%qc_eps=0
   config%probe=0.4
   call move_line(solv,config%smd_solvent)
   call move_line(solv//".cosmo",config%csm_solvent)
   call move_line("solute.cosmo",config%csm_solute)
   call move_line("crs",config%model)

   ex=.false.


   Call get_variable("USER",user)
   inquire(file="/home/"//user//"/cpcmx.toml", exist=ex)

   if (ex) then
      config%config_path="/home/"//user//"/cpcmx.toml"
   else
  
      inquire(file=home//"cpcmx.toml",exist=ex)

      if (.not. ex) then 
         !call fatal_error(error, "No configuration found in "//home)
         return
      else
         config%config_path=home//"cpcmx.toml"
      end if

   end if

   open(file=config%config_path,newunit=io)

   call toml_parse(config_table,io,config_error)
   close(io)
   if (allocated(config_table)) then
      call config_table%get_keys(list)
      do nconf=1,size(list)
         select case(list(nconf)%key)
         case("prog") 
               if (allocated(config%qc_calc)) cycle
               call get_value(config_table,list(nconf),line2)
               if (line2 .eq. "NONE") cycle
               call move_line(to_lower(line2),config%qc_calc)
         case("smd_h2o")
            if (solv .eq. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(line2,config%smd_param_path)
            end if
         case("smd_ot")
            if (solv .ne. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(line2,config%smd_param_path)
            end if
         case("crs_h2o")
            if (solv .eq. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(line2,config%sac_param_path)
            end if
         case("crs_ot")
            if (solv .ne. "water") then
               call get_value(config_table,list(nconf),line2)
               call move_line(line2,config%sac_param_path)
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
   type(configuration_type) :: config

   if (config%internal) then
      config%smd_param_path="internal SMD parameter"
      config%sac_param_path="internal CPCM-X parameter"
      config%csm_solvent="internal COSMO Database File"
   end if


   write(output_unit,'(10x,a)') &
      " ------------------------------------------------- ",&
      "|                 Initialization                  |",&
      " ------------------------------------------------- ", &
      ""

   if (allocated(config%config_path)) write(output_unit,'(5x,a,t35,a)') &
      "Configuration File used:", config%config_path
   write(output_unit,'(5x,a,t35,a)') &
      "SMD Parameter Path:", config%smd_param_path, &
      "CPCM-X Parameter Path:", config%sac_param_path, &
      "Solvent:", config%smd_solvent, &
      "Corresponding COSMO File:", config%csm_solvent

   if (allocated(config%qc_calc)) then 
      write(output_unit,'(5x,a,t35,a)') &
      "QC Program:", config%qc_calc
   else
      write(output_unit,'(5x,a,t35,a)') &
      "Solute COSMO File:", config%csm_solute
   end if
      
end subroutine echo_init

subroutine check_error(error)
   type(error_type), intent(in), allocatable :: error

   integer :: point
   
   if (allocated(error)) then
      write(error_unit,'(a)') ""
      point=index(error%message,".")
      if (point .ne. 0) then
         write(error_unit,'(5x,a,t20,a)') "[ERROR]",trim(error%message(:point))
         write(error_unit,'(5x,a,t20,a)') "",trim(error%message(point+1:))
      else
         write(error_unit,'(5x,a,t20,a)') "[ERROR]",error%message
      end if
      error stop
   end if

end subroutine check_error

end program CPCMX

