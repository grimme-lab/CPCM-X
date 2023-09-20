module cpx_c_api
   use, intrinsic :: iso_c_binding
   use mctc_env, only: wp
   use globals, only: to_lower
   use cpx

   type :: vcalc_type
      class(calculation_type), allocatable :: ptr
   end type vcalc_type

contains

   function newCalculation_api() result(vcalc) &
      & bind(c, name="cpx_newCalculation")
   type(vcalc_type), pointer :: calc
   type(c_ptr) :: vcalc

   allocate(calc)
   allocate(calc%ptr)
   vcalc=c_loc(calc)

   end function newCalculation_api

   subroutine deleteCalculation_api(vcalc) &
      & bind(c, name="cpx_deleteCalculation")
   type(c_ptr), intent(inout) :: vcalc
   type(vcalc_type), pointer :: calc

   if (c_associated(vcalc)) then
      call c_f_pointer(vcalc, calc)
      deallocate(calc%ptr)
      deallocate(calc)
      vcalc=c_null_ptr
   end if

   end subroutine deleteCalculation_api

   subroutine read_cosmo_file_api(vcalc,vsolvsol,vcosmo_file) &
      & bind(c, name="cpx_readCosmoFile")
   use mctc_env, only: error_type
   character(kind=c_char), dimension(*), intent(in) :: vcosmo_file
   character(kind=c_char), dimension(*), intent(in) :: vsolvsol
   type(c_ptr), intent(inout) :: vcalc
   type(vcalc_type), pointer :: calc
   type(error_type), allocatable :: error
   character(len=:, kind=c_char), allocatable :: cosmo_file
   character(len=:, kind=c_char), allocatable :: solvsol

   call c_f_character(vcosmo_file, cosmo_file)
   call c_f_character(vsolvsol, solvsol)
   call c_f_pointer(vcalc, calc)
   select case(trim(to_lower(solvsol)))
         case("solvent")
            call read_cosmo(cosmo_file, calc%ptr%solvent, 'NONE', error)
         case("solute")
            call read_cosmo(cosmo_file, calc%ptr%solute, 'NONE', error)
         case default
            vcalc=c_null_ptr
            return
   end select

   if (allocated(error)) then
      vcalc=c_null_ptr
   else
      vcalc=c_loc(calc)
   end if

   end subroutine read_cosmo_file_api


   !> Load parameter from internal database
   subroutine load_param_api(vmethod,vsolvent,vcalc) &
      & bind(c, name="cpx_loadparam")
   use mctc_env, only: error_type
   character(kind=c_char), dimension(*), intent(in) :: vmethod
   character(kind=c_char), dimension(*), intent(in) :: vsolvent
   character(len=:, kind=c_char), allocatable :: method, solvent
   type(c_ptr), intent(inout) :: vcalc
   type(vcalc_type), pointer :: calc
   type(calculation_type), allocatable :: dummy_calc
   type(error_type), allocatable :: error

   call c_f_character(vmethod, method)
   call c_f_character(vsolvent, solvent)
   call c_f_pointer(vcalc, calc)

   allocate(dummy_calc)
  
   call initialize_param(method, solvent, dummy_calc, error)
   if (allocated(error)) then
      vcalc=c_null_ptr
   else
      call move_alloc(dummy_calc,calc%ptr)
      vcalc=c_loc(calc)
   end if

   end subroutine load_param_api

   !> Read parameter from file
   subroutine read_param_api(vparam_file_crs,vparam_file_smd,vcalc) &
      & bind(c, name="cpx_readparam")
   use mctc_env, only: error_type, fatal_error
   character(kind=c_char), dimension(*), intent(in) :: vparam_file_crs
   character(kind=c_char), dimension(*), intent(in) :: vparam_file_smd
   character(len=:, kind=c_char), allocatable :: method, solvent, param_file_crs, param_file_smd
   type(c_ptr), intent(inout) :: vcalc
   type(vcalc_type), pointer :: calc
   type(calculation_type), allocatable :: dummy_calc
   type(error_type), allocatable :: error

   logical :: e
   integer :: nu, line, estat, lines

   call c_f_character(vparam_file_crs, param_file_crs)
   call c_f_character(vparam_file_smd, param_file_smd)
   call c_f_pointer(vcalc, calc)


   allocate(dummy_calc)

   !! Initialize crs_parameters
   call initialize_param(param_file_crs, dummy_calc%param,error)
   inquire(file=trim(param_file_smd), exist=e)
   if (e) then
      !! Initialize smd_parameters
      open(newunit=nu, file=param_file_smd, status='old', action='read')
      lines=0
      estat=0
      do
         read(nu,'(a)',iostat=estat)
         if (estat/=0) exit
         lines=lines+1
      end do
      rewind(nu)
      allocate(dummy_calc%smd_param(lines))
      do line=1,lines
         read(nu,'(a)') dummy_calc%smd_param(line)
      end do
      close(nu)
   else
      call fatal_error(error,"SMD parameter file does not exist")
   end if

   if (allocated(error)) then
      vcalc=c_null_ptr
   else
      call move_alloc(dummy_calc,calc%ptr)
      vcalc=c_loc(calc)
   end if

   end subroutine read_param_api

   

   !> Load solvent from internal database
   subroutine load_solvent_api(vsolvent,vcalc) &
      & bind(c, name="cpx_loadsolvent")
   use mctc_env, only: error_type
   character(kind=c_char), dimension(*), intent(in) :: vsolvent
   character(len=:, kind=c_char), allocatable :: solvent
   type(c_ptr), intent(inout) :: vcalc
   type(vcalc_type), pointer :: calc
   type(error_type), allocatable :: error

   call c_f_character(vsolvent, solvent)
   call c_f_pointer(vcalc,calc)
   call load_solvent(solvent, calc%ptr%solvent, error)
   if (allocated(error)) then
      vcalc=c_null_ptr
   else
      vcalc=c_loc(calc)
   end if

   end subroutine load_solvent_api

   !> Load solute from Database (only for testing purposes)
   subroutine load_solute_api(vsolute,vcalc) &
      & bind(c, name="cpx_loadsolute")
   use mctc_env, only: error_type
   character(kind=c_char), dimension(*), intent(in) :: vsolute
   character(len=:, kind=c_char), allocatable :: solute
   type(c_ptr), intent(out) :: vcalc
   type(vcalc_type), pointer :: calc
   type(error_type), allocatable :: error

   call c_f_character(vsolute, solute)
   call c_f_pointer(vcalc,calc)
   call load_solvent(solute, calc%ptr%solute, error)
   if (allocated(error)) then
      vcalc=c_null_ptr
   else
      vcalc=c_loc(calc)
   end if

   end subroutine load_solute_api

   subroutine calculate_cpcmx_api(vcalc,vmethod,vsolvent,egas,probe,T,conv_crit) &
      & bind(c, name="cpx_calculate")
   use mctc_env, only: error_type
   character(kind=c_char), dimension(*), intent(in) :: vmethod
   character(kind=c_char), dimension(*), intent(in) :: vsolvent
   real(c_double), intent(in), value :: egas
   type(c_ptr),intent(inout) :: vcalc
   real(c_double), intent(in), value :: probe
   real(c_double), intent(in), value :: T
   real(c_double), intent(in), value :: conv_crit
   type(vcalc_type), pointer :: calc
   type(error_type), allocatable :: error
   character(len=:, kind=c_char), allocatable :: method, solvent

   call c_f_character(vmethod, method)
   call c_f_character(vsolvent, solvent)

   if (c_associated(vcalc)) then
      call c_f_pointer(vcalc, calc)
      calc%ptr%solute%energy_gas=egas


      call calc%ptr%average_charge(error)
      call calc%ptr%init_bonding()
      call calc%ptr%solv(method,error,T,500,conv_crit)
      call calc%ptr%state_correction(density(solvent),atomicmass(calc%ptr%solvent%element),T)
      call calc%ptr%cds(probe,solvent)

      vcalc=c_loc(calc)
   else
      vcalc=c_null_ptr
      return
   end if

   
   end subroutine calculate_cpcmx_api


   subroutine c_f_character(rhs, lhs)
      character(kind=c_char), intent(in) :: rhs(*)
      !> Resulting Fortran string
      character(len=:, kind=c_char), allocatable, intent(out) :: lhs
   
      integer :: ii
   
      do ii = 1, huge(ii) - 1
         if (rhs(ii) == c_null_char) then
            exit
         end if
      end do
      allocate(character(len=ii-1) :: lhs)
      lhs = transfer(rhs(1:ii-1), lhs)
   
   end subroutine c_f_character

   subroutine get_energies_api(vcalc,energies) &
      & bind(c, name="cpx_getenergies")
   use mctc_env, only: wp
   type(c_ptr), intent(in) :: vcalc
   type(vcalc_type), pointer :: calc
   real(c_double), dimension(6), intent(out) :: energies

   energies=0.0_wp
   call c_f_pointer(vcalc, calc)
   energies(1) = calc%ptr%dG_is
   energies(2) = calc%ptr%dG_cc
   energies(3) = calc%ptr%dG_res
   energies(4) = calc%ptr%dG_smd
   energies(5) = calc%ptr%dG_ss
   energies(6) = calc%ptr%dG_shift

   end subroutine get_energies_api


end module cpx_c_api

