module cpxcalc
    use mctc_env, only: wp
    use type, only: molecule_data, parameter_type

    private

    public :: calculation_type

        !> COSMO calculation type
    type calculation_type
        ! COSMO solute type
        type(molecule_data) :: solute
        ! COSMO solvent type
        type(molecule_data) :: solvent
        !> Parameters
        type(parameter_type) :: param 
        !> SMD Parameter file as string array
        character(len=200), dimension(:), allocatable :: smd_param
        !> Temperature of the calculation
        real(wp) :: T
        !> Ideal Screening Energy
        real(wp) :: id_scr
        !> Maximum number of iterations
        integer :: max_cycle
        !> Convergence criterion
        real(wp) :: conv_crit
        !> Chemical potential of the solute in the solvent
        real(wp) :: chem_pot_sol

        !> restoring free energy in Eh
        real(wp) :: dG_res
        !> averaging correction in Eh
        real(wp) :: dG_cc
        !> ideal state energy in Eh
        real(wp) :: dG_is
        !> SMD correction in Eh
        real(wp) :: dG_smd
        !> Standard State correction in Eh
        real(wp) :: dG_ss
        !> Empircal Shift in Eh
        real(wp) :: dG_shift


    contains
        !> Total free energy in Eh
        procedure :: dG
        !> Avearge charge
        procedure :: average_charge => average_charge
        !> Calculate sigma profile
        procedure :: sigma => calc_sigma
        !> Setup bonding situation, hydrogen bond grouping and ring detection
        procedure :: init_bonding
        !> Split sigma profile depending on the bonding situation
        procedure :: split_sigma => sigma_splitting
        !> Calculate solvation energy
        procedure :: solv => calculate_solvation
        procedure :: cds_normal
        procedure :: cds_isodens
        procedure :: cds_default
        !> Calculate CDS part of the energy (SMD)
        generic :: cds => cds_normal, cds_isodens, cds_default
        !> Calculate the state correction
        procedure :: state_correction => calculate_state_correction
        
    end type calculation_type


    contains

    subroutine calculate_state_correction(self,density,mass,T)
        use crs, only: state_correction
        class(calculation_type), intent(inout) :: self
        !> Density of the solvent
        real(wp), intent(in) :: density
        !> Molecular mass of the solvent
        real(wp), intent(in) :: mass
        !> Temperature of the calculation
        real(wp), intent(in) :: T

        Call state_correction(density,mass,T,self%dG_ss)
    end subroutine calculate_state_correction

    pure function dG(self) 
        class(calculation_type), intent(in) :: self
        real(wp) :: dG

        !> sum(dG_res,dG_cc,dG_is,dG_smd,dG_ss,dG_shift)
        dG=self%dG_res+self%dG_cc+self%dG_is+self%dG_smd+self%dG_ss+self%dG_shift

    end function dG

    subroutine cds_default(self,probe,smd_solvent)
        use sdm, only: calculate_cds
        use sort, only: unique
        class(calculation_type), intent(inout) :: self
        !> Probe radius
        real(wp), intent(in) :: probe
        !> SMD solvent
        character(len=*), intent(in) :: smd_solvent


        Call calculate_cds(unique(self%solute%id),self%solute%element,self%solute%atom_xyz,probe,smd_solvent,&
        &self%dG_smd,self%smd_param)
    end subroutine cds_default

    subroutine cds_normal(self,probe,smd_solvent,smd_param_path,smd_default)
        use sdm, only: calculate_cds
        use sort, only: unique
        class(calculation_type), intent(inout) :: self
        !> Probe radius
        real(wp), intent(in) :: probe
        !> SMD solvent
        character(len=*), intent(in) :: smd_solvent
        !> SMD parameter path
        character(len=*), intent(in) :: smd_param_path
        !> Use default SMD parameters
        logical, intent(in), optional :: smd_default


        Call calculate_cds(unique(self%solute%id),self%solute%element,self%solute%atom_xyz,probe,smd_solvent,&
        &smd_param_path,self%dG_smd,smd_default)


    end subroutine cds_normal

    subroutine cds_isodens(self,probe,smd_solvent,smd_param_path,isodens_rad)
        use sdm, only: calculate_cds
        use sort, only: unique
        class(calculation_type), intent(inout) :: self
        !> Probe radius
        real(wp), intent(in) :: probe
        !> SMD solvent
        character(len=*), intent(in) :: smd_solvent
        !> SMD parameter path
        character(len=*), intent(in) :: smd_param_path
        !> Isodensity radius
        real(wp), dimension(:), allocatable, intent(in) :: isodens_rad


        Call calculate_cds(unique(self%solute%id),self%solute%element,self%solute%atom_xyz,probe,smd_solvent,&
        &smd_param_path,self%dG_smd,isodens_rad)

    end subroutine cds_isodens

    subroutine calculate_solvation(self,model,error,T,max_cycle,conv_crit)
        use mctc_env, only: error_type, fatal_error
        use globals, only: autokcal
        class(calculation_type), intent(inout) :: self
        !> Model to use ("cpcmx,sac,sac2010")
        character(len=*), intent(in) :: model
        !> Error handling
        type(error_type), intent(out), allocatable :: error
        !> Temperature
        real(wp), intent(in), optional :: T
        !> Maximum number of iterations
        integer, intent(in), optional :: max_cycle
        !> Convergence criterion
        real(wp), intent(in), optional :: conv_crit

        if (present(T)) then
            self%T=T
        else
            self%T=298.15_wp
        end if

        if (present(max_cycle)) then
            self%max_cycle=max_cycle
        else
            self%max_cycle=500
        end if

        if (present(conv_crit)) then
            self%conv_crit=conv_crit
        else
            self%conv_crit=1.0e-4_wp
        end if

        select case(model)
            case("cpcmx", "crs")
                call calculate_cpcmx(self,error)
            case("sac")
                Call fatal_error(error, "SAC model is not implemented yet.")
                return
            case("sac2010")
                Call fatal_error(error, "SAC2010 model is not implemented yet.")
                return
            case default
                call fatal_error(error,"Unknown model: "//model)
                return
        end select

        self%dG_res=(self%chem_pot_sol+self%param%omega*self%solute%near)/autokcal

    end subroutine calculate_solvation

    subroutine calculate_cpcmx(self,error)
        use mctc_env, only: error_type, fatal_error
        use data, only: density, AtomicMass
        use crs, only: calcgas, compute_solute, compute_solvent, state_correction
        class(calculation_type), intent(inout) :: self
        !> Error handling
        type(error_type), intent(out), allocatable :: error

        Call calcgas(self%solute%energy,self%solute%energy_gas,self%id_scr,self%solute%area,self%solute%sv,&
        &self%solute%su,self%solute%pot,self%solute%element,self%solute%id,self%T,self%param%eta,self%dG_is,self%dG_cc)
        Call compute_solvent(self%param,self%solvent%poti,self%solvent%sv,self%solvent%svt,&
        &self%solvent%area,self%T,self%max_cycle,self%conv_crit,error)
        if (allocated(error)) return
        Call compute_solute(self%param,self%solute%poti,self%solvent%poti,self%solute%sv,self%solute%svt,&
        &self%solvent%sv,self%solvent%svt,self%solute%area,self%solvent%area,self%T,self%chem_pot_sol)


    end subroutine calculate_cpcmx

    subroutine sigma_splitting(self,writesigma)
        use profile, only: split_sigma
        class(calculation_type), intent(inout) :: self
        logical, intent(in), optional :: writesigma

        logical :: sw

        if (present(writesigma)) then
            sw=writesigma
        else
            sw=.false.
        end if

        if (sw) then
            call split_sigma(self%solvent%sv, self%solvent%area, self%solvent%hb_group,self%solvent%id&
            &,self%solvent%element,self%solvent%sigma3, "solvent")
            call split_sigma(self%solute%sv, self%solute%area, self%solute%hb_group,self%solute%id&
            &,self%solute%element,self%solute%sigma3, "solute")
        else
            call split_sigma(self%solvent%sv, self%solvent%area, self%solvent%hb_group,self%solvent%id&
            &,self%solvent%element,self%solvent%sigma3)
            call split_sigma(self%solute%sv, self%solute%area, self%solute%hb_group,self%solute%id&
            &,self%solute%element,self%solute%sigma3)
        end if

    end subroutine sigma_splitting

    subroutine init_bonding(self)
        use bonding, only: det_bonds, hb_grouping, det_rings

        class(calculation_type), intent(inout) :: self

        Call det_bonds(self%solute%id,self%solute%atom_xyz,self%solute%element,self%solute%is_bonded&
        &,self%solute%noh,self%solute%nnh)
        Call det_bonds(self%solvent%id,self%solvent%atom_xyz,self%solvent%element,self%solvent%is_bonded&
        &,self%solvent%noh,self%solvent%nnh)
        Call hb_grouping(self%solute%id,self%solute%element,self%solute%is_bonded,self%solute%hb_group)
        Call hb_grouping(self%solvent%id,self%solvent%element,self%solvent%is_bonded,self%solvent%hb_group)
        Call det_rings(self%solute%id,self%solute%is_bonded,self%solute%ring,self%solute%near)


    end subroutine init_bonding

    subroutine average_charge(self, error)
        use sigma_av, only: average_one_charge, ortho_charge
        use mctc_env, only: error_type, fatal_error

        real(wp) :: rav
        class(calculation_type), intent(inout) :: self
        type(error_type), intent(out), allocatable :: error

        if (self%param%rav==0.0_wp) then
            call fatal_error(error,"rav is zero")
            return
        end if
        rav=self%param%rav

        call average_one_charge(rav, self%solvent%xyz,self%solvent%su,self%solvent%area,self%solvent%sv)
        call average_one_charge(rav, self%solute%xyz, self%solute%su, self%solute%area, self%solute%sv)
        call average_one_charge(rav*2.0_wp, self%solvent%xyz,self%solvent%su,self%solvent%area,self%solvent%sv0)
        call average_one_charge(rav*2.0_wp, self%solute%xyz, self%solute%su, self%solute%area, self%solute%sv0)
        Call ortho_charge(self%solvent%sv, self%solvent%sv0, self%solvent%svt)
        Call ortho_charge(self%solute%sv, self%solute%sv0, self%solute%svt)

    end subroutine average_charge

    subroutine calc_sigma(self,writesigma)
        use profile, only: single_sigma

        !> Calculation type
        class(calculation_type), intent(inout) :: self

        !> Write output files for solvent and solute
        logical, intent(in),optional :: writesigma

        logical :: sw

        if (present(writesigma)) then
            sw=writesigma
        else
            sw=.false.
        end if 

        if (sw) then
            call single_sigma(self%solvent%sv, self%solvent%area, self%solvent%sigma, "solvent")
            call single_sigma(self%solute%sv, self%solute%area, self%solute%sigma, "solute")
        else
            call single_sigma(self%solvent%sv, self%solvent%area, self%solvent%sigma)
            call single_sigma(self%solute%sv, self%solute%area, self%solute%sigma)
        end if

    end subroutine calc_sigma


end module cpxcalc