module test_calc
    use testdrive, only : new_unittest, error_type, unittest_type, check_ => check, test_failed
    use cpx, only: calculation_type, read_cosmo, initialize_param, load_solvent, density, atomicmass
    use mctc_env, only: wp, err_type => error_type
    implicit none

    real(wp), parameter :: thr = 1e-6_wp

contains

subroutine collect_calctests(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest("xtb_internal_aq", test_calculation_aqueos_xtb_internal), &
        new_unittest("xtb_internal_nonaq", test_calculation_nonaqueos_xtb_internal) &
        ]
end subroutine collect_calctests

subroutine test_calculation_aqueos_xtb_internal(error)
    type(error_type), allocatable, intent(out) :: error
    type(err_type), allocatable :: err
    type(calculation_type) :: calc

    character(len=*), parameter :: solvent='water'

    call initialize_param('xtb', solvent, calc, err)
    if (allocated(err)) call test_failed(error, err%message)
    call load_solvent(solvent,calc%solvent,err)
    if (allocated(err)) call test_failed(error, err%message)

    call read_cosmo('DB/xtb/nitroethane.cosmo',calc%solute,'NONE',err)
    if (allocated(err)) call test_failed(error, err%message)

    call calc%average_charge(err)
    if (allocated(err)) call test_failed(error, err%message)
    call calc%init_bonding()

    call calc%solv('crs',err,298.15_wp,500,0.0001_wp)
    call calc%state_correction(density(solvent),atomicmass(calc%solvent%element),298.15_wp)
    call calc%cds(0.3_wp,solvent)
    
    call check_(error,calc%dG_cc, -.1888025915125338E-03_wp,thr=thr)
    call check_(error,calc%dG_res, .9739686917154823E-02_wp,thr=thr)
    call check_(error,calc%dG_ss, -.6821274273689155E-02_wp,thr=thr)
    call check_(error,calc%dG_smd, -.2244166126383620E-02_wp,thr=thr)
    

end subroutine test_calculation_aqueos_xtb_internal

subroutine test_calculation_nonaqueos_xtb_internal(error)
    type(error_type), allocatable, intent(out) :: error
    type(err_type), allocatable :: err
    type(calculation_type) :: calc

    character(len=*), parameter :: solvent='acetonitrile'

    call initialize_param('xtb', solvent, calc, err)
    if (allocated(err)) call test_failed(error, err%message)
    call load_solvent(solvent,calc%solvent,err)
    if (allocated(err)) call test_failed(error, err%message)

    call read_cosmo('DB/xtb/nitroethane.cosmo',calc%solute,'NONE',err)
    if (allocated(err)) call test_failed(error, err%message)

    call calc%average_charge(err)
    if (allocated(err)) call test_failed(error, err%message)
    call calc%init_bonding()

    call calc%solv('crs',err,298.15_wp,500,0.0001_wp)
    call calc%state_correction(density(solvent),atomicmass(calc%solvent%element),298.15_wp)
    call calc%cds(0.3_wp,solvent)
    
    call check_(error,calc%dG_cc, -.4007613830056836E-03_wp,thr=thr)
    call check_(error,calc%dG_res, .2635940258746455E-02_wp,thr=thr)
    call check_(error,calc%dG_ss, -.5819412185087306E-02_wp,thr=thr)
    call check_(error,calc%dG_smd, .1061086051245725E-02_wp,thr=thr)
    

end subroutine test_calculation_nonaqueos_xtb_internal

end module test_calc

    

