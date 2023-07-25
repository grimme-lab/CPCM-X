module test_io
    use testdrive, only : new_unittest, error_type, unittest_type, check, test_failed
    use cpx, only: calculation_type, read_cosmo
    use mctc_env, only: wp, err_type => error_type
    implicit none

contains

subroutine collect_iotests(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest("read_cosmo_xtb", test_read_cosmo_xtb), &
        new_unittest("read_cosmo_tm", test_read_cosmo_tm), &
        new_unittest("read_cosmo_orca", test_read_cosmo_orca) &
        ]

end subroutine collect_iotests

subroutine test_read_cosmo_xtb(error)
    type(error_type), allocatable, intent(out) :: error
    type(err_type), allocatable :: err
    type(calculation_type) :: calc

    call read_cosmo('DB/xtb/water.cosmo',calc%solute,'NONE',err)

    if (allocated(err)) Call test_failed(error, err%message)
end subroutine test_read_cosmo_xtb

subroutine test_read_cosmo_tm(error)
    type(error_type), allocatable, intent(out) :: error
    type(err_type), allocatable :: err
    type(calculation_type) :: calc

    call read_cosmo('DB/tm/water.cosmo',calc%solute,'NONE',err)

    if (allocated(err)) Call test_failed(error, err%message)    
end subroutine test_read_cosmo_tm

subroutine test_read_cosmo_orca(error)
    type(error_type), allocatable, intent(out) :: error
    type(err_type), allocatable :: err
    type(calculation_type) :: calc
    
    call read_cosmo('DB/orca/water.cosmo',calc%solute,'NONE',err)
    
    if (allocated(err)) Call test_failed(error, err%message)
end subroutine test_read_cosmo_orca

end module test_io
