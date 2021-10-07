!> External QC Packages Driver Program

module qc_calc
    use mctc_env, only : wp
    interface qc_cal
        module procedure :: turbomole
    end interface qc_cal

contains

    !> turbomole subroutine - needs control file for gas phase calculation and coord file in 
    subroutine turbomole(epsilon, cosmo_out)
        !> Dielectric Constant for COSMO Calculation
        real(wp), intent(in) :: epsilon
        !> Output File for COSMO Calculation
        character(len=*), intent(in) :: cosmo_out

        !> Control File exists?
        logical :: ex

        INQUIRE(file='control', exist=ex)

        if (.not. ex) error stop 'Turbomole driver mode specified, but no control file found.'
        write(*,*) 'Turbomole Driver initialized.'
        write(*,*) 'Gas phase single point calculation.'
        Call execute_command_line('ridft > gas.out', WAIT=.true.)
        Call execute_command_line('kdg end')
        open(11, file='control', access='append')
        write(11,*) '$cosmo'
        if (epsilon .ne. 0) then 
            write(11,*) '   epsilon=',epsilon
        else 
            write(11,*) '   epsilon=infinity'
        end if
        write(11,*) '$cosmo_out file=',cosmo_out
        write(11,*) '$end' 
        if (epsilon .ne. 0) then 
            write(*,*) 'COSMO Calculation with epsilon=',epsilon
        else
            write(*,*) 'COSMO Calculation with epsilon=infinity'
        end if
        Call execute_command_line('ridft > solv.out', WAIT=.true.)
    end subroutine turbomole

end module qc_calc
