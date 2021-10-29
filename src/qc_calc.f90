!> External QC Packages Driver Program

module qc_calc
    use mctc_env, only : wp
    private
    public :: qc_cal
    interface qc_cal
        module procedure :: turbomole
    end interface qc_cal

contains

    !> turbomole subroutine - needs control file for gas phase calculation and coord file in 
    subroutine turbomole(epsilon, cosmo_out, solvent)
        !> Dielectric Constant for COSMO Calculation
        real(wp), intent(inout) :: epsilon
        !> Output File for COSMO Calculation
        character(len=*), intent(in) :: cosmo_out
        !> SMD Solvent for default Epsilon (optional)
        character(len=*), intent(in), optional :: solvent

        !> Control File exists?
        logical :: ex
        !> Reading the control file
        character(len=100) :: line
        !> Charge of molecule from control file
        integer :: charge
        !> Is Molecule charged?
        logical :: ion

        INQUIRE(file='control', exist=ex)

        if (.not. ex) error stop 'Turbomole driver mode specified, but no control file found.'
        if (epsilon .lt. 0) epsilon=minnesota_eps(solvent)
        write(*,*) 'Turbomole Driver initialized.'
        write(*,*) 'Gas phase single point calculation.'
        Call execute_command_line('ridft > gas.out', WAIT=.true.)
        Call execute_command_line('kdg end')
        ion=.false.
        INQUIRE(file='.CHRG', exist=ex)
        if (ex) then
            open(11, file='.CHRG')
            read(11,*) charge
            if (charge .ne. 0) ion=.true.
            close(11)
        end if
        open(11, file='control', access='append')
        write(11,'(A)') '$cosmo'
        if (epsilon .ne. 0) then 
            write(11,'(A11, F0.2, A4)')'   epsilon=',epsilon, merge(' ion','    ',ion)
        else 
            write(11,'(A19, A4)')'   epsilon=infinity', merge(' ion','    ',ion)
        end if
        write(11,'(A16,A)') '$cosmo_out file=',cosmo_out
        write(11,'(A4)') '$end' 
        if (epsilon .ne. 0) then 
            write(*,'(A,F0.2,A),A') 'COSMO Calculation with epsilon=',epsilon, merge(' ion','    ',ion)
        else
            write(*,'(A,A)') 'COSMO Calculation with epsilon=infinity', merge(' ion','    ',ion)
        end if
        Call execute_command_line('ridft > solv.out', WAIT=.true.)
    end subroutine turbomole

    !> Get default dielectric constant from Minnesota Solvation Database
    function minnesota_eps(solvent) result(epsilon)
        character(len=*), intent(in) :: solvent
        real(wp):: epsilon

        select case(solvent)
        case('2methylpyridine')
            epsilon=9.9533_wp
        case('4methyl2pentanone')
            epsilon=12.8871
        case('aceticacid')
            epsilon=6.2528
        case('acetonitrile')
            epsilon=35.6881
        case('acetophenone')
            epsilon=17.44
        case('aniline')
            epsilon=6.8882
        case('anisole')
            epsilon=4.2247
        case('benzene')
            epsilon=2.2706
        case('benzonitrile')
            epsilon=25.592
        case('benzylalcohol')
            epsilon=12.4569
        case('bromobenzene')
            epsilon=5.3954
        case('bromoethane')
            epsilon=9.01
        case('bromoform')
            epsilon=4.2488
        case('bromooctane')
            epsilon=5.0244
        case('butanol')
            epsilon=17.3323
        case('butanone')
            epsilon=18.2457
        case('butylacetate')
            epsilon=4.9941
        case('butylbenzene')
            epsilon=2.36
        case('carbondisulfide')
            epsilon=2.6105
        case('carbontet')
            epsilon=2.228
        case('chlorobenzene')
            epsilon=5.6968
        case('chloroform')
            epsilon=4.7113
        case('chlorohexane')
            epsilon=5.9491
        case('cyclohexane')
            epsilon=2.0165
        case('cyclohexanone')
            epsilon=15.6186
        case('decalin')
            epsilon=2.196
        case('decane')
            epsilon=1.9846
        case('decanol')
            epsilon=7.5305
        case('dibromoethane')
            epsilon=4.9313
        case('dibutylether')
            epsilon=3.0473
        case('dichloroethane')
            epsilon=10.125
        case('diethylether')
            epsilon=4.24
        case('diisopropylether')
            epsilon=3.38
        case('dimethylacetamide')
            epsilon=37.7807
        case('dimethylformamide')
            epsilon=37.219
        case('dimethylpyridine')
            epsilon=7.1735
        case('dimethylsulfoxide')
            epsilon=46.826
        case('dodecane')
            epsilon=2.006
        case('ethanol')
            epsilon=24.852
        case('ethoxybenzene')
            epsilon=4.1797
        case('ethylacetate')
            epsilon=5.9867
        case('ethylbenzene')
            epsilon=2.4339
        case('fluorobenzene')
            epsilon=5.42
        case('fluoroctane')
            epsilon=3.89
        case('heptane')
            epsilon=1.9113
        case('heptanol')
            epsilon=11.321
        case('hexadecane')
            epsilon=2.0402
        case('hexadecyliodide')
            epsilon=3.5338
        case('hexane')
            epsilon=1.8819
        case('hexanol')
            epsilon=12.5102
        case('iodobenzene')
            epsilon=4.547
        case('isobutanol')
            epsilon=16.7766
        case('isooctane')
            epsilon=1.9358
        case('isopropanol')
            epsilon=19.2645
        case('isopropylbenzene')
            epsilon=2.3712
        case('isopropyltoluene')
            epsilon=2.2322
        case('mcresol')
            epsilon=12.44
        case('mesitylene')
            epsilon=2.265
        case('methoxyethanol')
            epsilon=17.2
        case('methylenechloride')
            epsilon=8.93
        case('methylformamide')
            epsilon=181.5619
        case('nitrobenzene')
            epsilon=34.8091
        case('nitroethane')
            epsilon=28.2896
        case('nitromethane')
            epsilon=36.5623
        case('nonane')
            epsilon=1.9605
        case('nonanol')
            epsilon=8.5991
        case('octane')
            epsilon=1.9406
        case('octanol')
            epsilon=9.8629
        case('odichlorobenzene')
            epsilon=9.9949
        case('onitrotoluene')
            epsilon=25.6692
        case('pentadecane')
            epsilon=2.0333
        case('pentane')
            epsilon=1.8371
        case('pentanol')
            epsilon=15.13
        case('perfluorobenzene')
            epsilon=2.029
        case('phenylether')
            epsilon=3.73
        case('propanol')
            epsilon=20.5237
        case('pyridine')
            epsilon=12.9776
        case('secbutanol')
            epsilon=15.9436
        case('secbutylbenzene')
            epsilon=2.3446
        case('tbutylbenzene')
            epsilon=2.3447
        case('tetrachloroethene')
            epsilon=2.268
        case('tetrahydrofuran')
            epsilon=7.4257
        case('tetrahydrothiophenedioxide')
            epsilon=43.9622
        case('tetralin')
            epsilon=2.771
        case('toluene')
            epsilon=2.3741
        case('tributylphosphate')
            epsilon=8.1781
        case('triethylamine')
            epsilon=2.3832
        case('trimethylbenzene')
            epsilon=2.3653
        case('undecane')
            epsilon=1.991
        case('water','h2o')
            epsilon=78.36_wp
        case('xylene')
            epsilon=2.3879
        case('benzene-water')
            epsilon=2.2706
        case('carbontet-water')
            epsilon=2.228
        case('chlorobenzene-water')
            epsilon=5.6968
        case('chloroform-water')
            epsilon=4.7113
        case('cyclohexane-water')
            epsilon=2.0165
        case('dibromoethane-water')
            epsilon=4.9313
        case('dibutylether-water')
            epsilon=3.0473
        case('dichloroethane-water')
            epsilon=10.125
        case('diethylether-water')
            epsilon=4.24
        case('ethylacetate-water')
            epsilon=5.9867
        case('heptane-water')
            epsilon=1.9113
        case('hexane-water')
            epsilon=1.8819
        case('nitrobenzene-water')
            epsilon=34.8091
        case('octanol-water')
            epsilon=9.8629
        case('methanol')
            epsilon=32.613
        end select
    end function minnesota_eps


end module qc_calc
