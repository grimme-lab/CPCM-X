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

!> External QC Packages Driver Program

module qc_calc
    use mctc_env, only : wp, error_type, fatal_error
    use, intrinsic :: iso_fortran_env, only : output_unit, file_storage_size

    implicit none
    private
    public :: qc_cal
    interface qc_cal
        module procedure :: turbomole
        module procedure :: orca
    end interface qc_cal

contains

    !> Orca subroutine, only supports eps=infinity, needs .xyz input file
    subroutine orca(input,error,new_functional,new_basis)
        !> Input Coordinates for Solute
        character(len=*), intent(in) :: input
        !> Error Handling
        type(error_type), allocatable, intent(out) :: error
        !> Changing default functional/basis
        character(len=*), intent(inout), optional :: new_functional, new_basis

        !> Default Functional, Basis
        character(len=:), allocatable :: functional, basis
        !> Charge and Multiplicity
        integer :: charge, multi
        !> Solute is Input without .xyz
        character(:), allocatable :: solute
        !> Input file existing?
        logical :: ex
        !> Reading Output
        character(len=100) :: line
        character(len=25) :: words, dummy1, dummy2, dummy3, dummy4
        real(wp) :: E_gas, E_solv
        !> I/O error
        integer :: io_error, file_size

        !> Default Functional is r2scan-3c
        if (present(new_basis)) then
            functional=new_functional
            basis=new_basis
        else
            functional="r2scan-3c"
            basis=""
        end if

        allocate(character((len(input)-4)) :: solute)
        solute=input(1:len(input)-4)

        INQUIRE(file=input, exist=ex)

        if (.not. ex) then
            Call fatal_error(error,"Input File for Orca Driver Mode not found.")
            return
        end if
        if (input(len(input)-3:len(input)) .ne. ".xyz") then
            Call fatal_error(error,"Orca Driver Mode chosen, but input file does not look like a .xyz file.")
            return
        end if
        write(output_unit,'(a)') ""
        write(output_unit,'(10x,a)') &
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
            " ------------------------------------------------- ",&
            "|                   ORCA Driver                   |",&
            " ------------------------------------------------- "
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!

        charge=0
        multi=1
        
        INQUIRE(file=".CHRG",exist=ex)
        if (ex) then
            open(11,file=".CHRG")
            read(11) charge
            close(11)
        end if

        INQUIRE(file=".UHF",exist=ex)
        if (ex) then
            open(11,file=".UHF")
            read(11) multi
            close(11)
        end if

        INQUIRE(file=solute//".inp",exist=ex)

        if (ex) then
            write(output_unit,'(a)') ""
            write(output_unit,'(5x,a, t20, a)') &
                "[WARNING]", "Found a orca input file in the working directory.", &
                "","This Driver does not support custom input files.", &
                "","Overwriting the found Input. Check your results."
            write(output_unit,'(a)') ""
        end if

        write(output_unit,'(5x,a,3x,a)') "Functional for Orca QC calculations:", functional
        if (basis .ne. "") write(output_unit,'(5x,a,3x,a)')&
                            "Basisset for Orca QC calculations:", basis
        !> Writing an gas phase Orca File
        open(11,file=solute//".inp",status="unknown")
        write(11,'(A1,1x,A,1x,A)') "!", functional, basis
        write(11,'(a)') "!DEFGRID3"
        write(11,'(A1,1x,A7,1x,I1,1x,I1,1x,A)') "*","xyzfile",charge,multi,input
        close(11)

        write(output_unit,'(5x, A)') 'Starting gas phase single point calculation.'
        Call execute_command_line("orca "//solute//".inp > gas.out 2>error.out",WAIT=.true.)
        INQUIRE(file="error.out",size=file_size)
        if (file_size .gt. 0) then
            Call fatal_error(error,"Orca signals an error in the gas phase calculation."&
            &//new_line('a')//"Check error.out for more information.")
            return
        end if
        write(output_unit,'(5x, A)') &
            "Done! Gas phase calculation terminated normally.", &
            ""
        
        open(11,file="gas.out")
        do while (.true.)
            read(11,'(a)') line
            read(line,'(a25)') words
            if (words .eq. "FINAL SINGLE POINT ENERGY") then
                read(line,*) dummy1,dummy2,dummy3,dummy4,E_gas
                exit
            end if
        end do
        close(11)
        open(11,file="gas.energy")
        write(11,*) E_gas
        close(11)


        !> Writing an COSMO Infinity Orca File
        open(11,file=solute//".inp",status="OLD")
        write(11,'(A1,1x,A,1x,A)') "!", functional, basis
        write(11,'(a)') "!DEFGRID3"
        write(11,'(A6)') "! CPCM"
        write(11,'(A1,1x,A7,1x,I1,1x,I1,1x,A)') "*","xyzfile",charge,multi,input
        close(11)

        write(output_unit,'(5x, A)') 'Starting CPCM calculation with epsilon=infinity.'
        Call execute_command_line("orca "//solute//".inp > solvent.out 2>error.out",WAIT=.true.)
        INQUIRE(file="error.out",size=file_size)
        if (file_size .gt. 0) then
            Call fatal_error(error,"Orca signals an error in the solvent phase calculation."&
            &//new_line('a')//"Check error.out for more information.")
            return
        end if
        write(output_unit,'(5x, A)') &
            "Done! CPCM calculation terminated normally.", &
            ""
        open(11,file="solvent.out")
        do while (.true.)
            read(11,'(a)') line
            read(line,'(a25)') words
            if (words .eq. "FINAL SINGLE POINT ENERGY") then
                read(line,*) dummy1,dummy2,dummy3,dummy4,E_solv
                exit
            end if
        end do
        close(11)
        Call orcatocosmo(solute, E_solv)
    end subroutine orca
        

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
        !> Catches some Errors.
        integer :: io_error
        !> Necessary for reading gas phase energy.
        integer :: lines, i
        real(wp) :: E_gas

        INQUIRE(file='control', exist=ex)
        if (.not. ex) error stop 'Turbomole driver mode specified, but no control file found.'
        if (epsilon .lt. 0) epsilon=minnesota_eps(solvent)
        write(output_unit,'(a)') ""
        write(output_unit,'(10x,a)') &
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
            " ------------------------------------------------- ",&
            "|                TURBOMOLE Driver                 |",&
            " ------------------------------------------------- "
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
        io_error=0
        open(11, file="control" ,status='old')
        do while (io_error .ge. 0)
            read(11,'(A)',iostat=io_error) line
            if (index(line,"cosmo") .ne. 0) then
                write(output_unit,'(a)') ""
                write(output_unit,'(5x,a, t20, a)') &
                    "[WARNING]", "Found a COSMO command in the control file.", &
                    "","Deleting the whole COSMO Block for the Gas Phase calculation.", &
                    "","Carefully check your results."
                write(output_unit,'(a)') ""
                Call execute_command_line('kdg cosmo', WAIT=.true.)
                exit
            end if
        end do
        close(11)
        write(output_unit,'(5x, A)') 'Starting Gas phase single point calculation.'
        Call execute_command_line('ridft > gas.out 2>error', WAIT=.true.)
        open(11,file="error")
        read(11,'(a)') line
        if (index(line,"abnormal") .ne. 0) error stop "Gas phase calculation stopped abnormally."
        close(11)

        write(output_unit,'(5x, A)') &
            "Done! Gas phase calculation terminated normally.", &
            ""

        io_error=0
        open(11, file="energy", status="OLD")
        if (io_error .ne. 0) error stop 'Error while reading the gas phase energy. Check your control file.'
        io_error=0
        lines=0
        do while (.TRUE.)
            read(11,*,iostat=io_error)
            if (io_error .lt. 0) exit
            lines=lines+1
        end do
        rewind(11)
        do i=1,lines-2
            read(11,*)
        end do
        read(11,*) i, E_gas
        close(11)
        open(11,file="gas.energy")
        write(11,*) E_gas
        close(11)

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
            write(output_unit,'(5x,A,F0.2,A),A') 'Starting COSMO Calculation with epsilon=',epsilon, merge(' ion','    ',ion)
        else
            write(output_unit,'(5x,A,A)') 'Starting COSMO Calculation with epsilon=infinity', merge(' ion','    ',ion)
        end if
        Call execute_command_line('ridft > solv.out 2>error', WAIT=.true.)
        if (index(line,"abnormal") .ne. 0) error stop "Gas phase calculation stopped abnormally."
        write(output_unit,'(5x, A)') &
            "Done! COSMO phase calculation terminated normally.", &
            ""

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



    subroutine orcatocosmo(oc_inp,energy)
        use globals, only: BtoA
        !> CPCM File (without .cpcm)
        character(len=:),intent(in), allocatable :: oc_inp
        !> CPCM Solvent Energy
        real(wp), intent(in) :: energy

        integer :: segments, atoms
        character(len=100) :: line
        character(len=2), allocatable :: symbols(:)
        real(wp), allocatable :: xyz(:,:), segment_info(:,:)
        real(wp) :: dummy1, dummy2, dummy3, volume, area

        integer :: i


        open(11,file=oc_inp//".xyz")
        read(11,*) atoms
        allocate(symbols(atoms))
        allocate(xyz(atoms,3))
        read(11,*)
        do i=1,atoms
            read(11,*) symbols(i), xyz(i,1), xyz(i,2), xyz(i,3)
        end do
        xyz = xyz/BtoA
        close(11)
        open(11,file=oc_inp//".cpcm")
        read(11,*)
        read(11,*) segments
        allocate(segment_info(segments,9))
        do while (.true.)
            read(11,'(a)') line 
            if (index(line,"Volume") .ne. 0) read(line,*) volume
            if (index(line,"Area") .ne. 0) read(line,*) area
            if (index(line,"SURFACE POINTS (A.U.)") .ne. 0) exit
        end do
        read(11,*)
        read(11,*)
        do i=1,segments
            segment_info(i,1)=i
            read(11,*) segment_info(i,3), segment_info(i,4), segment_info(i,5),&
            & segment_info(i,7),segment_info(i,9),segment_info(i,6), dummy1, dummy2, &
            & dummy3, segment_info(i,2)
            segment_info(i,2)=segment_info(i,2)+1
            segment_info(i,8)=segment_info(i,6)/segment_info(i,7)
        end do
        close(11)

        open(11,file="solute.cosmo")
        write(11,'(a)') "#This is a custom COSMO File generated by CPCM-X"
        write(11,*) "area=",area
        write(11,*) "volume=",volume
        write(11,'(a)') "#atom"
        do i=1,atoms
            write(11,*) i, xyz(i,1), xyz(i,2), xyz(i,3), symbols(i)
        end do
        write(11,'(a)') "$cosmo_energy"
        write(11,'(a)') "#CPCM Singlepoint Energy taken from Orca"
        write(11,'(a,F15.9)') "Final CPCM Single Point energy [a.u.] = ",energy
        write(11,'(a)') ""
        write(11,'(a)') "$segment_information"
        write(11,'(a)') "#Reordered segment information from "//oc_inp//".cpcm"
        do i=1,segments
            write(11,('(I5,3x,I3,3x,F15.9,3x,F15.9,3x,&
            &F15.9,3x,F15.9,3x,F15.9,3x,F15.9,3x,F15.9)'))&
            & int(segment_info(i,1)), int(segment_info(i,2)), segment_info(i,3),&
            & segment_info(i,4), segment_info(i,5), segment_info(i,6), segment_info(i,7),&
            & segment_info(i,8), segment_info(i,9)
        end do

    end subroutine orcatocosmo


end module qc_calc  