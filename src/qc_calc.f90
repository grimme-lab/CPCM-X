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
    use mctc_io, only: to_number
    use data, only: minnesota_eps
    use numsa, only: get_vdw_rad_cosmo
    use globals, only: rename, autokcal, to_lower, BtoA
    use, intrinsic :: iso_fortran_env, only : output_unit, file_storage_size

    implicit none
    private
    public :: qc_cal, orcatocosmo
    interface qc_cal
        module procedure :: turbomole
        module procedure :: orca
        module procedure :: gtb
        module procedure :: xtb
    end interface qc_cal

contains

    !> turbomole subroutine - needs control file for gas phase calculation and coord file in
    subroutine xtb(solvent, level, error)
        !> Dielectric Constant for COSMO Calculation
        character(len=*), intent(in) :: solvent
        !> Error Handling
        type(error_type), intent(out), allocatable :: error
        !> xtb level used
        character(len=*), intent(in) :: level

        !> File handling
        logical :: ex
        integer :: io_error
        
        !> Charge and multiplicity of molecule (from control file)
        integer :: charge, multi
        !> Necessary for reading gas phase energy.
        character(len=200) :: line
        character(len=:), allocatable :: xtb_bin
        integer :: lines, i
        real(wp) :: E_gas

        charge=0
        multi=0

        xtb_bin='xtb_dev'

        write(output_unit,'(a)') ""
        write(output_unit,'(10x,a)') &
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
            " ------------------------------------------------- ",&
            "|                   XTB Driver                    |",&
            " ------------------------------------------------- "
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
        write(output_unit,'(a)') ""
        io_error=0

        INQUIRE(file='coord', exist=ex)
        if (.not. ex) then
            Call fatal_error(error,'There is no coord file in the working directory.')
            write(output_unit,'(a)') ""
            return
        end if

        Call execute_command_line(xtb_bin//" coord --gfn "//level//" --norestart > gas.out 2>error", WAIT=.true.)

        write(output_unit,'(5x, A)') &
            "Done! Gas phase calculation terminated normally.", &
            ""

        io_error=0
        open(11, file="gas.out", status="OLD")
        if (io_error .ne. 0) then
            Call fatal_error(error,'Error while reading the gas phase energy.')
            return
        end if
        io_error=0
        do while (.TRUE.)
            read(11,'(a)',iostat=io_error) line
            if (io_error .ne. 0) then
                Call fatal_error(error,'Could not find gas phase energy in gas.out.')
                return
            end if
            if (index(line,'TOTAL ENERGY') .ne. 0) exit
        end do

        read(line(index(line,'TOTAL ENERGY',back=.true.)+13:),*) E_gas
        close(11)
        open(11,file="gas.energy")
        write(11,*) E_gas
        close(11)

        Call execute_command_line(xtb_bin//" coord --gfn "//level//" --cosmo "//solvent//" --norestart &
        &> solv.out 2>error", WAIT=.true.)

        write(output_unit,'(5x, A)') &
            "Done! Solvent phase calculation terminated normally.", &
            ""

        call rename('xtb.cosmo','solute.cosmo',io_error)

    end subroutine xtb    

    subroutine gtb(damp,scale,error)

        !> Damping for SCF
        real(wp), intent(in) :: damp

        !> Extrapolation scaling for SCF Energies
        real(wp), intent(in) :: scale

        !> Error Handling
        type(error_type), allocatable, intent(out) :: error

        !> File Handling
        logical :: ex
        integer :: io_error

        !> Charge Handling
        integer :: charge

        logical :: ion

        !> Determine UHF
        integer :: nocc, nalpha,nbeta, nopen

        !> Energy reading
        real(wp), dimension(2) :: energies
        integer :: dummyi, i

        character(len=100) :: line, dummy

        !> Check if coord file is in directory
        INQUIRE(file="coord",exist=ex)

        if (.not. ex) then
            Call fatal_error(error,'gTB mode requested but no coord file in working directory.') 
            return
        end if

        !> Clean up, if there are already mos
        open(5741,file="mos",iostat=io_error,status="old")
        if (io_error .eq. 0) close(5741, status="delete")

        write(output_unit,'(a)') ""
        write(output_unit,'(10x,a)') &
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
            " ------------------------------------------------- ",&
            "|                   P-gTB Driver                   |",&
            " ------------------------------------------------- "
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
        write(output_unit,'(a)') ""
        write(output_unit,'(5x, A)') 'Running P-gTB Density Calculation.'
        Call execute_command_line('gtb coord -stda > gtb.out', WAIT=.true.)
        write(output_unit,'(5x, A)') 'Done.'
        write(output_unit,'(5x, A)') 'Setting up Turbomole control file.'

        nocc=0
        nalpha=0
        nbeta=0
        nopen=0
        open(11,file="gtb.out")

        read(11,'(a)',iostat=io_error) line
        do while (io_error .eq. 0)
            if (index(line,'nalpha') .ne. 0) then
                read(line,*) dummy, nalpha
                read(11,*) dummy, nbeta
                read(11,*) dummy, nocc
                read(11,*) dummy, nopen
                nocc=nocc/2
                exit
            end if
            read(11,'(a)',iostat=io_error) line
        end do

        close(11)

        if (nopen .eq. 0) then
            call rename("mos.tmp", "mos",io_error)
        else
            call rename("alpha.tmp","alpha",io_error)
            call rename("beta.tmp","beta",io_error)
        end if

        if (io_error .ne. 0) then
            call fatal_error(error,"Could not rename MO files.")
            return
        end if

        !> If there is a control file, delete it
        open(11,file="control",status="old",iostat=io_error)
        if (io_error .eq. 0) close(11, status="delete")

        open(12,file="control.tmp",status="old",iostat=io_error)

        if (io_error .ne. 0) then
            Call fatal_error(error, "Could not find control.tmp. This should be created by gtb.")
            return
        end if

        open(11, file="control", status="new")
            write(11,'(a)') &
                "$symmetry c1", &
                "$coord file=coord", &
                "$energy file=energy", &
                "$grad file=gradient", &
                "$atoms"
        read(12,'(a)', iostat=io_error) line

        do while (io_error .eq. 0)
            write(11,'(a)') line
            read(12,'(a)', iostat=io_error) line
        end do

        close(12)

        write(11,'(a)') "$scfiterlimit 2"
        write(11,'(a,F6.2,a)') "$scfdamp start=",damp," step=0.0"
        write(11,'(a)') &
            "$scfdiis maxiter=0", &
            "$dft", &
            " functional libxc 117", &
            " functional libxc add 1 130", &
            " gridsize 1", &
            "$ricore 64000", &
            "$rij", &
            "$scfconv 5"!, &
            !"$disp4 --param 1.0 0.36 0.627 2.836 0.0"
        if (nopen .eq. 0) then
            write(11,'(a)') &
                "$scfmo file=mos", &
                "$closed shells"
            write(11,'(a,I4,a)') " a 1-",nocc,"   (2)"
        else
            write(11,'(a)') &
                "$uhfmo_alpha   file=alpha", &
                "$uhfmo_beta     file=beta", &
                "$uhf", &
                "$alpha shells"
            write(11,'(a,I4,a)') " a 1-",nalpha,"   (1)"
            if (nbeta .ne. 0) then
                write(11,'(a)') "$beta shells"
                write(11,'(a,I4,a)') " a 1-",nbeta,"   (1)"
            end if
        end if

        write(11,'(a)') "$end"
        close(11)

        write(output_unit,'(5x, A)') 'Running gas phase single point on P-gTB density.'
        Call execute_command_line("ridft > gas.out 2>/dev/null", WAIT=.true.)
        write(output_unit,'(5x, A)') 'Done.'

        open(11,file="gas.out")

        read(11,'(a)', iostat=io_error) line

        i=1
        energies(:)=0.0_wp

        do while (io_error .eq. 0)
            if (index(line,"ITERATION") .ne. 0) then
                read(11,*) dummyi, energies(i)
                i=i+1
            end if
            if (i .eq. 3) exit
            read(11,'(a)', iostat=io_error) line
        end do
        
        close(11)
        if (i .ne. 3) then
            Call fatal_error(error, "Could not read two iteration energies from gas.out.")
            RETURN
        end if


        open(11,file="gas.energy")
            write(11,*) energies(2)+((energies(2)-energies(1))*scale)
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
        write(11,'(A19, A4)')'   epsilon=infinity', merge(' ion','    ',ion)
        write(11,'(a)')'   force_cosmowrite'
        write(11,'(a)') '$cosmo_out file=solute.cosmo'
        write(11,'(A4)') '$end'

        write(output_unit,'(5x, A)') 'Running single point on P-gTB density with epsilon=infinity.'
        Call execute_command_line("ridft > solvent.out 2>/dev/null", WAIT=.true.)
        write(output_unit,'(5x, A)') 'Done.'

        open(11,file="solvent.out")

        read(11,'(a)', iostat=io_error) line

        i=1
        energies(:)=0.0_wp

        do while (io_error .eq. 0)
            if (index(line,"ITERATION") .ne. 0) then
                read(11,*) dummyi, energies(i)
                i=i+1
            end if
            if (i .eq. 3) exit
            read(11,'(a)', iostat=io_error) line
        end do
        
        close(11)
        if (i .ne. 3) then
            Call fatal_error(error, "Could not read two iteration energies from solvent.out.")
            RETURN
        end if

        open(11,file="solute.energy")
        write(11,*) energies(2)+((energies(2)-energies(1))*scale)
        close(11)
        write(output_unit,'(a)') ""


    end subroutine gtb
        





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
        character(len=:), allocatable :: solute
        !> Input file existing?
        logical :: ex
        !> Reading Output
        character(len=100) :: line
        character(len=25) :: words, dummy1, dummy2, dummy3, dummy4
        real(wp) :: E_gas, E_solv
        !> I/O error
        integer :: io_error, file_size

        !> Needed for determining right radii
        integer :: atoms, i
        character(len=2), allocatable :: symbols(:)
        real(wp), allocatable :: rad(:)

        !> Default Functional is r2scan-3c
        if (present(new_basis)) then
            functional=new_functional
            basis=new_basis
        else
            functional="r2scan-3c"
            basis=""
        end if

        solute=input(1:len(input)-4)

        INQUIRE(file=input, exist=ex)

        if (.not. ex) then
            Call fatal_error(error,"Input File for Orca Driver Mode not found.")
            return
        end if
        
        open(11,file=input)
        read(11,*) atoms
        allocate(symbols(atoms))
        read(11,*)
        do i=1,atoms
            read(11,*) symbols(i)
        end do
        close(11)

        rad = get_vdw_rad_cosmo(symbols)
        
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
            read(11,*) charge
            close(11)
        end if

        INQUIRE(file=".UHF",exist=ex)
        if (ex) then
            open(11,file=".UHF")
            read(11,*) multi
            !> .UHF is unpaired Electrons (not multiplicity)
            multi=multi+1
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
        write(11,'(A1,1x,A7,1x,I4,1x,I4,1x,A)') "*","xyzfile",charge,multi,input
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
        if (charge .ne. 0) then
            write(11,'(A6)') "! CPCM"
        else
            write(11,'(A7)') "! CPCMC"
        end if
        write(11,'(a)') "%cpcm"
        do i=1,size(symbols)
            write(11,'(3x,a,i0,a,F4.2,a)') "AtomRadii(",i-1,",",rad(i)*BtoA,")"
        end do
        write(11,'(3x,a)') "end"
        write(11,'(A1,1x,A7,1x,I4,1x,I4,1x,A)') "*","xyzfile",charge,multi,input
        close(11)

        if (charge .ne. 0) then
            write(output_unit,'(5x, A)') 'Starting CPCM calculation with epsilon=infinity.'
        else
            write(output_unit,'(5x, A)') 'Starting CPCMC calculation with epsilon=infinity.'
        end if
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
    subroutine turbomole(epsilon, cosmo_out, error, isodens, solvent)
        !> Dielectric Constant for COSMO Calculation
        real(wp), intent(inout) :: epsilon
        !> Output File for COSMO Calculation
        character(len=*), intent(in) :: cosmo_out
        !> Error Handling
        type(error_type), intent(out), allocatable :: error
        !> SMD Solvent for default Epsilon (optional)
        character(len=*), intent(in), optional :: solvent
        !> Isodens Cavity?
        logical, intent(in) :: isodens

        !> Control File exists?
        logical :: ex
        !> Reading the control file
        character(len=100) :: line, line2
        !> Charge and multiplicity of molecule (from control file)
        integer :: charge, multi
        !> Is Molecule charged?
        logical :: ion
        !> Catches some Errors.
        integer :: io_error
        !> Necessary for reading gas phase energy.
        integer :: lines, i
        real(wp) :: E_gas

        charge=0
        multi=0
        if (epsilon .lt. 0) epsilon=minnesota_eps(solvent,error)

        write(output_unit,'(a)') ""
        write(output_unit,'(10x,a)') &
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
            " ------------------------------------------------- ",&
            "|                TURBOMOLE Driver                 |",&
            " ------------------------------------------------- "
            !< < < < < < < < < < < < < > > > > > > > > > > > > >!
        write(output_unit,'(a)') ""
        io_error=0

        Call prepTM('r2scan-3c','def2-mTZVPP',error)
        if (allocated(error)) return

        INQUIRE(file='control', exist=ex)
        if (.not. ex) then
            Call fatal_error(error,'There is no control file in the working directory.&
            &This should not happen.')
            write(output_unit,'(a)') ""
        end if

        open(11, file="control" ,status='old')
        do while (io_error .ge. 0)
            read(11,'(A)',iostat=io_error) line
            if (index(line,"cosmo") .ne. 0) then
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
        read(11,'(a)',iostat=io_error) line
        do while (io_error .eq. 0)
            if (index(line,"abnormal") .ne. 0) then
                Call fatal_error(error,"Gas phase calculation stopped abnormally.&
                &You should maybe try setting up a custom control file.") 
                return
            end if
            read(11,'(a)',iostat=io_error) line
        end do
        close(11)

        write(output_unit,'(5x, A)') &
            "Done! Gas phase calculation terminated normally.", &
            ""

        io_error=0
        open(11, file="energy", status="OLD")
        if (io_error .ne. 0) then
            Call fatal_error(error,'Error while reading the gas phase energy.&
            &Check your control file.')
            return
        end if
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
        if (isodens) then
            write(11,'(a)') &
            " $cosmo_isodens", &
            "   dx=0.1"
        else 
            write(11,'(a)') &
                " cavity closed", &
                " use_contcav", &
                " nspa=272", &
                " nsph=162"
        end if
        write(11,'(A16,A)') '$cosmo_out file=',cosmo_out
        write(11,'(A4)') '$end'
        if (epsilon .ne. 0) then
            write(output_unit,'(5x,A,F0.2,A)') 'Starting COSMO Calculation with epsilon=',epsilon, merge(' ion','    ',ion)
        else
            write(output_unit,'(5x,A,A)') 'Starting COSMO Calculation with epsilon=infinity', merge(' ion','    ',ion)
        end if
        Call execute_command_line('ridft > solv.out 2>error', WAIT=.true.)
        open(11,file="error")
        read(11,'(a)',iostat=io_error) line
        do while (io_error .eq. 0)
            if (index(line,"abnormal") .ne. 0) then
                if (isodens) then
                    open(12,file="solv.out")
                    read(12,'(a)',iostat=io_error) line
                    do while (io_error .eq. 0) 
                        if (index(line,"min. dist. between &
                        &inner and outer cavity smaller than 0.5*rsolv") .ne. 0) then
                            write(output_unit,'(5x,a, t20, a)') &
                            "","",&
                            "[WARNING]", "Failure in outer charge correction while using isodens cavity.", &
                            "","Trying again with larger routf value (routf = 1.1).", &
                            "","Carefully check your results.", &
                            "",""
                            Call execute_command_line("kdg end")
                            Call execute_command_line("kdg cosmo")
                            open(13,file="control", access="append")
                            write(13,'(a)'), &
                            "$cosmo", &
                            "   epsilon=infinity", &
                            "   routf=1.1", &
                            "$cosmo_isodens", &
                            "   dx=0.1", &
                            "$cosmo_out file=solute.cosmo",&
                            "$end"
                            close(13)
                            Call execute_command_line('ridft > solv.out 2>error', WAIT=.true.)
                            open(13, file="error")
                            read(13,'(a)',iostat=io_error) line2
                            do while (io_error .eq. 0)
                                if (index(line2,"abnormal") .ne. 0) then
                                    Call fatal_error(error,"Failure in outer charge corrections while using isodens cavity.")
                                    return
                                end if
                                read(13,'(a)',iostat=io_error) line2
                            end do
                            exit
                        end if
                        read(12,'(a)',iostat=io_error) line
                    end do
                    close(12)
                else
                    Call fatal_error(error,"Solvent phase calculation stopped abnormally.") 
                    return
                end if
                exit
            end if
            read(11,'(a)',iostat=io_error) line
        end do
        close(11)
        write(output_unit,'(5x, A)') &
            "Done! COSMO phase calculation terminated normally.", &
            ""

    end subroutine turbomole

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
            write(11,'(3x,I0,3x,F17.10,3x,F17.10,3x,F17.10,3x,a)') i, xyz(i,1), xyz(i,2), xyz(i,3), symbols(i)
        end do
        write(11,'(a)') "$cosmo_energy"
        write(11,'(a)') "#CPCM Singlepoint Energy taken from Orca"
        write(11,'(a,F15.9)') "Final CPCM Single Point energy [a.u.] = ",energy
        write(11,'(a)') ""
        write(11,'(a)') "$segment_information"
        write(11,'(a)') "#Reordered segment information from "//oc_inp//".cpcm"
        do i=1,segments
            write(11,('(3x,I0,3x,I3,3x,F15.9,3x,F15.9,3x,&
            &F15.9,3x,F15.9,3x,F15.9,3x,F15.9,3x,F15.9)'))&
            & int(segment_info(i,1)), int(segment_info(i,2)), segment_info(i,3),&
            & segment_info(i,4), segment_info(i,5), segment_info(i,6), segment_info(i,7),&
            & segment_info(i,8), segment_info(i,9)
        end do

    end subroutine orcatocosmo

    subroutine prepTM(functional,basis, error)
        use sort, only: unique
        !> Functional and Basis for control file
        character(len=*), intent(in) :: functional, basis
        !> Error Handling
        type(error_type), allocatable, intent(out) :: error

        !> Elements of the molecule
        character(len=2), allocatable :: elements(:)
        !> Unique elements
        character(len=2), allocatable :: unique_elements(:)

        !> File I/O Handling
        character(len=150) :: line
        integer :: ele_count
        integer :: io_error
        logical :: ex
        
        !> Charge and Multiplicity of the System
        integer :: charge, multi

        !> Dummy Variables and Loop Variable
        integer :: i, j
        real(wp) :: d1, d2, d3

        character(len=2), parameter :: ecp28(32) = ['rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cu', &
                                    &'in','sn','sb','te','i ','xe','ce','pr','nd','pm','sm','eu',&
                                    &'gd','tb','dy','ho','er','tm','yb','lu'], &
                                &ecp46(3)=['cs','ba','la'],&
                                &ecp60(15)=['hf', 'ta', 'w ', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi',&
                                    &'po', 'at', 'rn']

                        

        
        Inquire(file='coord', exist=ex)

        if (.not. ex) then
            Call fatal_error(error, 'There is no coord file in the working directory.')
            return
        end if
        
        Inquire(file='control', exist=ex)
        if (ex) then
            write(output_unit,'(5x,a, t20, a)') &
                "[WARNING]", "Found a control file in the working directory.", &
                "","CPCM-X will use this control file instead of writing its own.", &
                "","Carefully check your results."
            write(output_unit,'(a)') ""
            return
        end if

        write(output_unit,'(a)') ""
        write(output_unit,'(5x,a)') "Using internal prepTM routine to write a control file."
        write(output_unit,'(5x,a,a,a,a)') "Using functional: ", functional, " and basis: ", basis

        open(11,file='coord')
        read(11,*,iostat=io_error) line
        ele_count=0

        do while (io_error .eq. 0)
            read(11,'(a)',iostat=io_error) line
            if (index(line,'$') .ne. 0) exit
            ele_count=ele_count+1
        end do

        rewind(11)
        
        allocate(elements(ele_count))
        read(11,('(a)')) line
        do i=1,ele_count
            read(11,*) d1, d2, d3, elements(i)
            elements(i)=to_lower(elements(i))
        end do
        close(11)

        charge=0
        multi=0

        INQUIRE(file=".CHRG",exist=ex)
        if (ex) then
            open(11,file=".CHRG")
            read(11,*) charge
            close(11)
        end if

        INQUIRE(file=".UHF",exist=ex)
        if (ex) then
            open(11,file=".UHF")
            read(11,*) multi
            close(11)
        end if
    
        unique_elements=unique(elements)

        open(11,file="control", status='new', iostat=io_error)
        if (io_error .ne. 0) then
            Call fatal_error(error,'Something went wrong while writing the control file.')
            return
        end if

        write(11,'(a)') &
        "$symmetry c1", &
        "$coord file=coord", &
        "$energy file=energy", &
        "$atoms"
        do i=1,size(unique_elements)
            write(11,'(a,1x, a,1x, a)') unique_elements(i), positions(unique_elements(i),elements), '\'
            write(11,'(a,a,1x,a)') '  basis =',unique_elements(i), basis
            if (any(ecp28 .eq. unique_elements(i))) then
                write(11,'(a,a,a)') '  jbas =',unique_elements(i),' universal-ecp-28'
            else if (any(ecp46 .eq. unique_elements(i))) then
                write(11,'(a,a,a)') '  jbas =',unique_elements(i),' universal-ecp-46'
            else if (any(ecp60 .eq. unique_elements(i))) then
                write(11,'(a,a,a)') '  jbas =',unique_elements(i),' universal-ecp-28'
            else
                write(11,'(a,a,a)') '  jbas =',unique_elements(i),' universal'
            end if
        end do
        write(11,'(a)') &
        "$rij", &
        "$dft"
        write(11,'(a,a)') " functional ",functional
        write(11,'(a)') &
        " gridsize m4", &
        " radsize 8", &
        "$scfconv 6", &
        "$escfiterlimit 250", &
        "$ricore    16000"
        if (multi .ne. 0) then
            write(11,'(a,I0,a,I0)') "$eht charge=",charge," unpaired=", multi
        else
            write(11,'(a,I0)') "$eht charge=",charge
        end if
        write(11,'(a)') "$end"
        write(output_unit,'(5x,a)') "Done."
        
        write(output_unit,'(a)') ""

    end subroutine prepTM

    function positions(unique_element,element_array) result(short_pos)
        !> Element Array for which the positions should be determined
        character(len=*), allocatable, intent(in) :: element_array(:)
        !> Unique Element in Element Array
        character(len=*), intent(in) :: unique_element
        !> Short Positional character string for control file (e.g. 1-4,13)
        character(:), allocatable :: short_pos

        integer, allocatable :: position(:)

        integer, allocatable :: tmp_pos(:)
        
        character(len=3) :: tmp

        integer :: i, pos_num, j, break

        pos_num=0

        do i=1,size(element_array)
            if (unique_element .eq. element_array(i)) then
                pos_num=pos_num+1
                if (allocated(position)) deallocate(position)
                allocate(position(pos_num))
                do j=1,pos_num-1
                    position(j)=tmp_pos(j)
                end do
                position(pos_num)=i
                if (allocated(tmp_pos)) deallocate(tmp_pos)
                allocate(tmp_pos(pos_num))
                tmp_pos=position
            end if
        end do

        break=1
        short_pos=''
        write(tmp,'(I0)') position(1)
        short_pos=trim(tmp)
        do i=2,size(position)
            if (position(i) .eq. position(i-1)+1) then
                if (i .eq. size(position)) then
                    write(tmp,'(i0)') position(i)
                    short_pos=short_pos//'-'//trim(tmp)
                    exit
                end if
                if (position(i) .ne. position(i+1)-1) then
                    write(tmp,'(i0)') position(i)
                    short_pos=short_pos//'-'//trim(tmp)
                end if
            cycle
            end if
            short_pos=short_pos//','
            if (len(short_pos)/break .gt. 50) then
                short_pos=short_pos//'\'//new_line('a')
                break=break+1
            end if
            write(tmp,'(i0)') position(i)
            short_pos=short_pos//trim(tmp)
        end do
    end function positions


end module qc_calc  
