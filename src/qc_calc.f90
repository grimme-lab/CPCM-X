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
    use data, only: minnesota_eps
    use globals, only: rename, autokcal
    use, intrinsic :: iso_fortran_env, only : output_unit, file_storage_size

    implicit none
    private
    public :: qc_cal
    interface qc_cal
        module procedure :: turbomole
        module procedure :: orca
        module procedure :: gtb
    end interface qc_cal

contains

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