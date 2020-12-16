module initialize_cosmo

   contains
      subroutine read_cosmo(compound,ident,xyz,charges,area,pot,volume)
         use globals
         implicit none
         character(*), intent(in) :: compound
         character(30) :: line 
         real(8), dimension(:), allocatable,intent(out) :: charges,&
            &area,pot
         real(8), dimension(:,:), allocatable, intent(out) :: xyz
         character(2), allocatable, dimension(:),intent(out) :: ident
         real(8), intent(out) :: volume
         integer :: i, io_error, dummy1, dummy2, num
         real(8) :: dummy3, dummy4, dummy5
         real(8), dimension(:), allocatable :: dummy_ident
         character(2) :: element
         
         logical :: exists

         INQUIRE(file=compound,EXIST=exists)
   
         if (.NOT. exists) then
            write(*,*) "COSMO file does not exists."
            stop
         end if

         open(1,file=compound,iostat=io_error)

         if (io_error .NE. 0) then
            write(*,*) "Problem while reading COSMO input. &
               &Check your data."
            stop
         end if

         io_error=0 
         do while (io_error .GE. 0)
            read(1,*,iostat=io_error) line
         end do
         read(line,*) num
         allocate(charges(num))
         allocate(ident(num))
         allocate(dummy_ident(num))
         allocate(area(num))
         allocate(xyz(num,3))
         allocate(pot(num))
         rewind(1)
         ident(:)="XX"
         do while (line .NE. "$segment_information")
            read(1,*) line
         end do
         
         num=1
         dummy4=0
         do while (.TRUE.)
            read(1,'(A1)',advance='no',iostat=io_error) line
            if (line=="#") then
               read(1,*)
               cycle
            else if (io_error .LT. 0) then
               exit
            else
               read(1,*) dummy1,dummy_ident(num),xyz(num,1),xyz(num,2),&
                  &xyz(num,3),dummy3,area(num),charges(num),pot(num)
            !   charges(num)=anint(charges(num)*1000)/1000
               pot(num)=pot(num)*BtoA
               num=num+1
            end if
         end do
         rewind(1)
         do while (line .NE. "#atom")
            read(1,*) line
         end do

         do while (.TRUE.)
            read(1,'(A1)',advance='no',iostat=io_error) line
            if (line .NE. "$") then
               read(1,*) dummy1, dummy3, dummy4, dummy5, element
               do i=1,size(dummy_ident)
                  if (dummy1==dummy_ident(i)) then
                     ident(i)=element
                  end if
               end do
            else
               exit
            end if
         end do

         rewind(1)
         do while (line .NE. "area=")
            read(1,*) line
         end do

         read(1,*) line, volume
         
         close(1)
         deallocate(dummy_ident)
      end subroutine read_cosmo

      subroutine initialize_param(param,r_cav,disp_con)
         use element_dict
         implicit none

         real(8), dimension(10) :: param
         type(DICT_STRUCT), pointer, intent(inout) :: r_cav, disp_con

         type(DICT_DATA) :: data1
         character(len=2) :: symbol
         
         !Setting global COSMO Parameters - TO DO: read from file

         param(1)=0.40_8
         param(2)=0.8_8
         param(3)=1515.0_8
         param(4)=2.802_8
         param(5)=7400.0_8
         param(6)=0.00854_8
         param(7)=7.62_8
         param(8)=0.129_8
         param(9)=-0.217_8
         param(10)=-9.910_8

         !Setting cavity radii - to do: read from file

         data1%param = 1.30_8
         call dict_create(r_cav, 'h', data1)
         data1%param = 2.00_8
         call dict_add_key(r_cav, 'c', data1)
         data1%param = 1.72_8
         call dict_add_key(r_cav, 'o', data1)

         !Setting dispersion contant -- need to be read from file

         data1%param = -0.041_8
         call dict_create(disp_con, 'h', data1)
         data1%param = -0.037_8
         call dict_add_key(disp_con, 'c', data1)
         data1%param = -0.042_8
         call dict_add_key(disp_con, 'o', data1)



      end subroutine initialize_param
      subroutine getargs(solvent,solute,T)
      
         integer :: i
         character(len=20),intent(out) ::solvent, solute
         real(8),intent(out) :: T

         character(len=20) :: arg

         T=298.15_8
         solute=''
         solvent=''

         do i=1,command_argument_count()

         Call get_command_argument(i, arg) 
            Select case (arg)
               case ("--c")
                  Call get_command_argument(i+1,arg)
                  solute=arg
               case ("--s")
                  Call get_command_argument(i+1,arg)
                  solvent=arg
               case ("--T")
                  Call get_command_argument(i+1,arg)
                  read(arg,*) T
               case default
                  cycle
            end select
         end do
      end subroutine getargs
end module initialize_cosmo
