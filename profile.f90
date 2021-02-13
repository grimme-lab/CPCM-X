module profile

   !! This Module contains algorythms for reading sigma profiles from file
   !! or creating Sigma Profiles from COSMO Files.

   contains

      subroutine read_singlesig(sigma,nam,volume)

         !! Read a Single Sigma Profile from File, uses the format of the UD Database
         !! Input:
         !! nam: path to the sigma file
         !! Output:
         !! Sigma: Sigma Profile array ordered from -0.25 to +0.25
         !! Volume: Volume of the Compound from COSMO Calculation

         character(len=*), intent(in) :: nam
         real(8), intent(out) :: volume
         real(8), dimension(:), intent(out) :: sigma (0:50)

         character(len=20) :: dummy1, dummy2
         integer :: io_error,i
         real(8) :: dummy3

         open(unit=2,file=nam)

         read(2,*)
         read(2,*)
         read(2,*) 
         read(2,*) dummy1,dummy2,volume
     	 io_error=0
         sigma(:) = 0.0_8
         i=0
     	 do i=0,50
     	    read(2,*,iostat=io_error) dummy3, sigma(i)
     	 end do
         close(2)
     	
      end subroutine read_singlesig

      subroutine read_triplesig(sigma3,nam,volume)

         !! Read Splitted Sigma Profile from File, uses the format of the UD Database
         !! Input:
         !! nam: path to the sigma file
         !! Output:
         !! Sigma3: Splitted Sigma Profile array ordered from -0.25 to +0.25
         !! Volume: Volume of the Compound from COSMO Calculation

         character(len=*), intent(in) :: nam
         real(8), intent(out) :: volume
         real(8), dimension(:), intent(out) :: sigma3 (3,0:50)

         character(len=20) :: dummy1, dummy2
         integer :: io_error,i,t
         real(8) :: dummy3

         open(unit=2,file=nam)

         read(2,*)
         read(2,*)
         read(2,*) 
         read(2,*) dummy1,dummy2,volume
     	 io_error=0
         sigma3 = 0.0_8
         i=0
         do t=1,3
     	    do i=0,50
     	       read(2,*,iostat=io_error) dummy3, sigma3(t,i)
     	    end do
         end do
         close(2)
     	
      end subroutine read_triplesig


      subroutine single_sigma(sv,area,sigma,nam)

         !! Creates a single Sigma Profile from averaged charge densities.
         !! Input:
         !! sv: charge densities of the segmens (array)
         !! area: area of the segments (array)
         !! nam: Path where the sigma file should be written
         !! Output:
         !! sigma: Single Sigma Profile for the Compound


         real(8), dimension(:), allocatable,intent(in) :: sv,area

         real(8), dimension(:),intent(out) :: sigma(0:50)

         character(len=*), intent(in), optional :: nam

         integer :: sigma_min, sigma_max, i, j,tmp

         real(4) :: punit

         real(4), parameter :: sig_width=0.025_8
         integer, parameter :: n_sig=51
         real(4) :: counter(0:n_sig-1)

         real(8) :: profile(0:n_sig-1), chdval(0:n_sig-1), temp

         punit=0.001


         profile(:)=0.0_8
         counter(:)=0.0
         
         do i=0,n_sig-1
            profile(i) = 0.0_8
            counter(i) = -sig_width+punit*i
         end do
         do i= 1, size(sv)
            temp = sv(i)
           
            tmp = int((temp-counter(0))/punit)
 
            if (tmp<0) tmp=0
            if (tmp>n_sig-1) tmp=n_sig-1
            profile(tmp) = profile(tmp)+area(i)*(counter(tmp+1)-temp)/punit
            profile(tmp+1) = profile(tmp+1)+area(i)*(temp-counter(tmp))/punit
         end do
         if (present(nam)) then
            open(unit=2,file=nam//"_sigma.txt",action="write",status="replace")
         
            do i=0,size(profile)-1
               write(2,*) counter(i),";", profile(i)!/sum(area)
            end do
            close(2)
         end if

         
         sigma(:)=profile(:)

      end subroutine

      subroutine split_sigma(sv,area,hb_group,ident,elements,sigma3,nam)
      
         !! This Routine splits the charge densities by the Hydrogen Bonding groups
         !! and creates a seperate Sigma Profile for each independent group.
         !! Input:
         !! sv: charge densities of the segments
         !! area: area of the segments
         !! nam: path where the triple sigma profile should be written
         !! hb_group: Hydrogen Bonding Group of the Segments (OH, OT, NH)
         !! Output:
         !! sigma3: Sigma Profiles for each HB Group (1: NH, 2: OH, 3: OT)

         real(8), dimension(:), allocatable,intent(in) :: sv,area
         character(2), dimension(:), allocatable, intent(in) :: hb_group
         character(len=*), intent(in), optional :: nam
         real(8), dimension(:), allocatable, intent(in) :: ident
         character(2), dimension(:), allocatable, intent(in) :: elements

         real(8), dimension(0:50) :: prob_hb !hb bond probability function
         real(8), dimension(3,0:50),intent(out) :: sigma3

         character(2), dimension(:), allocatable :: profile_group
         real(8), dimension(:), allocatable :: sv_oh, sv_ot, sv_nh, area_oh, area_ot, area_nh
         real(8) :: s_0, punit, max_sig,save1,save2,save3
         integer :: oh_count, ot_count, nh_count, i

         allocate(profile_group(size(hb_group)))
         profile_group="NH"
         s_0=0.007 !sigma_0 for the probability based function
         punit=0.001
         max_sig=0.025
         prob_hb=0.0

         ! Calculate the probability for each sigma to form hydrogen bonds

         do i=0,50
            prob_hb(i)=1.0_8-dexp(-((punit*i-max_sig)**2.0_8/(2.0_8*s_0**2.0_8)))
         end do

         oh_count=0
         ot_count=0
         nh_count=0

         ! Choose to which profile (OH,OT,NH) each Segment belongs.

         do i=1,size(sv)
            select case (elements(int(ident(i))))
               case ("o", "n", "f")
                  select case (hb_group(i))
                     case ("OH")

                        if (sv(i) .GT. 0) then
                           oh_count=oh_count+1
                           profile_group(i)="OH"
                        else
                           nh_count=nh_count+1
                           profile_group(i)="NH"
                        end if

                     case ("OT")

                        if (sv(i) .GT. 0) then
                           ot_count=ot_count+1
                           profile_group(i)="OT"
                        else
                           nh_count=nh_count+1
                           profile_group(i)="NH"
                        end if

                     case default

                        nh_count=nh_count+1
                        profile_group(i)="NH"
                  end select

               case ("h")
                  select case (hb_group(i))
                     case ("OH")

                        if (sv(i) .LT. 0) then
                           oh_count=oh_count+1
                           profile_group(i)="OH"
                        else
                           nh_count=nh_count+1
                           profile_group(i)="NH"
                        end if

                     case ("OT")

                        if (sv(i) .LT. 0) then
                           ot_count=ot_count+1
                           profile_group(i)="OT"
                        else
                           nh_count=nh_count+1
                           profile_group(i)="NH"
                        end if

                     case default
                        
                        nh_count=nh_count+1
                        profile_group(i)="NH"
                  end select

               case default

                  nh_count=nh_count+1
                  profile_group(i)="NH"

            end select
        !    write(*,*) hb_group(i)
        !    write(*,*) profile_group(i)
         end do

         ! Allocate the three profile array according to the 
         ! number of segments in each profile
         
         allocate(sv_oh(oh_count))
         allocate(area_oh(oh_count))
         allocate(sv_ot(ot_count))
         allocate(area_ot(ot_count))
         allocate(sv_nh(nh_count))
         allocate(area_nh(nh_count))

         oh_count=0
         ot_count=0
         nh_count=0

         ! Sort the Segments into the accordings profiles
         
         do i=1,size(sv)
            select case (profile_group(i))
               case ("OH")
                  oh_count=oh_count+1
                  sv_oh(oh_count)=sv(i)
                  area_oh(oh_count)=area(i)
               case ("OT")
                  ot_count=ot_count+1
                  sv_ot(ot_count)=sv(i)
                  area_ot(ot_count)=area(i)
               case ("NH")
                  nh_count=nh_count+1
                  sv_nh(nh_count)=sv(i)
                  area_nh(nh_count)=area(i)
               case default
                  cycle
            end select
         end do

         ! Create Split Sigma Profiles for each Profile Group
         
         Call single_sigma(sv_nh,area_nh,sigma3(1,:))
         Call single_sigma(sv_oh,area_oh,sigma3(2,:))
         Call single_sigma(sv_ot,area_ot,sigma3(3,:))
        
         ! Scale Profiles with probability to form hydrogen bond

         do i=0,50
            save1=sigma3(1,i)+((sigma3(2,i)+sigma3(3,i))*(1-prob_hb(i)))
            save2=sigma3(2,i)*prob_hb(i)
            save3=sigma3(3,i)*prob_hb(i)
            sigma3(1,i)=save1
            sigma3(2,i)=save2
            sigma3(3,i)=save3
         end do
         ! If choosen, write Profiles into file

         if (present(nam)) then

            open(unit=2,file=nam//"_sigma3.txt",action="write",status="replace")
         
            do i=0,50
                
               write(2,999) ((i*punit)-max_sig), sigma3(1,i)
            end do

            do i=0,50
               write(2,999) ((i*punit)-max_sig), sigma3(2,i)
            end do

            do i=0,50
               write(2,999) ((i*punit)-max_sig), sigma3(3,i)
            end do

            ! do i=0,50 !check if sigma+sigma+sigma=single_sigma
           !    write(2,*) ((i*0.001)-0.25),";",sigma3(1,i)+sigma3(2,i)+sigma3(3,i)
           ! end do
            close(2)
         end if

         deallocate(sv_oh,area_oh,sv_ot,area_ot,sv_nh,area_nh)
         999 FORMAT(F6.3,5x,F10.6)
      end subroutine split_sigma

end module profile
