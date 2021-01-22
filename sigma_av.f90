module sigma_av

   ! Module contains charge averaging and orthogonalizing algorythms
   !
   ! Averaging Uses:
   ! r_av : averaging radius
   ! xyz : Array of the coordinates of the segments
   ! charges : Array of Charge/Area of the Segments
   ! Area : Array of the Area of the Sements
   ! Results:
   ! av_charge : Averaged charge/area over the averaging radius r_av
   !
   ! Orthogonalizing Uses:
   ! v = Averaged charge/area over the averaging radius r_av
   ! v0 = Averaged charge/area over the averaging radius 2*r_av
   ! Results:
   ! vt = Orthogonalized charge/area
   


   contains

      subroutine average_charge (r_av, xyz, charges, area,av_charge)
         use globals
         implicit none
         real(8), dimension(:), allocatable, intent(in) :: charges, area
         real(8), dimension(:,:), allocatable, intent(in) :: xyz
         real(8), intent(in) :: r_av
         real(8), dimension(:), allocatable, intent(out) :: av_charge
         real(8) :: tmpcharge, tmpcounter, tmpdenominator, r_u2, r_av2
         integer :: num, i, j

         num = size(charges)
         allocate(av_charge(num))
         r_av2=r_av**2.0_8
         r_u2=0.0_8
         tmpcharge=0.0_8
         tmpcounter=0.0_8
         tmpdenominator=0.0_8
         do i=1,num
            do j=1,num
               r_u2=(area(j)/pi)
               tmpcounter=tmpcounter+(charges(j)*((r_u2*r_av2)/(r_u2+r_av2))*&
                         &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_8)/(r_u2+r_av2))))
               tmpdenominator=tmpdenominator+(((r_u2*r_av2)/(r_u2+r_av2))*&
                             &exp(-((distance(xyz(j,:),xyz(i,:))**2.0_8)/(r_u2+r_av2))))
               r_u2=0.0_8
            end do
            tmpcharge=tmpcounter/tmpdenominator
            tmpcounter=0.0_8
            tmpdenominator=0.0_8
            av_charge(i)=tmpcharge
         end do
      end subroutine average_charge


      subroutine ortho_charge (v,v0,vt)
         implicit none
         real(8), dimension(:), allocatable, intent(in) :: v,v0
         real(8), dimension(:), allocatable, intent(out) :: vt

         integer :: i

         allocate(vt(size(v)))

         do i=1,size(v)
            vt(i)=v0(i)-0.816_8*v(i)
         end do

      end subroutine
end module sigma_av
