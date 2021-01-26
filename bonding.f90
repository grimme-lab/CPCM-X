module bonding
   
   ! This module is used to determine the Bonding Situation of the Compounds

   contains

   ! This Subroutine determines covalent bonds between Atoms/Segments 
   ! It uses the Positions of the Segments and the hard coded covalent Radii
   ! A Covalent Bond is assumed when the distance is lower than the sum of the convalent Radii*1.15
   subroutine det_bonds(ident,xyz,elements,is_bonded)
      use globals
      use element_dict
      implicit none

      real(8), dimension(:), allocatable, intent(in) :: ident
      character(2), dimension(:), allocatable, intent(in) :: elements
      real(8), dimension(:,:), allocatable, intent(in) :: xyz

      logical, dimension(:,:), allocatable, intent(out) :: is_bonded

      type(DICT_DATA) :: radius
      integer :: i,j
      real(8) :: cov_a, cov_b

      allocate(is_bonded(size(ident),size(ident)))
      is_bonded(:,:) = .FALSE.

      cov_a=0.0_8
      cov_b=0.0_8
      do i=1,int(maxval(ident))
         radius=dict_get_key(cov_r,elements(i))
         cov_a=radius%param
         do j=i+1,int(maxval(ident))
            radius=dict_get_key(cov_r,elements(i))
            cov_b=radius%param
            if (distance(xyz(i,:),xyz(j,:)) .LE. ((cov_a+cov_b)*1.15)) then
               is_bonded(i,j)=.TRUE.
               is_bonded(j,i)=.TRUE.
            end if
            cov_b=0.0_8
         end do
         cov_a=0.0_8
      end do


      ! Check points, uncomment the following:
     ! do i=1,int(maxval(ident))
     !    write(*,*) i, elements(i)
     !    do j=1,int(maxval(ident))
     !       write(*,*) "to ",j,is_bonded(i,j)
     !    end do
     ! end do
         
   end subroutine

   subroutine hb_grouping(ident,elements,is_bonded,hb_group)
      implicit none
      
      real(8), dimension(:), allocatable, intent(in) :: ident
      character(2), dimension(:), allocatable, intent(in) :: elements
      logical, dimension(:,:), allocatable, intent(in) :: is_bonded

      character(2), dimension(:), allocatable, intent(out) :: hb_group

      
      integer :: i,j

      
      allocate(hb_group(size(ident)))

      hb_group(:) = 'NH'

      do i=1,size(ident)
         if (elements(int(ident(i))) .EQ. 'h') then
            do j=1,size(ident)
               if (is_bonded(int(ident(i)),int(ident(j)))) then
                  if (elements(int(ident(j))) .EQ. 'o') then
                     hb_group(i)='OH'
                     hb_group(j)='OH'
                     exit
                  else if ((elements(int(ident(j))) .EQ. 'n') .OR. &
                     &(elements(int(ident(j))) .EQ. 'f')) then
                     hb_group(i)='OT'
                     hb_group(j)='OT'
                     exit
                  end if
               end if
            end do
         end if

         if (elements(int(ident(i))) .EQ. 'o') then
            do j=1,size(ident)
               if (is_bonded(int(ident(i)),int(ident(j)))) then
                  if (elements(int(ident(j))) .EQ. 'h') then
                     hb_group(i)='OH'
                     hb_group(j)='OH'
                     exit
                  end if
               end if
            end do
            if (hb_group(i) .EQ. 'NH') then
               hb_group(i) = 'OT'
            end if
         end if 

         if (elements(int(ident(i))) .EQ. 'n') hb_group(i)='OT' 
         if (elements(int(ident(i))) .EQ. 'f') hb_group(i)='OT'

       ! write(*,*) i, elements(int(ident(i))), hb_group(i)
      end do

!      write(*,*) hb_group(66),size(ident)
!      do i=1,size(ident)
!         write(*,*) i,elements(int(ident(i))), is_bonded(i,66), is_bonded(66,i)
!      end do

   end subroutine hb_grouping
      

      
end module bonding
