module disp

   !Peng-Robinson Equation of State model for disperive interactions

   contains


      
      subroutine pr2018(area,elements,ident,oh_count,nh_count,n_ear)
         use globals
         use element_dict
         implicit none
         !New version of PR EOS 
         real(8), dimension(:), allocatable, intent(in) :: area, ident
         character(2), dimension(:), allocatable, intent(in) :: elements
         integer, intent(in) :: oh_count, nh_count, n_ear

         character(2) :: element
         real(8), dimension(:), allocatable :: atom_area
         real(8) :: dG_hb, dG_ring, dG_vdw
         integer :: i
         type(DICT_DATA) :: A, B
        
         !write(*,*) n_ear
         dG_hb = 0.0_8
         dG_ring = 0.0_8
         dG_vdw = 0.0_8
         allocate(atom_area(int(maxval(ident))))
         atom_area=0.0_8
         do i=1,size(area)
            atom_area(int(ident(i)))=atom_area(int(ident(i)))+area(i)
         end do

        ! do i=1,int(maxval(ident))
        !    write(*,*) i, atom_area(i)
        ! end do

        do i=1,int(maxval(ident))
         element=elements(i)
         A=dict_get_key(A_dsp,element)
         B=dict_get_key(B_dsp,element)
        ! write(*,*) atom_area(i), A%param, B%param
         dG_vdw=dG_vdw+atom_area(i)*(A%param*log(SysTemp)+B%param)
        end do
       !  write(*,*) oh_count
         dG_hb=oh_count*(pr_param(3)*(log(SysTemp)/SysTemp)+pr_param(4))&
            &+nh_count*(pr_param(5)*(log(SysTemp)/SysTemp)+pr_param(6))
!write(*,*) dG_hb
         dG_ring=n_ear*(pr_param(1)*log(SysTemp)+pr_param(2))

         dG_vdw=dG_vdw*jtokcal
         
         dG_disp = dG_vdw + dG_hb + dG_ring
        ! write(*,*)  dG_vdw, dG_hb, dG_ring, dG_disp
        if (ML) then
           open(5,file='ML.pr')
           write(5,'(F0.2,A,I0,A,I0,A,I0,A,I0)') SysTemp,",",oh_count,",",nh_count,",",n_ear,",",int(maxval(ident))
           close(5)
         end if


      end subroutine pr2018

end module disp
