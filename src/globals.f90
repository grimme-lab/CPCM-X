module globals
   use element_dict
   implicit none
   
   ! Global Parameters
   real(8), parameter :: autokcal=627.509_8, R=8.314_8, jtokcal=0.000239006
   real(8), parameter :: k_b=1.380649E-23,N_A=6.02214076E23, BtoA=0.529177
   real(8), parameter :: pi = 4*atan(1.0_8)
   real(8), dimension(10) :: param, pr_param !Parameter for several methods
   real(8) :: SysTemp !System Temperature

   real(8) :: dG_is, dG_cc, dG_res, dG_disp !contributions to free energy of solvation
   character(7) :: model ! Chosen COSMO model
   logical :: onlyprof ! Only creates a sigma profile
   logical :: ML !Data Preperation for ML

   ! Covalent Radii of Elements

   type(DICT_STRUCT), pointer :: cov_r

   ! PR Parameters

   type(DICT_STRUCT), pointer :: A_dsp, B_dsp, r_pr


   contains 
      ! Functions that are used in several modules
      function distance(xyz1,xyz2)
         implicit none
         real(8), intent(in) :: xyz1(:), xyz2(:)
         real(8) :: distance

         distance=(sqrt((xyz1(1)-xyz2(1))**2.0_8+(xyz1(2)&
                  &-xyz2(2))**2.0_8+(xyz1(3)-xyz2(3))**2.0_8))
         
      end function distance

   Pure Function to_upper (str) Result (string)

   !   ==============================
   !   Changes a string to upper case
   !   ==============================

      Implicit None
      Character(*), Intent(In) :: str
      Character(LEN(str))      :: string

      Integer :: ic, i

      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

   !   Capitalize each letter if it is lowecase
       string = str
       do i = 1, LEN_TRIM(str)
         ic = INDEX(low, str(i:i))
           if (ic > 0) string(i:i) = cap(ic:ic)
       end do

   End Function to_upper

   Pure Function to_lower (str) Result (string)

   !   ==============================
   !   Changes a string to upper case
   !   ==============================

      Implicit None
      Character(*), Intent(In) :: str
      Character(LEN(str))      :: string

      Integer :: ic, i

      Character(26), Parameter :: low = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       Character(26), Parameter :: cap = 'abcdefghijklmnopqrstuvwxyz'

   !   Capitalize each letter if it is lowecase
       string = str
       do i = 1, LEN_TRIM(str)
         ic = INDEX(low, str(i:i))
           if (ic > 0) string(i:i) = cap(ic:ic)
       end do

   End Function to_lower
end module globals
