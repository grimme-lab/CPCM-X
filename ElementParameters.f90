module EleData_Module
   implicit none
   integer, parameter :: DICT_KEY_LENGTH = 2

   type EleData
      real(8) :: param
   end type EleData

   type(EleData), parameter :: DICT_NULL = eledata(0)
end module

module Element_Dict
   use EleData_Module, DICT_DATA => EleData
   implicit none

   include "dictionary.f90"
end module Element_Dict
