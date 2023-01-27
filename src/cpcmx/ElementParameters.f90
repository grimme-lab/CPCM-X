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

module EleData_Module
   use mctc_env, only : wp
   implicit none
   integer, parameter :: DICT_KEY_LENGTH = 2

   type EleData
      real(wp) :: param
   end type EleData

   type(EleData), parameter :: DICT_NULL = eledata(0)
end module

module Element_Dict
   use EleData_Module, DICT_DATA => EleData
   implicit none

type LIST_DATA
    character(len=DICT_KEY_LENGTH) :: key
    type(DICT_DATA)                :: value
end type LIST_DATA

type HASH_LIST
    type(LINKED_LIST), pointer :: list
end type HASH_LIST

type DICT_STRUCT
    private
    type(HASH_LIST), pointer, dimension(:) :: table
end type DICT_STRUCT

!
! We do not want everything to be public
!
private :: LIST_DATA
private :: HASH_LIST
private :: LINKED_LIST
private :: list_create
private :: list_destroy
private :: list_count
private :: list_next
private :: list_insert
private :: list_insert_head
private :: list_delete_element
private :: list_get_data
private :: list_put_data
private :: dict_get_elem
private :: dict_hashkey

integer, parameter, private :: hash_size  = 4993
integer, parameter, private :: multiplier = 31

type LINKED_LIST
    type(LINKED_LIST), pointer :: next
    type(LIST_DATA)            :: data
end type LINKED_LIST

!
! define a private (!) interface to prevent
! mistakes with ordinary assignment
!
!interface assignment(=)
!    module procedure list_assign
!end interface
!private :: list_assign

!
! Define the subroutines and functions
!
contains

! list_assign
!     Subroutine to prevent errors with assignment
! Arguments:
!     list_left   List on the left-hand side
!     list_right  List on the right-hand side
!
! NOTE:
!     This does not work because of a private/public
!     conflict
!
!subroutine list_assign( list_left, list_right )
!    type(LINKED_LIST), INTENT(OUT)  :: list_left
!    type(LINKED_LIST), INTENT(IN)   :: list_right
!   !type(LINKED_LIST), pointer      :: list_left
!   !type(LINKED_LIST), pointer      :: list_right
!
!    !
!    ! Note the order!
!    !
!    stop 'Error: ordinary assignment for lists'
!    list_left%next => null()
!end subroutine list_assign

! list_create --
!     Create and initialise a list
! Arguments:
!     list       Pointer to new linked list
!     data       The data for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use list_destroy first to
!     destroy up an old list.
!
subroutine list_create( list, data )
    type(LINKED_LIST), pointer  :: list
    type(LIST_DATA), intent(in) :: data

    allocate( list )
    list%next => null()
    list%data =  data
end subroutine list_create

! list_destroy --
!     Destroy an entire list
! Arguments:
!     list       Pointer to the list to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!
subroutine list_destroy( list )
    type(LINKED_LIST), pointer  :: list

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: elem

    elem => list
    do while ( associated(elem) )
        current => elem
        elem => current%next
        deallocate( current )
    enddo
end subroutine list_destroy

! list_count --
!     Count the number of items in the list
! Arguments:
!     list       Pointer to the list
!
integer function list_count( list )
    type(LINKED_LIST), pointer  :: list

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next

    if ( associated(list) ) then
        list_count = 1
        current => list
        do while ( associated(current%next) )
            current => current%next
            list_count = list_count + 1
        enddo
    else
        list_count = 0
    endif
end function list_count

! list_next
!     Return the next element (if any)
! Arguments:
!     elem       Element in the linked list
! Result:
!
function list_next( elem ) result(next)
    type(LINKED_LIST), pointer :: elem
    type(LINKED_LIST), pointer :: next

    next => elem%next

end function list_next

! list_insert
!     Insert a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!
subroutine list_insert( elem, data )
    type(LINKED_LIST), pointer  :: elem
    type(LIST_DATA), intent(in) :: data

    type(LINKED_LIST), pointer :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data
end subroutine list_insert

! list_insert_head
!     Insert a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!
subroutine list_insert_head( list, data )
    type(LINKED_LIST), pointer  :: list
    type(LIST_DATA), intent(in) :: data

    type(LINKED_LIST), pointer :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem
end subroutine list_insert_head

! list_delete_element
!     Delete an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be
!                removed
!
subroutine list_delete_element( list, elem )
    type(LINKED_LIST), pointer  :: list
    type(LINKED_LIST), pointer  :: elem

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: prev

    if ( associated(list,elem) ) then
        list => elem%next
        deallocate( elem )
    else
        current => list
        prev    => list
        do while ( associated(current) )
            if ( associated(current,elem) ) then
                prev%next => current%next
                deallocate( current ) ! Is also "elem"
                exit
            endif
            prev    => current
            current => current%next
        enddo
    endif
!    allocate(next)
!
!    next%next => elem%next
!    elem%next => next
!    next%data =  data
end subroutine list_delete_element

! list_get_data
!     Get the data stored with a list element
! Arguments:
!     elem       Element in the linked list
!
function list_get_data( elem ) result(data)
    type(LINKED_LIST), pointer :: elem

    type(LIST_DATA)            :: data

    data = elem%data
end function list_get_data

! list_put_data
!     Store new data with a list element
! Arguments:
!     elem       Element in the linked list
!     data       The data to be stored
!
subroutine list_put_data( elem, data )
    type(LINKED_LIST), pointer  :: elem
    type(LIST_DATA), intent(in) :: data

    elem%data = data
end subroutine list_put_data


!
! Routines and functions specific to dictionaries
!

! dict_create --
!     Create and initialise a dictionary
! Arguments:
!     dict       Pointer to new dictionary
!     key        Key for the first element
!     value      Value for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use dict_destroy first to
!     destroy up an old list.
!
subroutine dict_create( dict, key, value )
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) ::  key
    type(DICT_DATA), intent(in)  :: value

    type(LIST_DATA)              :: data
    integer                      :: i
    integer                      :: hash

    allocate( dict )
    allocate( dict%table(hash_size) )

    do i = 1,hash_size
        dict%table(i)%list => null()
    enddo

    data%key   = key
    data%value = value

    hash = dict_hashkey( trim(key ) )
    call list_create( dict%table(hash)%list, data )

end subroutine dict_create

! dict_destroy --
!     Destroy an entire dictionary
! Arguments:
!     dict       Pointer to the dictionary to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!
subroutine dict_destroy( dict )
    type(DICT_STRUCT), pointer  :: dict

    integer                     :: i

    do i = 1,size(dict%table)
        if ( associated( dict%table(i)%list ) ) then
            call list_destroy( dict%table(i)%list )
        endif
    enddo
    deallocate( dict%table )
    deallocate( dict )

end subroutine dict_destroy

! dict_add_key
!     Add a new key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for the new element
!     value      Value for the new element
! Note:
!     If the key already exists, the
!     key's value is simply replaced
!
subroutine dict_add_key( dict, key, value )
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key
    type(DICT_DATA), intent(in)  :: value

    type(LIST_DATA)              :: data
    type(LINKED_LIST), pointer   :: elem
    integer                      :: hash

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
        elem%data%value = value
    else
        data%key   = key
        data%value = value
        hash       = dict_hashkey( trim(key) )
        if ( associated( dict%table(hash)%list ) ) then
            call list_insert( dict%table(hash)%list, data )
        else
            call list_create( dict%table(hash)%list, data )
        endif
    endif

end subroutine dict_add_key

! dict_delete_key
!     Delete a key-value pair from the dictionary
! Arguments:
!     dict       Dictionary in question
!     key        Key to be removed
!
subroutine dict_delete_key( dict, key )
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key

    type(LINKED_LIST), pointer   :: elem
    integer                      :: hash

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
        hash = dict_hashkey( trim(key) )
        call list_delete_element( dict%table(hash)%list, elem )
    endif
end subroutine dict_delete_key

! dict_get_key
!     Get the value belonging to a key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for which the values are sought
!
function dict_get_key( dict, key ) result(value)
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key
    type(DICT_DATA)              :: value

    type(LIST_DATA)              :: data
    type(LINKED_LIST), pointer   :: elem

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
        value = elem%data%value
    else
        value = DICT_NULL
    endif
end function dict_get_key

! dict_has_key
!     Check if the dictionary has a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!
function dict_has_key( dict, key ) result(has)
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key
    logical                      :: has

    type(LINKED_LIST), pointer   :: elem

    elem => dict_get_elem( dict, key )

    has = associated(elem)
end function dict_has_key

! dict_get_elem
!     Find the element with a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!
function dict_get_elem( dict, key ) result(elem)
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key

    type(LINKED_LIST), pointer   :: elem
    integer                      :: hash

    hash = dict_hashkey( trim(key) )

    elem => dict%table(hash)%list
    do while ( associated(elem) )
        if ( elem%data%key .eq. key ) then
            exit
        else
            elem => list_next( elem )
        endif
    enddo
end function dict_get_elem

! dict_hashkey
!     Determine the hash value from the string
! Arguments:
!     key        String to be examined
!
integer function dict_hashkey( key )
    character(len=*), intent(in) :: key

    integer                      :: hash
    integer                      :: i

    dict_hashkey = 0

    do i = 1,len(key)
        dict_hashkey = modulo( multiplier * dict_hashkey + ichar(key(i:i)), hash_size )
    enddo

    dict_hashkey = 1 + modulo( dict_hashkey-1, hash_size )
end function dict_hashkey

end module Element_Dict
