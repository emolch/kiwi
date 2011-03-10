! 
!    Copyright 2011 Sebastian Heimann
! 
!    Licensed under the Apache License, Version 2.0 (the "License");
!    you may not use this file except in compliance with the License.
!    You may obtain a copy of the License at
! 
!        http://www.apache.org/licenses/LICENSE-2.0
! 
!    Unless required by applicable law or agreed to in writing, software
!    distributed under the License is distributed on an "AS IS" BASIS,
!    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!    See the License for the specific language governing permissions and
!    limitations under the License.
!

module unit

    ! This module provides two subroutines: claim_unit( u ) and release_unit( u ) to
    ! get an unused unit number and set it free again.

    ! The unit numbers obtained by claim_unit( u ) are (currently) in the range [10000,90000]
    ! so the calling program should never use any unit number in this range explicitly.

    ! It also provides the constants stderr, stdin and stdout, wich should be set to the compiler's
    ! error, in and out preconnected units.
    
    ! LIMITATIONS
    ! (should be OK for all normal purposes, 
    !  hardly anybody needs 10000 open files at once, and chance is, that
    !  the OS won't allow it anyway...)
    ! This will fail, after 10001 units have been claimed (but not released).
    ! Performance may become bad when lots of units have been claimed (but not released).
    
    use util
    
    implicit none
    
    ! access specifiers
    
    private
    public :: claim_unit, release_unit 
    
    ! type defs
    
    type, private :: unit_t
        private
        integer :: number
        type(unit_t), pointer :: next
    end type
    
    ! globals constants
        
    integer, dimension(2), private, parameter :: unit_range = (/10000,20000/)

    ! global variables
    
    integer, private :: current = 10000
    type(unit_t), pointer :: used_units => null()
    
  contains
    
    subroutine claim_unit(u)
        
        integer, intent(out) :: u
        integer :: count
        
        ! get an unused unit number
        
        count = 0
        do while ( isin( used_units, current ) )
            current = current + 1
            if (current > unit_range(2)) current = unit_range(1)
            count = count+1
            if ( count > (unit_range(2)-unit_range(1)+1) ) then
                stop 'unit.f90: too many claimed units'
            end if
        end do
        call push( used_units, current )
        u = current
        
    end subroutine claim_unit

    subroutine release_unit(u)
      
        integer, intent(in) :: u
      
        ! free a unit number
        
        call remove_by_val( used_units, u )
 
    end subroutine release_unit


! stuff below here should be in separate module
! cause it is only a very basic (and slow) implementation of a list container
! should be provided by some standard library

    subroutine push( list, i )
        
        type(unit_t), pointer :: list
        integer, intent(in) :: i
        type(unit_t), pointer :: current, new
        
        ! add an item to the end of the list
        
        allocate( new )
        new%number = i
        new%next => null()

        if (.not. associated( list ) ) then
            list => new
        else 
            current => list
            do while (associated(current%next))
                current => current%next
            end do
            current%next => new
        end if
        
    end subroutine push

    subroutine remove_by_val( list, i )
    
        type(unit_t), pointer :: list
        integer, intent(in) :: i
        
        ! remove all items with value i from the list
        
        type(unit_t), pointer :: current, del
        
        if ( .not. associated(list) ) return
        
        do while (list%number == i)
            current => list
            list => current%next
            deallocate(current)
            if ( .not. associated(list) ) exit
        end do
        
        if ( .not. associated(list) ) return
        
        current => list
        do while (associated(current%next))
            if (current%next%number == i) then
                del => current%next
                current%next => current%next%next
                deallocate( del )
            else
                current => current%next
            end if
        end do
    
    end subroutine remove_by_val
    
    function isin( list, i ) result(yn)
    
        type(unit_t), pointer :: list
        integer, intent(in) :: i
        logical :: yn
            
        ! check if i is in list
         
        type(unit_t), pointer :: current
        
        current => list
        yn = .false.
        
        do while ( associated(current) )
            if ( current % number == i ) then
                yn = .true.
                return
            end if
            current => current%next
        end do 
        
    end function isin
    
end module unit
