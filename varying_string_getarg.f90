! $Id: varying_string_getarg.f90 681 2007-11-30 15:28:54Z sebastian $ 
! ------------------------------------------------------------------------------
! 
!    Copyright 2007 Sebastian Heimann
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

module varying_string_getarg

    use better_varying_string
    ! use f90_unix_env
    
    implicit none
    
    private
    public vs_getarg
    public vs_getenv
    
    interface vs_getenv
        module procedure vs_getenv_c
        module procedure vs_getenv_vs
    end interface
    
    public real_getarg
    public int_getarg
    
   integer, parameter :: max_arg_len = 100000

  contains
    
    subroutine vs_getarg( ipos, vs )
        
        integer, intent(in)               :: ipos
        type(varying_string), intent(inout) :: vs
    
    ! get argument at position ipos as varying_string
    ! argument 0 is the command
        
        integer                    :: leng
        character(len=max_arg_len) :: dummy
        
        ! get length of the argument
        call getarg( ipos, dummy )
        leng = len_trim(dummy)
        
        ! extract it to vs
        call extract_arg( ipos, leng, vs )
        
    end subroutine
    
    subroutine vs_getenv_c( key, vs )
        
        character(len=*), intent(in)    :: key
        type(varying_string), intent(inout) :: vs
    
    ! get argument at position ipos as varying_string
    ! argument 0 is the command
        
        integer                    :: leng
        character(len=max_arg_len) :: dummy
        
        ! get length of the argument
        call getenv( key, dummy )
        leng = len_trim(dummy)
        
        ! extract it to vs
        call extract_env( key, leng, vs )
        
    end subroutine

    subroutine vs_getenv_vs( key, vs )
        
        type(varying_string), intent(in)    :: key
        type(varying_string), intent(inout) :: vs
    
        call vs_getenv( char(key), vs )
        
    end subroutine

            
    subroutine extract_arg( ipos, length, vs )
        
        integer, intent(in)   :: ipos, length
        type(varying_string), intent(inout) :: vs
        character(len=length) :: buffer
        
        call getarg( ipos, buffer )
        vs = buffer
        
    end subroutine
    
    subroutine extract_env( key, length, vs )
        
        character(len=*), intent(in)    :: key
        integer, intent(in)                 :: length
        type(varying_string), intent(inout) :: vs
        character(len=length) :: buffer
        
        call getenv(key, buffer )
        vs = buffer
        
    end subroutine
    
    subroutine real_getarg( iarg, mini, maxi, val )
    
        integer, intent(in) :: iarg
        real, intent(in) :: mini, maxi
        real, intent(out) :: val
    
    ! reads argument iarg from the argument list and converts it to real number val
    ! it is checked, that val is in the range [min, max]
        
        type(varying_string) :: vs
        
        call vs_getarg( iarg, vs )
        val = vs
        if (val < mini) val = mini
        if (val > maxi) val = maxi
        
    end subroutine
    
    subroutine int_getarg( iarg, mini, maxi, val )
    
        integer, intent(in) :: iarg
        integer, intent(in) :: mini, maxi
        integer, intent(out) :: val
    
    ! reads argument iarg from the argument list and converts it to real number val
    ! it is checked, that val is in the range [min, max]
        
        type(varying_string) :: vs
        
        call vs_getarg( iarg, vs )
        val = vs
        
        if (val < mini) val = mini
        if (val > maxi) val = maxi
        
    end subroutine
    
end module
