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

module util

    use better_varying_string
    
    implicit none
    private
    
  ! globals constants
    
  ! numbers of the preconnected units. 
  ! these are at least valid for the intel fortran compiler.
    integer, public, parameter :: stderr = 0
    integer, public, parameter :: stdin  = 5
    integer, public, parameter :: stdout = 6
    
  ! iostat conditions
    integer, public, parameter :: iostat_eof = -1
    integer, public, parameter :: iostat_eor = -2

    integer, public, parameter :: sig_ignore = 1
    integer, public, parameter :: sig_default = 0
    integer, public, parameter :: sig_int = 2
    
  ! global variables
  
  ! this flag indicates that a test has failed
    integer, public :: g_fail = 0
  
  ! set by test_begin
    real, public :: g_time_begin = 0.0
 
  ! the name of the running program should be stored in here:
  ! this is for the die() and test_*() subs
    type(varying_string), public :: g_pn
    
  ! a usage string should be put in here, so that usage() works
    type(varying_string), public :: g_usage
    
  ! this holds a string describing the last error
    type(varying_string), public :: g_errstr
  
  ! this turns on printing of inform() messages
    logical, public :: g_verbose = .false.
    
    public die, usage, cleanup, count_words, count_non_comment_lines
    public whitespace, resize, error, warn, inform, skip_comments
    public test_begin, test_end, test_fail

    interface resize
        module procedure resize_r
        module procedure resize_d
        module procedure resize_i
        module procedure resize_c
        module procedure resize_l
    end interface
    
    interface die
        module procedure die_vs
        module procedure die_c
        module procedure die_void
    end interface
    
    interface warn
        module procedure warn_vs
        module procedure warn_c
        module procedure warn_void
    end interface
    
    interface inform
        module procedure inform_vs
        module procedure inform_c
    end interface
    
    interface test_begin
        module procedure test_begin_vs
        module procedure test_begin_c
    end interface
    
    interface test_fail
        module procedure test_fail_vs
        module procedure test_fail_c
    end interface
    
    interface error
        module procedure error_vs
        module procedure error_c
    end interface
    
    
contains 

    subroutine error_vs( str )
        type(varying_string),intent(in) :: str
        call error_c( char(str) )
    end subroutine
    
    subroutine error_c( c )
        character(len=*), intent(in) :: c
        g_errstr = c
    end subroutine
    
    subroutine die_vs( str )
        type(varying_string),intent(in) :: str
        call die_c( char(str) )
    end subroutine
    
    subroutine die_c( c )
        character(len=*), intent(in) :: c
        write (stderr,fmt="(a,a,a)") char(g_pn), ": ", c        
        stop 1
    end subroutine
    
    subroutine die_void()
        call die_vs( g_errstr )    
    end subroutine
    
    subroutine warn_vs( str )
        type(varying_string),intent(in) :: str
        call warn_c( char(str) )
    end subroutine
    
    subroutine warn_c( c )
        character(len=*), intent(in) :: c
        write (stderr,fmt="(a,a,a)") char(g_pn), "(warn): ", c
    end subroutine
    
    subroutine warn_void()
        call warn_vs( g_errstr )    
    end subroutine
    
    subroutine inform_vs( str )
        type(varying_string),intent(in) :: str
        call inform_c( char(str) )
    end subroutine
    
    subroutine inform_c( c )
        character(len=*), intent(in) :: c
        if (g_verbose) then
            write (stderr,fmt="(a,a,a)") char(g_pn), "(info): ", c
        end if
    end subroutine
    
    subroutine usage()
        write (stderr, fmt="(a)") char(g_usage)
        stop 1
    end subroutine

    subroutine cleanup()
        call delete(g_pn)
        call delete(g_usage)
        call delete(g_errstr)
    end subroutine
    
    subroutine test_begin_vs( str )
        type(varying_string),intent(in) :: str
        call test_begin_c( char(str) )
    end subroutine
    
    subroutine test_begin_c( c )
        
        character(len=*), intent(in) :: c
        
        call cpu_time( g_time_begin )
        g_pn = c
        g_fail = 0
        
    end subroutine
    
    subroutine test_fail_vs( str )    
        type(varying_string),intent(in) :: str
        call test_fail_c( char(str) )
    end subroutine
    
    subroutine test_fail_c( c )
        
        character(len=*), intent(in) :: c
        
        write (stderr,fmt="(a,a,a,a)") "fail: ",char(g_pn),": ",c
        
        g_fail = g_fail+1
        
    end subroutine
    
    subroutine test_end()
                
        real :: time_end
        call cpu_time( time_end )
        
        if (g_fail >0) then
            if (g_fail > 1) then
                write (stderr,fmt="(a)") char(g_pn // ": " // g_fail //" tests failed." )
            else 
                write (stderr,fmt="(a)") char(g_pn // ": one test failed.")
            end if
        else 
            write (stderr,fmt="(a)") char(g_pn //": OK, time: " //  (time_end-g_time_begin) // " s" )
        end if
        
    end subroutine
    
    pure function count_words( buffer ) result( nwords )
        character(len=*), intent(in) :: buffer
        integer                      :: nwords
        
        logical :: atword
        integer :: i
                
        nwords = 0
        atword = .false.
        do i=1,len_trim(buffer)
            if ( .not. atword ) then
                if ( .not. whitespace(buffer(i:i)) ) then
                    ! at the beginning of a word
                    atword = .true.
                    nwords = nwords + 1
                end if
            else
                if ( whitespace( buffer(i:i) ) ) then
                    atword = .false.
                end if
            end  if   
        end do
        

    end function count_words
    
    function count_non_comment_lines( iunit, iostat ) result( nlines )
    
        integer, intent(in)  :: iunit
        integer, intent(out) :: iostat
        integer :: nlines, nskip
    
      ! count non comment lines in file iunit
      ! it rewinds the file after counting
        
        nlines = 0
        iostat = 0
        
        line_loop : do
            
            nskip = skip_comments(iunit,iostat)
            if (iostat == IOSTAT_EOF .or. nskip < 0) exit line_loop
            if (iostat /= 0) return
        
            read (iunit,*,iostat=iostat)
            if (iostat == IOSTAT_EOF) exit line_loop
            if (iostat /= 0) return
      
            nlines = nlines+1
       
        end do line_loop
        
        rewind( iunit, iostat=iostat )
    
    end function
    
    
    function skip_comments( iunit, iostat )
    
        integer, intent(in)  :: iunit
        integer, intent(out) :: iostat
        integer              :: skip_comments
        
      ! this skips comment lines (empty lines and lines starting with #) on iunit
      ! the number of skipped lines is returned
      ! the last io commands iostat is returned via iostat
      !
      ! example usage:
      ! 
      ! line_loop : do
      !     nskip  = skip_comments(iunit,iostat)
      !     if (iostat == IOSTAT_EOF) exit line_loop
      !     if (iostat /= 0) stop "ioerror occured while skipping comments"
      !  
      !     read (iunit,*,iostat=iostat) data
      !     if (iostat == IOSTAT_EOF) exit line_loop
      !     if (iostat /= 0) stop "read failed"
      !     
      !      call do_sth_with_data()
      !  end do line_loop
      !
              
        character :: ch
        
        iostat = 0
        skip_comments = 0
        
        line_loop : do
            
            char_loop : do
                read (iunit,'(a)',iostat=iostat,advance='no') ch
                if (iostat == iostat_eof) return 
                if (iostat == iostat_eor) then
                   ! print *, "skipping a comment line (empty line)"
                    skip_comments = skip_comments+1
                    cycle line_loop
                end if
                if (iostat /= 0) return
            
                if (ch == '#') then
                   ! print *, "skipping a comment line"
                    skip_comments = skip_comments+1
                    read (iunit,*, iostat=iostat) ! move to next record
                    if (iostat == iostat_eof) return
                    cycle line_loop
                end if
                if (.not. whitespace(ch)) then
                    backspace(iunit,iostat=iostat) ! go back to beginning of current record
                    return
                end if
            end do char_loop
        
        end do line_loop
        
    end function

    pure function whitespace( c ) result( isws )
        character(len=1), intent(in) :: c
        logical                      :: isws
        isws = c == " " .or. c == char(9) .or. c == char(10) .or. c == char(13) .or. c == char(0)
    end function whitespace
    
    pure subroutine resize_r( array, offset, length )
        
        real, intent(inout), dimension(:), allocatable :: array
        integer, intent(in) :: offset, length
        
        if ( .not. allocated( array ) ) then
            if (length == 0) return
            allocate( array(offset:offset+length-1) )
            return
        end if
        
        if ( size(array) /= length .or. lbound(array,1) /= offset ) then
            deallocate( array )
            if (length /= 0) then
                allocate( array(offset:offset+length-1) )
            end if
        end if
        
    end subroutine
    
    pure subroutine resize_d( array, offset, length )
        
        real(kind=8), intent(inout), dimension(:), allocatable :: array
        integer, intent(in) :: offset, length
        
        if ( .not. allocated( array ) ) then
            if (length == 0) return
            allocate( array(offset:offset+length-1) )
            return
        end if
        
        if ( size(array) /= length .or. lbound(array,1) /= offset ) then
            deallocate( array )
            if (length /= 0) then
                allocate( array(offset:offset+length-1) )
            end if
        end if
        
    end subroutine
    
    pure subroutine resize_i( array, offset, length )
        
        integer, intent(inout), dimension(:), allocatable :: array
        integer, intent(in) :: offset, length
        if ( .not. allocated( array ) ) then
            if (length == 0) return
            allocate( array(offset:offset+length-1) )
            return
        end if
        
        if ( size(array) /= length .or. lbound(array,1) /= offset ) then
            deallocate( array )
            if (length /= 0) then
                allocate( array(offset:offset+length-1) )
            end if
        end if
        
    end subroutine
    
    pure subroutine resize_c( array, offset, length )
        
        complex, intent(inout), dimension(:), allocatable :: array
        integer, intent(in) :: offset, length
        
        if ( .not. allocated( array ) ) then
            if (length == 0) return
            allocate( array(offset:offset+length-1) )
            return
        end if
        
        if ( size(array) /= length .or. lbound(array,1) /= offset ) then
            deallocate( array )
            if (length /= 0) then
                allocate( array(offset:offset+length-1) )
            end if
        end if
        
    end subroutine
    
    pure subroutine resize_l( array, offset, length )
        
        logical, intent(inout), dimension(:), allocatable :: array
        integer, intent(in) :: offset, length
        
        if ( .not. allocated( array ) ) then
            if (length == 0) return
            allocate( array(offset:offset+length-1) )
            return
        end if
        
        if ( size(array) /= length .or. lbound(array,1) /= offset ) then
            deallocate( array )
            if (length /= 0) then
                allocate( array(offset:offset+length-1) )
            end if
        end if
        
    end subroutine


end module

