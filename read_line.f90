! $Id: read_line.f90 658 2007-08-03 12:48:49Z sebastian $ 
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



module read_line
    
  ! read arbitrarily long lines of input into character string
  
  ! skips comment lines starting with #
  ! # may be preceeded by whitespace
  
  ! because fortran can not directly handle character strings of varying length,
  ! callbacks are used such that a routine may be specified, wich will be called
  ! from read_line() with a character string of sufficient length, containing the current line.
  
  ! has been tested with ifort (Version 9.0) and F (Fortran Company/NAG F compiler Release 20031017)

  ! with ifort lines with more than a million characters work,
  ! a bug in the current version of F limits the line length to about 10000 characters.
  
    use util
  
    implicit none
    
    integer, parameter, private :: INITIAL_BUFFER_LEN = 1024
    integer, parameter, private :: MAX_READ_LEN       = 512
    integer, parameter, private :: BUFFER_MULT        = 4
     
    public  :: readline
    private :: readline_, det_buf_len, fill_buffer, is_comment
    
    
  contains
    
    subroutine readline( readsub, iostat, ok, unit )
    
        interface
            subroutine readsub( buffer, ok )
                character(len=*), intent(in) :: buffer
                logical, intent(out)         :: ok
            end subroutine readsub
        end interface
                
        integer, intent(out)            :: iostat
        logical, intent(out)            :: ok
        integer, intent(in), optional   :: unit

      ! read one line from * or unit
      ! eof condition is returned in iostat.
      ! readsub is a user function that will be called with an appropriatly large buffer
      ! to hold the complete line.
      ! readsub may return ok=true if a read from the string has worked and ok=false if not
      ! this will be passed back to the calling function
   
        if ( present(unit) ) then
            call readline_( "", readsub, iostat, ok, unit )
        else
            call readline_( "", readsub, iostat, ok )
        end if
      
    end subroutine readline
    
    recursive subroutine readline_( prev_buffer, readsub, iostat, ok, unit  )
        
        interface
            subroutine readsub( buffer, ok )
                character(len=*), intent(in) :: buffer
                logical, intent(out)         :: ok
            end subroutine readsub
        end interface
        
        character(len=*), intent(in)    :: prev_buffer
        integer, intent(out)            :: iostat
        logical, intent(out)            :: ok
        integer, intent(in), optional   :: unit
        
        character(len=det_buf_len(len(prev_buffer))) :: buffer
        integer                                      :: pblen, blen, n_chars_read
                     
        ok = .false.
        pblen = len(prev_buffer)
        blen = len(buffer)
        buffer(:pblen) = prev_buffer
        
        if ( present(unit) ) then
            call fill_buffer( buffer(pblen+1:blen), iostat, n_chars_read, unit )
        else
            call fill_buffer( buffer(pblen+1:blen), iostat, n_chars_read )
        end if
        
        ! if we are still not at the end try with an even bigger buffer
        if (iostat == 0) then
            if ( present(unit) ) then
                call readline_( buffer, readsub, iostat, ok, unit )
            else
                call readline_( buffer, readsub, iostat, ok )
            end if
            return
        end if
        if (iostat > 0) return
                
        if (iostat == IOSTAT_EOR) then
            if (is_comment( buffer, pblen+n_chars_read )) then
                ok = .true.
                return
            end if
            call readsub( trim(buffer), ok )
        end if 
    end subroutine readline_
    
    
    pure function det_buf_len( prev_buf_len ) result( buf_len )
        
        integer, intent(in) :: prev_buf_len
        integer             :: buf_len
        
        ! used by readline_ to determine new buffer length
        
        if (prev_buf_len == 0) then
            buf_len = INITIAL_BUFFER_LEN
        else
            buf_len = prev_buf_len*BUFFER_MULT
        end if
        
    end function det_buf_len
    
    
    subroutine fill_buffer( buffer, iostat, nfill, unit )
    
        character(len=*), intent(inout) :: buffer
        integer, intent(out)            :: iostat, nfill
        integer, intent(in), optional   :: unit
        
        ! fills the buffer chunk by chunk
        
        integer :: nwant, nchars
        
        nfill = 0
        nwant = len(buffer)
        fill_loop : do
            nchars = min(nwant,MAX_READ_LEN)
            if ( present(unit) ) then
                read (unit=unit, fmt="(A)", advance="NO", iostat=iostat, size=nchars) buffer(nfill+1:)
            else
                read (unit=*, fmt="(A)", advance="NO", iostat=iostat, size=nchars) buffer(nfill+1:)
            end if
            nfill = nfill + nchars
            nwant = nwant - nchars
            
            if (iostat > 0) iostat = 0
            if (nwant == 0) return
            if (iostat < 0) return 
        end do fill_loop
        
    end subroutine fill_buffer
    
    function is_comment( line, last ) result(comment)
        
        character(len=*), intent(in) :: line
        integer, intent(in)          :: last
        logical                      :: comment
        
      ! determine if line is a comment line
        
        character(len=1)             :: c
        integer                      :: i
        
        comment = .false.
        if (last == 0) then
            comment = .true.
            return
        end if
        c = " "
        ! go to first non-blank char
        char_loop : do i=1,last
            c = line(i:i)
            if (c /= " " .and. c /= char(9) .and. c /= char(10) .and. c /= char(13)) exit char_loop
        end do char_loop
        if (c == "#" .or. i-1 == last) comment = .true.
        return
        
    end function is_comment
        
end module read_line
