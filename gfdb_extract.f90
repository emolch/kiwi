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

module gfdb_extract_

    use util
    use gfdb
    use better_varying_string
    use sparse_trace
    use seismogram_io   
    
    implicit none
    
    private
    
    public getentry, set_database, close_database
    
    type(t_gfdb), save :: db
    integer :: traces_added = 0
  
  contains
  
    subroutine set_database( db_path, nipx, nipz )
        type(varying_string), intent(in)      :: db_path
        integer, intent(in) :: nipx, nipz
        call gfdb_init(db, db_path, nipx=nipx, nipz=nipz )
        
    end subroutine

    subroutine close_database()
        call gfdb_destroy( db )
    end subroutine
    
    subroutine getentry( buffer, ok )
    
      ! called by readline for every non-comment line found
      
        character(len=*), intent(in) :: buffer
        logical, intent(out)         :: ok
        
        character(len=len(buffer))   :: filenamebase
        type(t_trace), pointer       :: tracep
        integer :: ix, iz, ig, iostat, nerr
        real    :: x, z
        integer, dimension(2)        :: span
        type(varying_string)         :: filename
        type(t_strip)                :: continuous
        tracep => null()
   
        ok = .false.
        read (unit=buffer,fmt=*,iostat=iostat) x, z, ig, filenamebase
        if (iostat /= 0) return
        call gfdb_get_indices( db, x, z, ix, iz )
        call gfdb_get_trace( db, ix, iz, ig, tracep )
        if (.not. associated(tracep)) return
        if ( trace_is_empty(tracep)) return
        call trace_unpack( tracep, continuous )
        span = strip_span( continuous )
        filename = trim(filenamebase)
        call writeseismogram( filename, var_str("*"), continuous%data, dble(db%dt*(span(1)-1)),db%dt, &
            var_str(''), var_str(''), var_str(''), var_str(ig), &
            nerr )
        if (nerr == 0) ok = .true.

        call gfdb_uncache_trace( db, ix,iz,ig )
        traces_added = traces_added + 1

        ! periodically close the gfdb, so that hdf is forced to deallocate
        ! all it's memory
        if (traces_added > 8000) then
            call gfdb_close( db )
            traces_added = 0
        end if
        
    end subroutine
    
end module


program gfdb_extract

  ! This program is used to extract individual traces from a Greens function database created with gfdb_build.
  ! 
  ! usage: gfdb_extract database [ nipx nipz ] <<EOF
  ! x z ig 'filename'
  ! EOF
  !
  ! complete documentation available on 
  !
  !   http://kinherd.org/power/trac/wiki/GFDBExtractTool
  !

    use util
    use better_varying_string
    use varying_string_getarg
    use read_line
    use gfdb_extract_
    ! use f90_unix_env

    implicit none    
    
    type(varying_string) :: basefn
    integer              :: iostat, iline
    logical              :: ok
    character, parameter :: eol = char(10)
    integer              :: nipx, nipz
    
    g_pn = 'gfdb_extract'
    g_usage = 'usage: ' // g_pn // ' database [ nipx nipz ]<<EOF' // eol // & 
              "x z ig 'filename'" // eol // &
              "..." // eol // &
              "EOF" // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/GFDBExtractTool"

    if (iargc() /= 1 .and. iargc() /= 3) call usage()

    call vs_getarg( 1, basefn )
    nipx = 1
    nipz = 1
    if (iargc() == 3) then 
        call int_getarg( 2, 1, huge(nipx), nipx )
        call int_getarg( 3, 1, huge(nipz), nipz )
    end if
    
    call set_database( basefn, nipx, nipz )

    iline = 1
    line_loop : do
        call readline( getentry, iostat, ok )
        if (iostat == IOSTAT_EOF) exit line_loop
        if (.not. ok) then
            write (stdout,*) 'nok'
        else
            write (stdout,*) 'ok'
        end if
        call flush(stdout)
        iline = iline+1
    end do line_loop

    call close_database()

    call cleanup()
    
    call delete(basefn)
    
end program
  


