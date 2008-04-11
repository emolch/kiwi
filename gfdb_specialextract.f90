! $Id: gfdb_specialextract.f90 623 2007-04-25 14:49:07Z sebastian $ 
! ------------------------------------------------------------------------------
module gfdb_specialextract_

    use util
    use gfdb
    use better_varying_string
    use sparse_trace
    use seismogram_io   
    
    implicit none
    
    private
    
    public getentry, set_database, close_database
    
    type(t_gfdb), save :: db
  
  contains
  
    subroutine set_database( db_path )
        type(varying_string), intent(in)      :: db_path
        
        call gfdb_init(db, db_path)
        
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
        integer :: ix, iz, ig, iostat, nerr, istrip
        real    :: x, z
        integer, dimension(2)        :: span
        type(varying_string)         :: filename
        type(t_strip)                :: continuous
   
        read (unit=buffer,fmt=*,iostat=iostat) x, filenamebase
        if (iostat /= 0) return
        ok = .true.
        z = 0.0
        call gfdb_get_indices( db, x, z, ix, iz )
        tracep => null()
        do iz=1,db%nz
            do ig=1,db%ng
        
                call gfdb_get_trace( db, ix, iz, ig, tracep )
                if (.not. associated(tracep)) cycle
                if ( trace_is_empty(tracep)) cycle
                print *, ig, iz
                do istrip=1,tracep%nstrips
                    tracep%strips(istrip)%data = tracep%strips(istrip)%data**2
                end do

                call trace_multiply_add( tracep, continuous )
                call gfdb_uncache_trace( db, ix,iz,ig )
            end do
        end do
        
        continuous%data = sqrt( continuous%data )
        span = strip_span( continuous )
        filename = trim(filenamebase)
        call writeseismogram( filename, var_str("*"), continuous%data, db%dt*(span(1)-1),db%dt,nerr )
                
    end subroutine
    
end module


program gfdb_specialextract

  ! extract sum over depths of sqrt(g_1^2+ ... + g_8^2)
  ! 
  ! this is not something physically meaningful, 
  ! i only want to use this to pick time windows
  !
  ! usage: gfdb_specialextract database
  
  ! it will then take look at stdin and scan for lines like this and save
  ! the so selected seismogram to that file
  
  !   x 'filenamebase'

  ! where ix is the receiver number
  ! filenamebase is used to generate output files.
  !              one file for every time window
  

    use util
    use better_varying_string
    use varying_string_getarg
    use read_line
    use gfdb_specialextract_
    ! use f90_unix_env

    implicit none    
    
    type(varying_string) :: basefn
    integer              :: iostat, iline
    logical              :: ok
    character, parameter :: eol = char(10)
    
    g_pn = 'gfdb_specialextract'
    g_usage = 'usage: ' // g_pn // ' database <<EOF' // eol // & 
              "x ig 'filename'" // eol // &
              "..." // eol // &
              "EOF" // eol // eol // &
              "This program is not documented and I will never write any documentation on it."

    if (iargc() /= 1) call usage()

    call vs_getarg( 1, basefn )
    
    call set_database( basefn )

    iline = 1
    line_loop : do
        call readline( getentry, iostat, ok )
        if (iostat == IOSTAT_EOF) exit line_loop
        if (.not. ok) call die( g_pn // ": parsing of line "// iline &
                                 //" from standard input failed."// &
                                 "expected: x ig 'filename'" )
        iline = iline+1
    end do line_loop

    call close_database()

    call cleanup()
    
end program


