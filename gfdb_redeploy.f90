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

module gfdb_redeploy_

    use util
    use gfdb
    use better_varying_string
    use sparse_trace
    
    implicit none
      
    private  
    public copyentry
    public set_database_in, set_database_out_nocreate
    public close_databases
    public redeploy_cleanup
    public set_mapping
    
    type(t_gfdb), save :: in, out
    integer :: traces_added = 0
    integer, dimension(2) :: last_span = (/1,1/)
    integer, dimension(:), allocatable :: mapping
    
  contains
  
    

    subroutine set_database_in( basefn, nipx, nipz )
        
        type(varying_string), intent(in) :: basefn
        integer, intent(in) :: nipx, nipz
        integer :: ig
        
        call gfdb_init(in, basefn, nipx=nipx, nipz=nipz )
        allocate(mapping(in%ng))
         do ig=1,in%ng
            mapping(ig) = ig
        end do

    end subroutine
  
    subroutine set_mapping(igs, igt)
        integer, intent(in) :: igs, igt
        if (igs < 1 .or. igs > size(mapping)) then
            call die("invalid mapping "//igs//" -> "//igt)
        end if 
        mapping(igs) = igt
    end subroutine
  
! not needed anymore, I think the -nocreate option was too dangerous
! so now, gfdb_redeploy will never try to create the database.    
!     subroutine set_database_out( basefn )
!         
!         type(varying_string), intent(in) :: basefn
!        
!         call gfdb_init(out, basefn, in%nchunks, in%nx,in%nz,in%ng, &
!                            in%dt,in%dx,in%dz, in%firstx, in%firstz )
!         
!     end subroutine
    
    subroutine set_database_out_nocreate( basefn )
        
        type(varying_string), intent(in) :: basefn
       
        call gfdb_init(out, basefn)
        
    end subroutine
    
    subroutine close_databases()
        call gfdb_destroy( in )
        call gfdb_destroy( out )
    end subroutine

    subroutine copyentry( buffer, ok )
         
        character(len=*), intent(in) :: buffer
        logical, intent(out)         :: ok
        
        type(t_trace), pointer       :: tracep
        integer :: ix, iz, ig, iostat, ixo, izo, nargs
        real :: x,z, tbeg, tend, factor
        type(t_strip) :: strip
        type(t_trace) :: short_trace
        integer, dimension(2) :: span
        tracep => null()
        

        factor = 1.0
        nargs = count_words( buffer )
        if (nargs == 2) then
            read (unit=buffer,fmt=*,iostat=iostat) x, z
        else if (nargs == 3) then
            read (unit=buffer,fmt=*,iostat=iostat) x, z, factor 
        else if (nargs == 4) then
            read (unit=buffer,fmt=*,iostat=iostat) x, z, tbeg, tend
            if (tbeg > tend) return
        end if

        if (iostat == 0) then
            ok = .true.
            
            do ig=1,in%ng
                call gfdb_get_indices( in, x, z, ix, iz )
                call gfdb_get_trace( in, ix, iz, ig, tracep )
                if (factor /= 1.0) then
                    call trace_scale( tracep, factor )
                end if

                if (associated(tracep)) then
                    if ( .not. trace_is_empty(tracep)) then
                        call gfdb_get_indices( out, x, z, ixo, izo )
                        if (nargs /= 4) then
                            if (mapping(ig) /= 0 .and. mapping(ig) <= out%ng) then
                                call gfdb_save_trace( out, ixo, izo, mapping(ig), tracep )
                            end if
                            last_span = tracep%span
                                
                        else 
                            
                            span(1) = max(floor(tbeg/out%dt), tracep%span(1))
                            span(2) = min(ceiling(tend/out%dt), tracep%span(2))
                            
                            call resize( strip%data, span(1), span(2)-span(1)+1 )
                            strip%data = 0.

                            call trace_multiply_add_nogrow( tracep, strip%data, span )

                            if (traces_added > 1) then
                                call trace_pack( strip, short_trace, last_span )
                            else 
                                call trace_pack( strip, short_trace )
                            end if
                            if (mapping(ig) /= 0 .and. mapping(ig) <= out%ng) then
                                call gfdb_save_trace( out, ixo, izo, mapping(ig), short_trace )
                            end if
                            last_span = short_trace%span
                            call strip_destroy( strip )
                        end if
                        traces_added = traces_added + 1
                    end if
                end if
                
                call gfdb_uncache_trace( in, ix,iz,ig )
            end do
            
            call trace_destroy( short_trace )
            
          ! periodically close the gfdb, so that hdf is forced to deallocate
          ! all it's memory
            if (traces_added > 1000) then
                call gfdb_close( out )
                call gfdb_close( in )
                traces_added = 0
            end if
        
        end if
        

    end subroutine
 
    subroutine redeploy_cleanup()
        call close_databases()
        if (allocated(mapping)) deallocate(mapping)
    end subroutine
    
end module


program gfdb_redeploy

  ! This program copies selected traces from one database into another.
  ! 
  ! usage: gfdb_redeploy input-database [ nipx nipz ] output-database <<EOF
  ! x z
  ! ...
  ! EOF
  !
  ! Complete documentation is available on
  !
  !   http://kinherd.org/power/trac/wiki/GFDBRedeployTool
  !
  
    use util
    use better_varying_string
    use varying_string_getarg
    use read_line
    use gfdb_redeploy_
    
    ! use f90_unix_env

    implicit none    
    
    type(varying_string) :: basefn_in, basefn_out
    integer              :: iostat, iline
    logical              :: ok
    character, parameter :: eol = char(10)
    integer              :: nipx, nipz, igt, ig
    
    g_pn = 'gfdb_redeploy'
    g_usage = 'usage: ' // g_pn // ' input-database [ nipx nipz [ g1_mapping g2_mapping ... ] ] output-database <<EOF' // eol // &
              'x z [ tbeg tend ]' // eol // &
              '...' // eol // &
              'EOF' // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/GFDBRedeployTool"
    
              
  ! this argument processing definitely needs cleanup.
  ! should implement a getopt(), but this is no fun in fortran.
    if (iargc() < 2) call usage()
    
    call vs_getarg( 1, basefn_in )
    nipx = 1
    nipz = 1
    if (iargc() >= 4) then 
        call int_getarg( 2, 1, huge(nipx), nipx )
        call int_getarg( 3, 1, huge(nipz), nipz )
    end if

    call set_database_in( basefn_in, nipx, nipz )
    
    if (iargc() > 4) then
        do ig=1, iargc()-4
            call int_getarg( ig+3, 0, huge(igt), igt)
            print *, ig, igt
            call set_mapping(ig,igt)
        end do
    end if

    call vs_getarg( iargc(), basefn_out )
    
    call set_database_out_nocreate( basefn_out )
    
    iline = 1
    line_loop : do
        call readline( copyentry, iostat, ok )
        if (iostat == IOSTAT_EOF) exit line_loop
        if (.not. ok) call die( "gfdb_redeploy: parsing of line "// iline &
                                 //" from standard input failed." )
        iline = iline+1
    end do line_loop

    call redeploy_cleanup()
    call cleanup()
    
end program
