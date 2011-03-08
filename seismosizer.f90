! $Id: seismosizer.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


program seismosizer
 
  ! This program calculates seismograms at a fixed set of receivers for possibly
  !  many different source parameterizations at a common origin.
  !
  ! usage: seismosizer database effective-dt origin-lat origin-lon receivers \
  !          output-base output-format info-base [ reference-base reference-format ] <<EOF
  ! source-type source-params ...
  ! ...
  ! EOF
  !
  ! Complete documentation is available on
  ! 
  !   http://kinherd.org/power/trac/wiki/SeismosizerTool
  !

  
    use constants
    use util
    use orthodrome
    use unit
    use gfdb
    use source
    use seismogram
    use seismogram_io
    use sparse_trace
    use receiver
    use better_varying_string
    use varying_string_getarg
    use read_table
    use comparator
    
    implicit none
    
    integer :: nerr
    type(t_psm)                 :: psm
    type(t_tdsm)                :: tdsm
    type(varying_string)        :: dbfnbase, receiversfn, outfnbase, infofnbase, outformat, refformat
    type(varying_string)        :: string, outfn, infofn
    type(varying_string)        :: reffnbase, reffn, sourcetypename
    type(t_gfdb)                 :: db
    type(t_strip),dimension(3)  :: displacement
    integer, dimension(2)       :: span
    integer                     :: iargc, iostat, ofile
    integer                     :: ireceiver, nreceivers
    integer                     :: icomponent, icentroid
    type(t_receiver), dimension(:), allocatable     :: rec
    real, dimension(:,:), allocatable               :: receivers_tab
    logical                                         :: reference_provided
    real                                            :: deltat, toffset
    type(t_strip)                                   :: strip
    type(t_probe), dimension(:,:), allocatable      :: ref_probes, syn_probes
    real, dimension(:), allocatable                 :: temp_seismogram
    real, dimension(:,:), allocatable               :: dist
    real, dimension(:), allocatable                 :: params
    real                                            :: effective_dt
    type(t_geo_coords)                              :: origin
    double precision                                :: ref_time
    integer                                         :: sourcetype, nparams
    character, parameter :: eol = char(10)
    
    g_pn = "seismosizer"
    
    g_usage = "usage: "//g_pn//" database effective-dt origin-lat origin-lon receivers \" // eol // &
              "            output-base output-format info-base [ reference-base reference-format ] <<EOF" // eol // &
              "source-type source-params ..." // eol // &
              "..." // eol // &
              "EOF" // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/SeismosizerTool"
    
    psm%ref_time = 0.
    
    psm%origin = t_geo_coords(0., 0.)
    
    if (iargc() /= 8 .and. iargc() /= 10) then
        call usage()
    end if
    
    call vs_getarg( 1, dbfnbase )
    call gfdb_init(db, dbfnbase)
    
    call vs_getarg( 2, string )
    effective_dt = string
    
    call vs_getarg( 3, string )
    origin%lat = string
    
    call vs_getarg( 4, string )
    origin%lon = string
    
    
    ref_time = 0.0
    call psm_set_origin_and_time(psm, d2r(origin), ref_time)

    call vs_getarg( 5, receiversfn )
    call readtable_file( receivers_tab, receiversfn, min_rows=2, min_cols=1 )

    call vs_getarg( 6, outfnbase )
    call vs_getarg( 7, outformat )
    call vs_getarg( 8, infofnbase )
    
    reference_provided = .false.
    if (iargc() == 10) then
        call vs_getarg( 9, reffnbase )
        call vs_getarg( 10, refformat )
        reference_provided = .true.
    end if
    
    
    nreceivers = size(receivers_tab,2)
    allocate( rec(nreceivers) )
    if (reference_provided) then
        allocate( ref_probes(3,nreceivers) )
        allocate( syn_probes(3,nreceivers) )
        allocate( dist(3,nreceivers) )
    end if

    do ireceiver=1,nreceivers
        rec(ireceiver)%origin%lat = d2r(receivers_tab(1,ireceiver))
        rec(ireceiver)%origin%lon = d2r(receivers_tab(2,ireceiver))
        if (reference_provided) then
            do icomponent=1,3
                reffn = reffnbase // "-" // ireceiver // "-" // icomponent // "." // refformat
                call readseismogram( char(reffn), "*", temp_seismogram, toffset, &
                                     deltat, nerr )
                if (nerr /= 0) call die()
                if (abs(deltat - db%dt)>db%dt/10000.) then
                    call die("sampling distance of seismogram from file "// reffn // &
                             " does not match gfdb sampling distance" )
                end if
                call seismogram_to_strip( temp_seismogram, toffset, deltat, &
                                          strip )
                call probe_init( ref_probes(icomponent,ireceiver), db%dt )
                call probe_init( syn_probes(icomponent,ireceiver), db%dt )
                call probe_set_array( ref_probes(icomponent,ireceiver), strip )
            end do
        end if
    end do
    if ( allocated(temp_seismogram) ) deallocate(temp_seismogram)
    
    stdinloop: do  ! until end-of-file on stdin
    
      ! det number of source params for given source model type
      ! adjust size of parameter vector
        call get(sourcetypename," " ,iostat=iostat)
        if (iostat /= 0) exit stdinloop
        call psm_get_source_id( sourcetypename, sourcetype )

        if (sourcetype == 0) then
            call warn("unknown source type name: " // sourcetypename)
            read (*,iostat=iostat)
            if (iostat /= 0) exit stdinloop
            print *, "fail"
            cycle stdinloop
        end if
        
        nparams = psm_get_n_source_params( sourcetype )
        call resize( params, 1, nparams )
    
      ! read parameters for parameterized source model
        read (unit=*,fmt=*,iostat=iostat) params
        call psm_set(psm, sourcetype, params )
        
        if (iostat /= 0) exit stdinloop
        
      ! convert to discrete source
        call psm_to_tdsm( psm, tdsm, effective_dt )
      
      ! write some debugging information  
        if (infofnbase /= "OFF") then
            infofn = infofnbase // "-psm.info"
            call psm_write_info_file( psm, infofn )
            
            infofn = infofnbase // "-tdsm.info"
            call tdsm_write_info_file( tdsm, infofn )
        
          ! dump discrete source centroids to file
            call claim_unit( ofile )
            infofn = infofnbase // "-dsm.table"
            open( unit=ofile, file=char(infofn), status='unknown', iostat=iostat )
            if (iostat /= 0) call die( "failed to open output file: " // infofn )
            do icentroid=1,tdsm%centroids%n
                write (unit=ofile, fmt=*) &
                    tdsm%centroids%north(icentroid), tdsm%centroids%east(icentroid), &
                    tdsm%centroids%depth(icentroid), tdsm%centroids%time(icentroid)
            end do
            close( ofile )
            call release_unit( ofile )
        end if
           
      ! make seismograms at every receiver
        do ireceiver=1,nreceivers
            call make_seismogram( tdsm, rec(ireceiver), db, displacement )
            do icomponent=1,3
                span = strip_span(displacement(icomponent))
                if (outfnbase /= "OFF") then
                    outfn = outfnbase // "-" // ireceiver // "-" // icomponent // "." // outformat
                    call writeseismogram( char(outfn), "*", &
                                          displacement(icomponent)%data, &
                                          '', '', '', '', &
                                          (span(1)-1)*db%dt, db%dt, nerr )
                    if (nerr /= 0) call die( "failed to write output file: " // outfn )
                end if
                if (reference_provided) then
                    call probe_set_array( syn_probes(icomponent,ireceiver), displacement(icomponent) )
                    dist(icomponent,ireceiver) = &
                       probes_norm( ref_probes(icomponent,ireceiver), syn_probes(icomponent,ireceiver), L2NORM )
                end if
            end do
        end do
        
      ! indicate that we have finished this step
        if (reference_provided) then
            print *, sqrt(sum(dist*dist))*db%dt
        else
            print *, "ok"
        end if
        call flush( stdout )
        
    end do stdinloop
    
    call gfdb_destroy(db)
    
    call psm_cleanup()
    call cleanup()

  contains
  
    subroutine seismogram_to_strip( seismogram, tbegin, deltat, strip )
    
        real, dimension(:), intent(in) :: seismogram
        real, intent(in) :: tbegin, deltat
        type(t_strip), intent(inout) :: strip
        
        integer :: ibeg,nlen
        
        ibeg = nint(tbegin/deltat)
        nlen = size(seismogram,1)
        if (abs(ibeg*deltat - tbegin) > deltat/100.) then
            print *, ibeg*deltat, tbegin, abs(ibeg*deltat - tbegin)
            call die( "time of first sample of seismogram not "// &
                       "divideable by sampling distance" )
        end if
        
        call strip_init( (/ibeg+1,ibeg+nlen/), seismogram, strip )
    
    end subroutine
    
end program
