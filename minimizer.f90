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

module minimizer_wrappers
  
  ! Comments starting with a double exclamation mark may be extracted and catenated
  ! by a simple script. I use Trac wiki syntax in these blocks, so that the wiki
  ! documentation on this program may be inlined here.
  
 !! = The `minimizer` tool =
  !
  ! 
  ! == Usage ==
  ! {{{
  ! > minimizer <<EOF
  ! [minimizer commands]
  ! ...
  ! EOF
  ! }}}
  ! 
  ! == Minimizer Commands ==
  ! 
    use minimizer_engine
    use better_varying_string
    use constants
    use util
    use unit
    use orthodrome
    use source
    use receiver
    use comparator
   ! use minimizer_engine
    
    implicit none

  ! wrappers doing command line processing
    public do_set_database
    public do_set_local_interpolation
    public do_set_spacial_undersampling
    public do_set_receivers
    public do_set_ref_seismograms
    public do_set_source_constraints
    public do_set_source_location
    public do_set_source_crustal_thickness_limit
    public do_get_source_crustal_thickness
    public do_set_source_params
    public do_set_source_params_mask
    public do_set_source_subparams
    public do_set_effective_dt
    public do_set_misfit_method
    public do_set_misfit_filter
    public do_set_misfit_taper
    public do_set_synthetics_factor
    public do_minimize_lm
    public do_output_seismograms
    public do_output_seismogram_spectra
    public do_output_source_model
    public do_get_source_subparams
    public do_get_global_misfit
    public do_get_misfits
    public do_get_peak_amplitudes
    public do_get_principal_axes
    public do_output_distances
    public do_output_cross_correlations
    public do_shift_ref_seismogram
    public do_autoshift_ref_seismogram
    public do_get_cached_traces_memory
    public do_set_cached_traces_memory_limit
    public do_set_verbose

  contains
  
    subroutine do_set_database( args, answer, ok )
        
     !! === {{{set_database dbpath [ nipx nipz ]}}} ===
      !
      ! Select Greens function database to use for the calculation of synthetic seismograms.
      !
      ! {{{dbpath}}} is the path to a Greens function database created with {{{gfdb_build}}}.
      ! This is the path without the filename extensions {{{.index}}} or {{{.chunk}}}.
      !
      ! {{{nipx}}} {{{nipz}}} turn on Gulunay's interpolation in the Greens function database
      ! if set to values other than one. A Greens function database opened this way will pretend to have {{{nipx}}} times the
      ! number of traces in the horizontal direction, inserting interpolated traces as needed.
      ! Same applies with {{{nipz}}} to the vertical.
      ! Gulunay's generalized FK interpolation is used to fill the interpolated traces.
      ! If either of {{{nipx}}} or {{{nipz}}} is set to one, a 2D interpolation (time-distance or time-depth)
      ! is performed. If both {{{nipx}}} and {{{nipz}}} are set to the same value,
      ! a 3D (time-distance-depth) interpolation is performed.
      ! If {{{nipx}}} and {{{nipz}}} are not set to the same value, first horizontal
      ! 2D interpolation is applied followed by a vertical 2D interpolation.
      !
      ! '''Note:''' Gulunay's interpolation works in the spectral domain and uses FFTs and 
      ! thus has cyclic properties.
      ! To prevent wrap-around artifacts, the interpolation is done block-wise with some overlap.
      ! At the boundaries of the database, repeating end points are used to gain a margin and
      ! enough traces for the interpolation. Nevertheless, this introduces errors 
      ! near to the surface and at the maximum depth, as well as at the ends of the 
      ! distance range of the database.
      
        type(varying_string), intent(in)  :: args
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        
        character(len=len(args)) :: buffer
        type(varying_string)  :: db_path
        integer :: nipx, nipz, iostat
        type(varying_string)   :: args_mut
        
        answer = ''
        ok = .false.
        
        nipx = 1
        nipz = 1
        if (3 == count_words( char(args) )) then
            args_mut = args
            call split(args_mut, db_path, ' ') 
            buffer = char(args_mut)
            read (unit=buffer,fmt=*,iostat=iostat) nipx, nipz
            if (iostat /= 0) then
                call error( "set_database: failed to parse arguments" )
                return
            end if
        else
            db_path = args
        end if
        
        if (nipx < 1 .or. nipz < 1) then
            call error( "set_database: nipx and nipz must be positive" )
            return
        end if
        
        call set_database( db_path, nipx, nipz, ok )
        
    end subroutine

    
    subroutine do_set_local_interpolation( arg, answer, ok )
    
      !! === {{{set_local_interpolation ( nearest_neighbor | bilinear )}}} ===
       ! 
       ! Set local interpolation method used during calculation of synthetic seismograms.
       !
    
        type(varying_string), intent(in)  :: arg
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
    
        ok = .false.
        answer = ''
        
        if (arg == 'nearest_neighbor') then
            call set_local_interpolation(.false.)
        else if (arg == 'bilinear') then
            call set_local_interpolation(.true.)
        else 
            call error( "set_local_interpolation: unknown interpolation method: "//arg)
            return
        end if
        ok = .true.
            
    end subroutine

    subroutine do_set_spacial_undersampling( line, answer, ok )

      !! === {{{set_spacial_undersampling nxunder nzunder }}} ===
       ! 
       ! Tell minimizer to use only a subset of the databases Green's functions.
       !
       !   nxunder: use every nxunder'th horizontal Green's function distance.
       !   nzunder: use every nzunder'th vertical Green's function depth. 

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok

        character(len=len(line)) :: buffer
        integer :: iostat, xunder, zunder
        
        answer = ''
        ok = .false.
        buffer = char(line)
        read (unit=buffer,fmt=*,iostat=iostat) xunder, zunder
        if (iostat /= 0) then
            call error( "set_spacial_undersampling: failed to parse arguments" )
            return
        end if
        call set_spacial_undersampling( xunder, zunder, ok )

    end subroutine
        
    subroutine do_set_receivers( receiversfn, answer, ok )
      
     !! === {{{set_receivers filename}}} ===
      !
      ! Read a list of receiver coordinates from three column (lat lon components) ascii file {{{filename}}}.
      !
      ! The file format is as follows:
      !
      !  * first column:   latitude in degrees
      !  * second column:  longitude in degrees
      !  * third column:   selected components of the station; for every component you want,
      !    add one of the characters below:
      !     * radial component:
      !        * a = positive is displacement away from source
      !        * c = positive is displacement coming towards source
      !     * transversal component:
      !        * r = positive is rightwards seen from source
      !        * l = positive is leftwards seen from source
      !     * vertical component:
      !        * d = positive is downwards
      !        * u = positive is upwards
      !     * horizontal component (north-south)
      !        * n = positive is north
      !        * s = positive is south
      !     * horizontal component (east-west)
      !        * e = positive is east
      !        * w = positive is west
      !
      ! Adding the same component more than once is not allowed, so at most 5 components may be given.
      ! Lines starting with a '#' are considered to be comment lines.
      !
      ! An example receivers file might look like this:
      ! {{{
      ! # ok:
      ! 42.35    13.4  ard
      ! 49.78   17.54  ard
      ! # north component broken:
      ! 45.49   25.95  ed
      ! # only vertical component available
      ! 47.92 19.89 d
      ! 35.87 14.52 d       
      ! }}}
      
        type(varying_string), intent(in)  :: receiversfn
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        call set_receivers( receiversfn, answer, ok )
        
    end subroutine

    subroutine do_switch_receiver( line, answer, ok )

     !! === {{{switch_receiver ireceiver ( on | off ) }}} ===
      !
      ! Turn receiver number ireceiver on or off.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        
        character(len=len(line)) :: buffer, onoff
        integer :: nerr, ireceiver
        logical :: state

        answer = ''
        ok = .true.

        buffer = char(line)
        read (unit=buffer,fmt=*, iostat=nerr) ireceiver, onoff
        if (nerr /= 0) then
            call error( "usage: switch_receiver ireceiver ( on | off )" )
            ok = .false.
            return
        end if
        
        state  = .true.
        if (onoff == 'on') then
            state = .true.
        else if (onoff == 'off') then
            state = .false.
        else 
            call error( "usage: switch_receiver ireceiver ( on | off )" )
            ok = .false.
            return
        end if
    
        call switch_receiver( ireceiver, state, ok )

    end subroutine
    
    subroutine do_set_ref_seismograms( line, answer, ok )
      
     !! === {{{set_ref_seismograms filenamebase format}}} ===
      !
      ! Read a set of reference seismograms.
      !
      ! For every component at every of the receivers which have 
      ! been set with {{{set_receivers}}} one file must be povided.
      !
      ! Currently the following formats are available:
      !
      !  * {{{table}}}: ASCII tables with two columns: time [s] and displacement [m].
      !  * {{{mseed}}}: Single trace "Data Only SEED Volume" 
      !    (Mini-SEED, http://www.iris.edu/manuals/SEED_appG.htm).
      !  * {{{sac}}}: [wiki:SacBinaryFile SAC binary file]. 
      !    Please note, that this file format is platform dependant.
      !
      ! The files are expected to be named using the following scheme:
      ! 
      !   {{{$filenamebase-$ReceiverNumber-$ComponentCharacter.$format}}}
      !
      ! where
      !
      !  * {{{$ReceiverNumber}}} is the number of the receiver, as defined by the ordering
      !    of receivers in the receiver file (see {{{set_receivers}}}).
      !  * {{{$ComponentCharacter}}} is one of the characters defining receiver components
      !    as described in {{{set_receivers}}}.
        
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        type(varying_string)  :: fformat, reffnbase
        
        answer = ''
        ok = .true.
        fformat = line
        call split(fformat, reffnbase, ' ')
        call set_ref_seismograms( reffnbase, fformat, ok )

    end subroutine

    subroutine do_shift_ref_seismogram( line, answer, ok )

     !! === {{{shift_ref_seismogram ireceiver shift}}} ===
      !
      ! Timeshift reference seismogram.
      !
      ! Shift reference seismogram at receiver number `ireceiver` by `shift` seconds.
      !

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok

        character(len=len(line)) :: buffer
        integer :: nerr, ireceiver
        real :: shift

        answer = ''
        ok = .true.

        buffer = char(line)
        read (unit=buffer,fmt=*, iostat=nerr) ireceiver, shift
        if (nerr /= 0) then
            call error( "usage: shift_ref_seismogram ireceiver shift" )
            ok = .false.
            return
        end if
        
        call shift_ref_seismogram( ireceiver, shift, ok )

    end subroutine

    subroutine do_autoshift_ref_seismogram( line, answer, ok )

     !! === {{{autoshift_ref_seismogram ireceiver min-shift max-shift}}} ===
      !
      ! Automatically timeshift reference seismogram.
      !
      ! Shift reference seismogram at receiver number `ireceiver` to where 
      ! the cross-correlation has a maximum in the interval [`min-shift`,`max-shift`].
      !
      ! If `ireceiver` is set to zero, all seismograms are auto-shifted. 

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok

        character(len=len(line)) :: buffer
        integer :: nerr, ireceiver
        real, dimension(2) :: shiftrange
        real, dimension(:), allocatable :: shifts
        
        answer = ''
        ok = .true.

        buffer = char(line)
        read (unit=buffer,fmt=*, iostat=nerr) ireceiver, shiftrange(1), shiftrange(2)
        if (nerr /= 0) then
            call error( "usage: autoshift_ref_seismogram ireceiver min-shift max-shift" )
            ok = .false.
            return
        end if
        call autoshift_ref_seismogram( ireceiver, shiftrange, shifts, ok )
        do ireceiver=1,size(shifts,1)
            answer = answer // " " // shifts(ireceiver)
        end do
        if (allocated(shifts)) deallocate(shifts)

    end subroutine
    
    subroutine  do_set_source_location( line, answer, ok )
      
     !! === {{{set_source_location latitude longitude reference-time}}} ===
      !
      ! Sets the source location and reference time.
      !
      !  * {{{latitude}}}, {{{longitude}}}: Geographical coordinates of source reference point in [degrees].
      !    All locations given in the source model description are measured relative 
      !    to this reference point.
      !  * {{{reference-time}}}: source reference time in seconds. 
      !    All times given in the source model description are measured relative 
      !    to this reference time.
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        integer :: nerr
        character(len=len(line)) :: buffer
        real :: lat, lon
        real(kind=8) :: ref_time
        
        answer = ''
        ok = .true.

        buffer = char(line)
        read (unit=buffer,fmt=*, iostat=nerr) lat, lon, ref_time
        if (nerr /= 0) then
            call error( "usage: set_source_location latitude longitude reference-time" )
            ok = .false.
            return
        end if
        call set_source_location(d2r(lat), d2r(lon), ref_time)

    end subroutine

    subroutine do_set_source_constraints( line, answer, ok )

     !! === {{{set_source_constraints px1 py1 pz1 nx1 ny1 nz1 ...}}} ===
      !
      ! Set constraining planes which affect source geometry for certain source models.
      ! 
      ! Each constraining plane is defined by a point and a normal vector.
      ! They are specified in the local carthesian coordinate system at the source, which has its principal
      ! axes pointing north, east, and downward, and whose origin is at the surface
      ! at the coordinates given with set_source_location.
      !
      !  * {{{px1 py1 pz1}}}: coordinates of point for plane number 1 in [m]
      !  * {{{nx1 ny1 nz1}}}: components of normal vector of plane number 1
      !
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok

        real, dimension(:,:), allocatable :: points, normals
        real, dimension(:), allocatable   :: numbers
        integer                           :: n, nplanes, iplane, iostat
        character(len=len(line))          :: buffer

        answer = ''
        ok = .true.
        buffer = char(line)
        n = count_words( buffer )
        
        if (mod(n,6) /= 0) then
            ok = .false.
            call error( "number of arguments is not divideable by 6" )
            return
        end if

        nplanes = n / 6
        allocate( numbers(n) )
        allocate( points(3,nplanes) )
        allocate( normals(3,nplanes) )

        read (unit=buffer,fmt=*,iostat=iostat) numbers(:)

        do iplane=1,nplanes
            points(:,iplane) = numbers((iplane-1)*6+1:(iplane-1)*6+3) 
            normals(:,iplane) = numbers((iplane-1)*6+4:(iplane-1)*6+6) 
            if (iostat > 0) then
                ok = .false.
                call error( "failed to parse numbers at plane " // iplane )
                return
            end if
        end do

        call set_source_constraints( points, normals )

        deallocate( points )
        deallocate( normals )
        deallocate( numbers )

    end subroutine

    subroutine do_set_source_crustal_thickness_limit( line, answer, ok )
      
     !! === {{{set_source_crustal_thickness_limit thickness-limit}}} ===
      !
      ! Limit crustal thickness at the source.
      !
      !  * {{{thickness-limit}}}: Maximal thickness of crust in [m].
      !
      ! Default values for the thickness are retrieved from the crust2x2 model.
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        integer :: nerr
        character(len=len(line)) :: buffer
        real :: thickness_limit
        
        answer = ''
        ok = .true.

        buffer = char(line)
        read (unit=buffer,fmt=*, iostat=nerr) thickness_limit
        if (nerr /= 0) then
            call error( "usage: set_source_crustal_thickness_limit thickness-limit" )
            ok = .false.
            return
        end if
        call set_source_crustal_thickness_limit(thickness_limit)

    end subroutine

    subroutine do_get_source_crustal_thickness( line, answer, ok )

     !! === {{{get_source_crustal_thickness}}} ===
      !
      ! Returns crustal thickness at the source in [m].
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok

        real :: thickness

        answer = ''
        ok = .false.
        
        call get_source_crustal_thickness( thickness, ok )
        
        if (ok) answer = thickness
        
    end subroutine
    
    subroutine do_set_source_params( line, answer, ok )
      
     !! === {{{set_source_params source-type source-params ...}}} ===
      !
      ! Sets the source type and parameters.
      !
      ! The available source types and a complete description of their parameters 
      ! are given in the [wiki:SourceTypes source type documentation].
      ! Short descriptions can be queried using the [wiki:SourceInfoTool source_info] tool.
      ! 
      ! This function detects if the same source parameters have already been 
      ! set, so that seismograms are not recalculated when the same source 
      ! is set several times.
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
    
        integer                             :: sourcetype, nparams, iostat
        real, dimension(:), allocatable     :: params
        type(varying_string)                :: sourcetypename, paramsstr
        character(len=len(line))            :: buffer
        
        answer = ''
        ok = .true.
   
        paramsstr = line
        call split( paramsstr, sourcetypename, " " )
        call psm_get_source_id( sourcetypename, sourcetype )
       
        if (sourcetype == 0) then
            call error("unknown source type name: " // sourcetypename)
            ok = .false.
            return
        end if
        
        nparams = psm_get_n_source_params( sourcetype )
        call resize( params, 1, nparams )
    
        buffer = char(paramsstr)

        if (count_words( buffer ) /= nparams) then
            call error("source of type '" // sourcetypename // "' requires "// nparams // " parameters.")
            ok = .false.
            return
        end if
        
        read (unit=buffer,fmt=*,iostat=iostat) params
        if (iostat > 0) then
            call error("failed to parse source params" )
            ok = .false.
            return
        end if
        
        call set_source_params( sourcetype, params, ok )
    
    end subroutine
    
    subroutine do_set_source_params_mask( line, answer, ok )
      
     !! === {{{set_source_params_mask mask ...}}} ===
      !
      ! Select inversion parameters for the next minimization with {{{minimize_lm}}}.
      !
      ! {{{mask}}} is built by giving a 'T' or 'F' for every source parameter of
      ! the source type that is currently in use.
      ! 'T' makes the corresponding parameter an actual inversion parameter, 
      ! 'F' fixes the corresponding parameter to its current value.
      ! The 'T's and 'F's must be separated by whitespace.
      !
      ! The values of the selected parameters can be set using {{{set_subparams}}}
      ! and queried using {{{get_subparams}}}.
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        logical, dimension(:), allocatable :: params_mask
        integer :: iostat
        character(len=len(line)) :: buffer
        
        integer :: nsubparams
        
        answer = ''
        ok = .true.
        buffer = char(line)
        nsubparams = count_words( buffer )
        allocate( params_mask(nsubparams) )
        
        read (unit=buffer,fmt=*,iostat=iostat) params_mask(:)
        if (iostat > 0) then
            ok = .false.
            call error( "failed to parse source parameter mask" )
            return
        end if
        
        call set_source_params_mask( params_mask, ok )
        
        deallocate( params_mask )
        
    end subroutine
    
    subroutine do_set_source_subparams( line, answer, ok )
    
     !! === {{{set_source_subparams subparams ...}}} ===
      !
      ! Assignes values to the currently selected inversion parameters.
      !
      ! This command expects one value for each parameter selected with {{{set_source_params_mask}}}.
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real, dimension(:), allocatable :: subparams
        
        character(len=len(line)) :: buffer
        integer :: iostat, nsubparams
        
        buffer = char(line)
        nsubparams  = count_words(buffer )
        allocate(subparams(nsubparams))
        ok = .true.
        answer = ''
        read (unit=buffer,fmt=*,iostat=iostat) subparams
        if (iostat > 0) then
            ok = .false.
            call error( "failed to parse source parameters" )
            return
        end if
        
        call set_source_subparams( subparams, ok )
        deallocate(subparams)
        
    end subroutine
    
    subroutine do_set_effective_dt( line, answer, ok )
      
     !! === {{{set_effective_dt effective_dt}}} ===
      !
      ! Sets the effective dt controlling the source parameterization.
      
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real :: effective_dt_
        integer :: iostat
        character(len=len(line)) :: buffer
        
        answer = ''
        ok = .true.
       
        buffer = char(line)
        read (unit=buffer,fmt=*,iostat=iostat) effective_dt_
        if (iostat > 0) then
            ok = .false.
            call error( "failed to parse effective dt" )
            return
        end if
        
        call set_effective_dt( effective_dt_ )
  
    end subroutine
    
    subroutine do_set_misfit_method( line, answer, ok )
      
     !! === {{{set_misfit_method ( l2norm | l1norm | ampspec_l2norm | ampspec_l1norm | scalar_product | peak )}}} ===
      !
      ! Set the misfit calculation method.
      !
      ! Available methods are:
      !  * {{{l2norm}}}: L2 norm is done on difference of time traces
      !  * {{{l1norm}}}: L1 norm is done on difference of time traces
      !  * {{{ampspec_l2norm}}}: L2 norm is done on difference of amplitude spectra
      !  * {{{ampspec_l1norm}}}: L2 norm is done on difference of amplitude spectra
      !  * {{{scalar_product}}}: instead of a norm, the scalar product is calculated
      !  * {{{peak}}}: instead of a norm, the peak amplitudes are returned
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
    
        integer :: id
        
        ok = .true.
        answer = ''
        call comparator_get_norm_id( line, id )
        if (id == 0) then
            ok = .false.
            call error( "unknown norm method: "//line )
        end if
        
        call set_misfit_method( id )      
    end subroutine
    
    subroutine do_set_misfit_filter( line, answer, ok )
      
     !! === {{{set_misfit_filter x0 y0  x1 y1  ...}}} ===
      !
      ! Defines a piecewise linear function which is multiplied to the
      ! spectra before calculating misfits in the frequency domain.
      !
      !  * {{{x0 y0  x1 y1 ...}}}: Control points with {{{xi}}}: frequency [Hz] and 
      !    {{{yi}}}: multiplicator amplitude.
      !
      ! The amplitude drops to zero before the first and after the last control point.
      !
      ! Example: {{{set_misfit_filter  0.2 1  0.5 1}}} defines a rectangular 
      ! window between 0.2 and 0.5 Hz.
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        integer :: i, n, nwords, iostat
        real, dimension(:), allocatable :: x,y
        character(len=len(line)) :: buffer
        
        answer = ''
        ok = .true.
        
        buffer = char(line)
        nwords = count_words( buffer )
        
        n = nwords/2
        allocate(x(n),y(n))
        
        buffer = char(line)
        read (unit=buffer,fmt=*,iostat=iostat) (x(i), y(i), i=1,n)
        if (iostat > 0) then
            ok = .false.
            call error( "failed to parse coordinates" )
            return
        end if
        
        call set_misfit_filter( x, y )
        
        deallocate(x,y)
                
    end subroutine
    
    subroutine do_set_misfit_taper( line, answer, ok )
      
     !! === {{{set_misfit_taper ireceiver x0 y0  x1 y1 ...}}} ===
      !
      ! Defines a piecewise linear function which is multiplied to seismogram
      ! traces before calculating spectra or misfits.
      !
      !  * {{{ireceiver}}}: Number of the receiver to which the taper shall be applied.
      !  * {{{x0 y0  x1 y1 ...}}}: Control points with {{{xi}}}: time [s] and 
      !    {{{yi}}}: multiplicator amplitude.
      !
      ! The amplitude drops to zero before the first and after the last control point.
      !
      ! Example: {{{set_misfit_taper  120 1  150 1}}} defines a rectangular 
      ! window between 120 and 150 s
        
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        integer :: i, n, nwords, iostat
        real, dimension(:), allocatable :: x,y
        character(len=len(line)) :: buffer
        integer :: ireceiver
        
        answer = ''
        ok = .true.
        
        buffer = char(line)
        nwords = count_words( buffer )
        
        n = (nwords-1)/2
        allocate(x(n),y(n))
        
        buffer = char(line)
        read (unit=buffer,fmt=*,iostat=iostat) ireceiver, (x(i), y(i), i=1,n)
        if (iostat > 0) then
            ok = .false.
            call error( "failed to parse values" )
            return
        end if
        
        call set_misfit_taper( ireceiver, x, y, ok )
        
        deallocate(x,y)
                
    end subroutine

    subroutine do_set_synthetics_factor( line, answer, ok )
      
     !! === {{{set_synthecics_factor factor}}} ===
      !
      ! Scale amplitude of synthecic seismograms by this factor during 
      ! misfit calculation.       
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real :: factor
        integer :: iostat
        character(len=len(line)) :: buffer
        
        answer = ''
        ok = .true.
       
        buffer = char(line)
        read (unit=buffer,fmt=*,iostat=iostat) factor
        if (iostat > 0) then
            ok = .false.
            call error( "failed to parse synthetics factor" )
            return
        end if
        
        call set_synthetics_factor( factor )
  
    end subroutine

    subroutine do_minimize_lm( line, answer, ok )
      
     !! === {{{minimize_lm}}} ===
      !
      ! Runs Levenberg-Marquardt minimization.
      !
      ! This tries to invert for the source parameters selected with {{{set_source_params_mask}}} 
      ! by searching for a minimum in the currently 
      ! selected misfit function, starting from the current source model parameterization.
      !
      ! This function makes use of {{{lmdif()}}} from MINPACK from Netlib.
      !
      ! '''Returns:''' {{{ info iterations misfit }}}
      ! 
      !  * {{{info}}}: Information on convergence, returned by lmdif(). See MINPACK documentation for details.
      !  * {{{iterations}}}: Number of function evaluations = number of source models tested.
      !  * {{{misfit}}}: Final global misfit.
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        integer :: info, iterations
        real :: misfit_
        
        ok = line /= '' ! get rid of warning, that line is not used
        
        call minimize_lm( info, iterations, misfit_, ok)
        if ( .not. ok ) then
            answer = '' 
            return
        end if
        
        answer = info // " " // iterations // " " // misfit_
        
    end subroutine
   
    subroutine do_output_source_model( line, answer, ok )
      
     !! === {{{output_source_model filenamebase}}} ===
      !
      ! Output information about the current source model.
          
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
    
        answer = ''
        ok = .true.
        
        call output_source_model(line, ok) ! dies on failure
        
    end subroutine
    
    subroutine do_output_seismogram_spectra( line, answer, ok )
        
     !! === {{{output_seismogram_spectra filenamebase (synthetics|references) (plain|filtered)}}} ===
      !
      ! Output the seismogram spectra which are used during misfit calculation.
      !
      !   * If the first argument is {{{references}}}, the spectra of the
      !     reference seismograms are outputted, if it is {{{synthetics}}},
      !     those of the synthetic seismograms for the current source model 
      !     are outputted.
      !   * If a filter has been set using {{{set_misfit_filter}}}, the 
      !     filtered spectra are written.
      !   * For every selected component at every defined receiver one
      !     file is written.
      ! 
      ! The files are named using the following scheme:
      !
      !     {{{$filenamebase-$ReceiverNumber-$ComponentCharacter.table}}}
      !
      ! where
      !
      !   * {{{$ReceiverNumber}}} is the number of the receiver, as defined by the ordering
      !     of receivers in the receiver file (see {{{set_receivers}}}).
      !   * {{{$ComponentCharacter}}} is one of the characters defining receiver components
      !     as described in {{{set_receivers}}}.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        type(varying_string) :: filenamebase, str_probe, str
        integer :: which_probe, which_processing
        
        answer = ''
        ok = .true.
        str = line
        call split( str, filenamebase, ' ' )
        call split( str, str_probe, ' ' )

        which_probe = SYNTHETICS
        if (str_probe .eq. 'references') which_probe = REFERENCES

        which_processing = PLAIN
        if (str .eq. 'filtered') which_processing = FILTERED

        call output_seismogram_spectra( filenamebase, which_probe, which_processing, ok )
            
    end subroutine
    
    subroutine do_output_seismograms( line, answer, ok )
    
     !! === {{{output_seismograms filenamebase fileformat (synthetics|references) (plain|tapered|filtered) }}} ===
      !
      ! Output current synthetic or reference seismograms.
      ! 
      !  * {{{filenamebase}}}: Stem for the creation of filenames, see below for details.
      !  * {{{fileformat}}}: Format of the ouputted files.
      !  * {{{[ tapered ]}}}: If this argument is present and a taper has been set using 
      !    {{{set_misfit_taper}}}, the tapered seismograms are written.
      !
      ! The files are named using the following scheme:
      !   
      !   {{{$filenamebase-$ReceiverNumber-$ComponentCharacter.$fileformat}}}
      !
      ! where
      !
      !  * {{{$ReceiverNumber}}} is the number of the receiver, as defined by the ordering
      !    of receivers in the receiver file (see {{{set_receivers}}}).
      !  * {{{$ComponentCharacter}}} is one of the characters defining receiver components
      !    as described in {{{set_receivers}}}.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        type(varying_string) :: fileformat, filenamebase, str_probe, str
        integer :: which_probe, which_processing

        answer = ''
        ok = .true.
        
        str = line
        call split( str, filenamebase, ' ' )
        call split( str, fileformat, ' ' )
        call split( str, str_probe, ' ' )

        which_probe = SYNTHETICS
        if (str_probe .eq. 'references') which_probe = REFERENCES

        which_processing = PLAIN
        if (str .eq. 'tapered') which_processing = TAPERED
        if (str .eq. 'filtered') which_processing = FILTERED

        call output_seismograms( filenamebase, fileformat, which_probe, which_processing, ok )
        
    end subroutine

    subroutine do_get_source_subparams( line, answer, ok )
    
     !! === {{{get_source_subparams}}} ===
      ! 
      ! Returns the current values of the source parameters selected with 
      ! {{{set_source_params_mask}}}.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real, dimension(:), allocatable :: subparams
        integer :: i
        
        ok = line /= '' ! get rid of warning, that line is not used
    
        answer = ''
        ok = .true.

        call get_source_subparams( subparams, ok )
        do i=1,size(subparams)
            answer = answer // subparams(i)
            if (i<size(subparams)) answer = answer // " "
        end do
    
    end subroutine
    
   
    subroutine do_get_global_misfit( line, answer, ok )
    
     !! === {{{get_global_misfit}}} ===
      ! 
      ! Returns the global misfit between the synthetic seismograms for the current 
      ! source model and the reference seismograms.
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real :: misfit_
        
        ok = line /= '' ! get rid of warning, that line is not used
        answer = ""
        call get_global_misfit( misfit_, ok )
        if (.not. ok) return
     
        answer =  misfit_
        
    end subroutine
    
    subroutine array2d_to_string( array, string )
    
        real, dimension(:,:), intent(in) :: array
        type(varying_string), intent(out) :: string
        
        character(len=size(array)*32) :: buffer
        integer :: ix,iy
        
        write (buffer,*) ((array(ix,iy), ix=1,size(array,1)), iy=1,size(array,2))
                
        string = var_str(trim(buffer))
    
    end subroutine
    
    subroutine do_get_misfits( line, answer, ok )
    
     !! === {{{get_misfits}}} ===
      ! 
      ! Returns the misfit and normalization factors between the synthetic 
      ! seismograms for the current source model and the reference seismograms.
      !
      ! Disabled stations are omitted in output list.
      !
      ! Returns: {{{misfit-receiver1-component1 normfactor-receiver1-component1 ...}}}
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real, dimension(:,:), allocatable :: misfits_
        
        ok = line /= '' ! get rid of warning, that line is not used
        answer = ""
        
        call get_misfits( misfits_, ok )
        if (.not. ok) return
     
        call array2d_to_string( misfits_, answer )
        if (allocated( misfits_)) deallocate(misfits_)
        
    end subroutine
    
    subroutine do_get_peak_amplitudes( line, answer, ok )
    
     !! === {{{get_peak_amplitudes}}} ===
      ! 
      ! Get the horizonal and vertical peak amplitudes of the synthetic traces. 
      !
      ! Disabled stations are omitted in output list.
      !
      ! Returns: {{{maxabs_receiver_1_horizontal maxabs_receiver_1_vertical ...}}}
      
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real, dimension(:,:), allocatable :: maxabs_
        
        ok = line /= '' ! get rid of warning, that line is not used
        answer = ""
        
        call get_peak_amplitudes( maxabs_, ok )
        if (.not. ok) return
     
        call array2d_to_string( maxabs_, answer )
        if (allocated(maxabs_)) deallocate(maxabs_)
        
    end subroutine

    subroutine do_get_principal_axes( line, answer, ok )
    
     !! === {{{get_principal_axes}}} ===
      ! 
      ! Get the orientation of the principal axes P and T of the current source model.
      ! 
      ! Returns: {{{p-axis-phi p-axis-theta t-axis-phi t-axis-theta}}}
      !
      !  * {{{p-axis-phi, p-axis-theta}}}: Spherical coordinates of the direction of P.
      !  * {{{t-axis-phi, t-axis-theta}}}: Spherical coordinates of the direction of T.
      !
      ! These are ordinary spherical coordinates based on the 
      ! local carthesian north-east-down coordinate system at the source.
       
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        real, dimension(2)    :: pax, tax
        
        ok = line /= '' ! get rid of warning, that line is not used
        answer = ''
        
        call get_principal_axes( pax, tax, ok )
        if (.not. ok) return
       
        answer = pax(1) // " " // pax(2) // " " // tax(1)  // " " // tax(2)
    
    end subroutine
    
    subroutine do_output_distances( line, answer, ok )
    
     !! === {{{output_distances filename}}} ===
      ! 
      ! Dump epicentral distances and azumiths to ascii file {{{filename}}}.
   
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        integer :: idist, unit, iostat
        real(kind=8), dimension(:), allocatable :: distances, azimuths
        
        ok = .true.
        answer = ''

        call get_distances( distances, azimuths, ok )
        if (.not. ok) return
        
        call claim_unit( unit )
        open( unit=unit, file=char(line), status='unknown', iostat=iostat )
        if (iostat > 0) then
            call error("failed to open file for output: "//line)
            ok = .false.
        end if
        
        do idist=1,size(distances)
            write (unit,*) r2d(distances(idist)), distances(idist)*earthradius, r2d(azimuths(idist))
        end do
        
        close( unit ) 
        call release_unit( unit )
    
        call resize( distances, 1, 0 )
        call resize( azimuths, 1, 0 )
        
    end subroutine
       
    subroutine do_output_cross_correlations( line, answer, ok )

     !! === {{{output_cross_correlations filenamebase shift-min shift-max}}} ===
      !
      ! Output cross-correlations between synthetics and references
      !
      !  * {{{filenamebase}}}: Stem for the creation of filenames, see below for details.
      !  * {{{shift-min shift-max}}}: Range of shifts for which cross-correlations are evaluated. (in [s]).
      !
      ! The files are named using the following scheme:
      !
      !   {{{$filenamebase-$ReceiverNumber-$ComponentCharacter.$fileformat}}}
      !
      ! where
      !
      !  * {{{$ReceiverNumber}}} is the number of the receiver, as defined by the ordering
      !    of receivers in the receiver file (see {{{set_receivers}}}).
      !  * {{{$ComponentCharacter}}} is one of the characters defining receiver components
      !    as described in {{{set_receivers}}}.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        type(varying_string) :: str_shiftmin, filenamebase, str
        real, dimension(2) :: r_shift

        answer = ''
        ok = .true.
        
        str = line
        call split( str, filenamebase, ' ' )
        call split( str, str_shiftmin, ' ' )

        r_shift(1) = str_shiftmin
        r_shift(2) = str


        call output_cross_correlations( filenamebase, r_shift, ok )

    end subroutine

    subroutine do_get_cached_traces_memory( line, answer, ok )
     
     !! === {{{get_cached_traces_memory}}} ===
      !
      ! Get memory usage by Green's function database cache
      ! 
      ! Returns number of bytes allocated for traces in the Greens function database cache.
      ! This number does not contain the overhead of header data in the traces, and index tables.
      ! It is the plain number of bytes used to hold the seismogram traces.
     
        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out) :: ok 

        integer(kind=8) :: nbytes
        character(len=128) :: buffer

        answer = ''
        ok = line /= '' ! get rid of warning, that line is not used

        call get_cached_traces_memory(nbytes, ok)
        write (buffer,*) nbytes
        answer = buffer

    end subroutine
    
    subroutine do_set_cached_traces_memory_limit( line, answer,  ok )

     !! === {{{set_cached_traces_memory_limit nbytes}}} ===
      !
      ! Set maximum memory usage by Green's function database cache
      ! 
      ! Sets the approximate maximum of memory used by the Greens function database cache.
      ! This limit does not include the overhead of header data in the traces, and index tables.
      ! It is the plain number of bytes which the Green's function database is allowed to use
      ! to cache seismogram traces.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out) :: ok 

        integer(kind=8)          :: nbytes_limit
        character(len=len(line)) :: buffer
        integer                  :: nerr

        answer = ''
        ok = .true.

        buffer = char(line)
        read (unit=buffer,fmt=*, iostat=nerr) nbytes_limit
        if (nerr /= 0) then
            call error( "usage: set_cached_traces_memory_limit limit" )
            ok = .false.
            return
        end if
        
        call set_cached_traces_memory_limit( nbytes_limit, ok )

    end subroutine

    subroutine do_set_verbose( line, answer, ok )

     !! === {{{set_verbose (T|F)}}} ===
      !
      ! Toggle verbose operation.
      !
      ! T turns on verbose operation.
      ! F turns it off.

        type(varying_string), intent(in)  :: line
        type(varying_string), intent(out) :: answer
        logical, intent(out)              :: ok
        
        answer = ''
        ok = .true.
        
        if (line == 'T') then
            g_verbose = .true.
        else if (line == 'F') then
            g_verbose = .false.
        else 
            call error("usage: set_verbose (T|F)")
            ok = .false.
        end if

    end subroutine
    
   !! == Example ==
    ! 
    ! To calculate synthetic seismograms at 11 receivers for the Izmit event:
    ! 
    ! {{{
    ! # setup receivers, indicating that for each receiver 
    ! # north, east and down components should be calculated.
    ! > cat >izmit-receivers.table <<EOF
    ! 42.350 13.400 ned
    ! 49.780 17.540 ned
    ! 45.490 25.950 ned
    ! 47.920 19.890 ned
    ! 35.870 14.520 ned
    ! 34.960 33.330 ned
    ! 35.280 24.890 ned
    ! 35.180 25.500 ned
    ! 49.630 22.710 ned
    ! 36.370 25.460 ned
    ! 42.620 23.240 ned
    ! EOF
    ! 
    ! > minimizer <<EOF
    ! set_database            /gfdb/gemini-prem/db
    ! set_effective_dt        0.5
    ! set_receivers           izmit-receivers.table
    ! set_source_location     40.75 29.86 0
    ! set_source_params       bilateral 0 0 0 10000 2e20  91 87 164  0  40000 20000 18000  3500 2
    ! calculate_seismograms
    ! output_seismograms      izmit-seismogram  mseed
    ! EOF
    ! }}}
    !

end module

program minimizer
  
    use constants
    use orthodrome
    use util
    use unit
    use source
    use better_varying_string
    use varying_string_getarg
    use minimizer_wrappers
    use crust2x2
    
    implicit none
        
    type(varying_string)           :: line, command, answer
    integer                        :: iostat
    character, parameter           :: eol = char(10)
    logical                        :: ok
    type(varying_string)           :: aux_path
    
    g_pn = "minimizer"
    
    g_usage = "usage: "//g_pn// eol // &
              "..." // eol // &
              "EOF" // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/MinimizerTool"
    
    call vs_getenv( 'KIWI_HOME', aux_path )
    if (aux_path == '') then
        call die('Environment variable KIWI_HOME not set.' //eol // eol // &
                 'KIWI_HOME should contain the path to the directory' // eol // &
                 'containing the "aux" directory where minimizer looks for crustal' // eol // &
                 'models etc.' // eol // &
                 'An initial "aux" directory resides in the directory where' // eol // &
                 'you have built the kiwi tools, so you can set' // eol // &
                 'KIWI_HOME to your build directory.')
    end if
    aux_path = aux_path // '/aux'
    
    call crust2x2_load( char(aux_path // '/crust2x2'), ok )
    if (.not. ok) then
        call die()
    end if
    
    stdinloop: do  ! until end-of-file on stdin
        
        call get(line,iostat=iostat)
        if (iostat > 0 .or. iostat == -1) exit stdinloop
        
        call do_command( line, answer, command, ok )
        if (ok) then
            if (answer == '') then
                call put_line(command // ": ok")
            else
                call put_line(command // ": ok >")
                call put_line(answer)
            end if
        else
            if (g_errstr == '') then
                call put_line(command // ": nok")
            else
                call put_line(command // ": nok >")
                call put_line(g_errstr)
                call error("")
            end if
        end if

        call flush( stdout )
        
    end do stdinloop

    call cleanup_minimizer()
    
    call cleanup()
    
  contains
  
    subroutine do_command( line, answer, command, ok )
    
        type(varying_string), intent(in) :: line
        type(varying_string), intent(out) :: answer
        type(varying_string), intent(out) :: command
        
        logical, intent(out) :: ok
        type(varying_string) :: arguments
        
        call reduce_whitespace(line,arguments)
        if ( char(arguments) == '' ) then
            ok = .true.
            return
        end if
        
        call split( arguments, command, " " )
        ok = .false.
        if (command == 'set_database') then
            call do_set_database( arguments, answer, ok )
        else if (command == 'set_local_interpolation') then
            call do_set_local_interpolation( arguments, answer, ok )
        else if (command == 'set_spacial_undersampling') then
            call do_set_spacial_undersampling( arguments, answer, ok )
        else if (command == 'set_receivers') then
            call do_set_receivers( arguments, answer, ok )
        else if (command == 'switch_receiver') then
            call do_switch_receiver( arguments, answer, ok )
        else if (command == 'set_ref_seismograms') then
            call do_set_ref_seismograms( arguments, answer, ok )
        else if (command == 'shift_ref_seismogram') then
            call do_shift_ref_seismogram( arguments, answer, ok )
        else if (command == 'autoshift_ref_seismogram') then
            call do_autoshift_ref_seismogram( arguments, answer, ok )
        else if (command == 'set_source_location') then
            call do_set_source_location( arguments, answer, ok )
        else if (command == 'set_source_constraints') then
            call do_set_source_constraints( arguments, answer, ok )
        else if (command == 'set_source_crustal_thickness_limit') then
            call do_set_source_crustal_thickness_limit( arguments, answer, ok )
        else if (command == 'get_source_crustal_thickness') then
            call do_get_source_crustal_thickness( arguments, answer, ok )
        else if (command == 'set_source_params') then
            call do_set_source_params( arguments, answer, ok )
        else if (command == 'set_source_params_mask') then
            call do_set_source_params_mask( arguments, answer, ok )
        else if (command == 'set_source_subparams') then
            call do_set_source_subparams( arguments, answer, ok )
        else if (command == 'set_effective_dt') then
            call do_set_effective_dt( arguments, answer, ok )
        else if (command == 'minimize_lm') then
            call do_minimize_lm( arguments, answer, ok )
        else if (command == 'output_source_model') then
            call do_output_source_model( arguments, answer, ok )
        else if (command == 'output_seismograms') then
            call do_output_seismograms( arguments, answer, ok )
        else if (command == 'output_seismogram_spectra') then
            call do_output_seismogram_spectra( arguments, answer, ok )
        else if (command == 'get_source_subparams') then
            call do_get_source_subparams( arguments, answer, ok )
        else if (command == 'get_global_misfit') then
            call do_get_global_misfit( arguments, answer, ok )
        else if (command == 'get_misfits') then
            call do_get_misfits( arguments, answer, ok )
        else if (command == 'get_peak_amplitudes') then
            call do_get_peak_amplitudes( arguments, answer, ok )
        else if (command == 'get_principal_axes') then
            call do_get_principal_axes( arguments, answer, ok )
        else if (command == 'output_distances') then
            call do_output_distances( arguments, answer, ok )
        else if (command == 'set_misfit_filter') then
            call do_set_misfit_filter( arguments, answer, ok )
        else if (command == 'set_misfit_taper') then
            call do_set_misfit_taper( arguments, answer, ok )
        else if (command == 'set_synthetics_factor') then
            call do_set_synthetics_factor( arguments, answer, ok )
        else if (command == 'set_misfit_method') then
            call do_set_misfit_method( arguments, answer, ok )
        else if (command == 'output_cross_correlations') then
            call do_output_cross_correlations( arguments, answer, ok )
        else if (command == 'get_cached_traces_memory') then
            call do_get_cached_traces_memory( arguments, answer, ok )
        else if (command == 'set_cached_traces_memory_limit') then
            call do_set_cached_traces_memory_limit( arguments, answer, ok )
        else if (command == 'set_verbose') then
            call do_set_verbose( arguments, answer, ok )
        else
            call error("unknown command: "//command)
        end if
        
    end subroutine
    
    subroutine reduce_whitespace(in, out)
    
        type(varying_string), intent(in) :: in
        type(varying_string), intent(out) :: out
        
        character(len=len(in)) :: buff
        logical :: ws
        integer :: j,i
        
        j = 1
        ws = .true.
        do i=1,len(in)
            if (at(in,i) .ne. ' ') then
                if (at(in,i) .eq. '#') exit
                buff(j:j) = at(in,i)
                j = j+1
                ws = .false.
            else
                if (.not. ws) then
                    buff(j:j) = ' '
                    j = j+1
                    ws = .true.
                end if
            end if
        end do
        do i=j,len(in)
            buff(i:i) = ' '
        end do
        
        out = trim(buff)
    end subroutine
    
end program
