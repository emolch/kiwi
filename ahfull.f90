program ahfull

    ! This program calculates homogeneous fullspace seismograms for a 
    ! superposition of sources at a number of receivers
    
    ! usage: ahfull sources receivers material source-time-function dt output-filename-base output-format 
    
    ! sources:
    !   t x y z mxx myy mzz mxy mxz myz fx fy fz

    ! receivers
    !   x y z nf ff
    
    ! material
    !  rho alpha beta
        
    ! stf (but time column is ignored; dt from command line is rules!)
    !   t f

    ! output-format: 'table' or 'sac' or 'mseed'
    !    'table' : one ascii table file per receiver, component, and window
    !              filenames are: $fnbase-$receiver-$component-$window.table
    !              columns are: <time> <displacement> 
    !    'sac'   : one sac seismogram file per receiver, component, and window
    !              filenames are: $fnbase-$receiver-$component-$window.sac
    !    'mseed' : one mseed seismogram file per receiver, component, and window
    !              filenames are: $fnbase-$receiver-$component-$window.sac
    
    use util
    use better_varying_string
    use varying_string_getarg
    use read_table
    use seismogram_io
    use elseis_oo
    use unit
    ! use f90_unix_env
    
    implicit none
    
    type source_t
        real                              :: tshift
        real, dimension(3)                :: location
        real, dimension(3,3)              :: mt
        real, dimension(3)                :: sf
    end type source_t
       
    type receiver_t
        real, dimension(3)                :: location
        logical                           :: nfflag
        logical                           :: ffflag
    end type receiver_t
   
    real, allocatable, dimension(:,:)           :: sources_tab, receivers_tab, material_tab, stf_tab
    real, allocatable, dimension(:)             :: stf
    real                                        :: dt
    integer                                     :: isource, ireceiver, nsources, nreceivers, nstfsamples
    integer                                     :: iwindow, nwindows
    integer                                     :: nsamples, n,p,q
    type(elseis_t)                              :: es
    type(source_t), allocatable, dimension(:)   :: sources(:)
    type(receiver_t), allocatable, dimension(:) :: receivers(:)
    real                                        :: rho, alpha, beta
    real, dimension(:,:), allocatable           :: seismograms
    real, dimension(3)                          :: rlocation, sourcelocation
    real, dimension(2)                          :: tbegin, tend
    real                                        :: d, tstf
    real                                        :: firstarrival_p, lastarrival_p, firstarrival_s, lastarrival_s 
    integer                                     :: nerr
    character, dimension(3), parameter          :: xyz = (/'x','y','z'/)
    type(varying_string)                        :: fn, basefn, offormat
    integer                                     :: ibeg, iend
    
     
    g_pn = 'ahfull'
    g_usage = 'usage: ahfull sources receivers material source-time-function dt output-filename-base output-format'

    ! process command line arguments and read input files
    if (iargc() /= 7) call usage()
    call readtable_arg( 1, 13, 1, sources_tab )
    call readtable_arg( 2, 5, 1, receivers_tab )
    call readtable_arg( 3, 3, 1, material_tab )
    call readtable_arg( 4, 2, 2, stf_tab )
    call readreal_arg( 5, tiny(dt), huge(dt), dt )
    call vs_getarg( 6, basefn )
    call vs_getarg( 7, offormat )
    

    if (offormat /= 'table' .and. offormat /= 'sac' .and. offormat /= 'mseed') &
        call die('unknown output fomat: ''' // offormat // '''; should be ''table'', ''sac'' or ''mseed''' )
    
    nstfsamples = size(stf_tab, 2)
    nsources = size(sources_tab, 2)
    nreceivers = size(receivers_tab, 2)

    ! prepare stf
    allocate( stf(nstfsamples) )
    stf(:) = stf_tab(2,:)
    tstf = (nstfsamples-1)*dt

    ! prepare sources
    allocate( sources(nsources) )
    sourcelocation = (/0,0,0/)
    do isource=1,nsources
        sources(isource) % tshift      = sources_tab(1,isource)
        sources(isource) % location(:) = sources_tab(2:4,isource)
        sources(isource) % mt(1,1)     = sources_tab(5,isource)
        sources(isource) % mt(2,2)     = sources_tab(6,isource)
        sources(isource) % mt(3,3)     = sources_tab(7,isource)
        sources(isource) % mt(1,2)     = sources_tab(8,isource)
        sources(isource) % mt(2,1)     = sources_tab(8,isource)
        sources(isource) % mt(1,3)     = sources_tab(9,isource)
        sources(isource) % mt(3,1)     = sources_tab(9,isource)
        sources(isource) % mt(2,3)     = sources_tab(10,isource)
        sources(isource) % mt(3,2)     = sources_tab(10,isource)
        sources(isource) % sf(:)       = sources_tab(11:13,isource)
        sourcelocation(:) =  sourcelocation(:) + sources(isource) % location(:)
    end do

    ! prepare receivers
    allocate( receivers(nreceivers) )
    do ireceiver=1,nreceivers
        receivers(ireceiver) % location(:) = receivers_tab(1:3,ireceiver)
        receivers(ireceiver) % nfflag      = (receivers_tab(4,ireceiver) /= 0.)
        receivers(ireceiver) % ffflag      = (receivers_tab(5,ireceiver) /= 0.)
    end do
 
    ! prepare material
    rho = material_tab(1,1)
    alpha = material_tab(2,1)
    beta = material_tab(3,1)

    call set_stf( es, stf, dt )
    call set_material( es, rho, alpha, beta )

    ! iterate over receivers
    do ireceiver=1,nreceivers
    
        ! det min and max distance between sources and receiver to calculate appropriate length
        firstarrival_p =  huge(firstarrival_p)
        firstarrival_s =  huge(firstarrival_s)
        lastarrival_p = -huge(lastarrival_p)
        lastarrival_s = -huge(lastarrival_s)
        
        do isource=1,nsources
            d = dist(sources(isource)%location,receivers(ireceiver)%location)
            
            firstarrival_p = min(firstarrival_p, &
                snapdown(d/alpha + sources(isource) % tshift,dt))
            lastarrival_p  = max(lastarrival_p, &
                snapup(d/alpha + sources(isource) % tshift + tstf,dt))
        
            firstarrival_s = min(firstarrival_s, &
                snapdown(d/beta + sources(isource) % tshift,dt))
            lastarrival_s  = max(lastarrival_s, &
                snapup(d/beta + sources(isource) % tshift + tstf,dt))
        end do
        
        if (lastarrival_p >= firstarrival_s) then   ! s and p time windows overlap
            nwindows = 1
            tbegin(1) = firstarrival_p
            tend(1) = lastarrival_s
        else                                        ! they are separated
            nwindows = 2
            tbegin(1) = firstarrival_p
            tend(1) = lastarrival_p
            tbegin(2) = firstarrival_s
            tend(2) = lastarrival_s
        end if
                
                
        ! iterate over time windows
        do iwindow=1,nwindows
        
            nsamples = (tend(iwindow) - tbegin(iwindow)) / dt + 1
            if (allocated( seismograms )) deallocate(seismograms)
            allocate( seismograms(3,nsamples) )
            seismograms(:,:) = 0.
            
            ! sum up elementary seismograms for every receiver-source combination
            do isource=1,nsources
                rlocation(:) = receivers(ireceiver)%location(:) - sources(isource)%location(:)
                call set_coords( es, rlocation, receivers(ireceiver)%nfflag, receivers(ireceiver)%ffflag )
                do n=1,3
                    do p=1,3
                        do q=1,3
                            call set_npq( es, n, p, q )
                            call add_elseis_mt( es, &
                                                tbegin(iwindow)-sources(isource)%tshift , &
                                                sources(isource) % mt(p,q), &
                                                seismograms(n,:) )
                        end do
                        call add_elseis_sf( es, &
                                            tbegin(iwindow)-sources(isource)%tshift , &
                                            sources(isource) % sf(p), &
                                            seismograms(n,:) )
                    end do
                end do
            end do
        
            do n=1,3
                ibeg = 1
                iend = nsamples
                
                fn = basefn // '-' // ireceiver // '-' // xyz(n) // '-' // iwindow // '.' // offormat
                call writeseismogram( fn, offormat, seismograms(n,ibeg:iend), real(tbegin(iwindow)+dt*(ibeg-1),8), dt, nerr )
                if (nerr /= 0) call die('failed to save seismogram to file: ' // fn )
                end do
        end do
        
        ! output filenames to stdout in a way, that gfdb_build understands
        do n=1,3
            fn = ""
            do iwindow=1,nwindows
                fn = fn // " " // basefn // '-' // ireceiver // '-' // xyz(n) // '-' // iwindow // '.' // offormat
            end do
            call put_line( fn )
            call flush(stdout)
        end do
        
    end do
    
    call cleanup()

  contains
    
    pure real function snapdown(t,dt) 
        real, intent(in) :: t, dt
        snapdown = floor(t/dt)*dt
    end function
    
    pure real function snapup(t,dt)
        real, intent(in) :: t, dt
        snapup = ceiling(t/dt)*dt
    end function
    
    function absmax_ind( array ) result(ind)
        real, dimension(:), intent(in) :: array
        integer                        :: ind,i
    
    ! get indice of first absmax
    
        real                           :: absmax
        
        absmax = 0.0
        ind = 1
        do i=1,size(array)
            if (abs(array(i)) > absmax) then 
                absmax=abs(array(i))
                ind = i
            end if
        end do
    
    end function absmax_ind
    
    subroutine readtable_arg( iarg, min_cols, min_rows, field )
    
        integer, intent(in) :: iarg, min_cols, min_rows
        real, dimension(:,:), allocatable, intent(out) :: field
        
    ! this sub reads a filename from argument iarg and then reads a table of real numbers contained in this file to field
    ! it is checked, that the table has at least min_cols columns and min_rows rows
        
        type(varying_string) :: fn
        integer :: ifile, iostat
        
        call vs_getarg( iarg, fn )
        call claim_unit( ifile )
        open( unit=ifile, file=char(fn), status='old', iostat=iostat )
        if (iostat /= 0) call die( 'can''t open file "' // fn // '"' )
        call readtable( field, iunit=ifile )
        close( ifile )
        call release_unit( ifile )
        if (size(field,1) <  min_cols) call die( 'expected at least ' // min_cols // ' column(s) in file "' // fn // '"' )
        if (size(field,1) <  min_rows) call die( 'expected at least ' // min_rows // ' row(s) in file "' // fn // '"' )

    end subroutine

    subroutine readreal_arg( iarg, min, max, val )
    
        integer, intent(in) :: iarg
        real, intent(in) :: min, max
        real, intent(out) :: val
    
    ! reads argument iarg from the argument list and converts it to real number val
    ! it is checked, that val is in the range [min, max]
        
        type(varying_string) :: vs
        
        call vs_getarg( iarg, vs )
        val = vs
        
        if (val < min) &
            call die( 'expected real value at argument ' // iarg  &
                 // ' in range [' // min // ',' // max // '], but got' // val // '.' )
        
    end subroutine readreal_arg

    
    pure function dist(a,b) result(r)
    
        real, intent(in), dimension(3) :: a,b
        real :: r
        
        r = sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)

    end function dist
    
    pure subroutine trans_xyz_lqt( origin, point, xyz )
        
        real, dimension(3), intent(in)       :: origin, point
        real, dimension(:,:), intent(inout)  :: xyz
        integer                              :: i
        real, dimension(3)                   :: er, eh, ep, rhp
        
    ! inplace tranform displacement vectors to local coordinate system 
    ! defined by radial (rel. to source), horizontal
    ! and a component perpendicular to the other two
        
        ! assuming size(xyz,2) == size(rhp,2)
        
        ! setup normal vectors for radial, horizontal and perpendicular components
        er(:) = point(:) - origin(:)
        call normalize(er)
        eh = cp(er,(/0.,0.,1./))
        if ( .not. any(eh /= 0.) ) eh = (/1.,0.,0./)
        call normalize(eh)
        ep = cp(er,eh)
        
        ! transform each element in xyz
        do i=1,size(xyz,2)
            rhp(1) = sp(xyz(:,i),er)
            rhp(2) = sp(xyz(:,i),eh)
            rhp(3) = sp(xyz(:,i),ep)
            xyz(:,i) = rhp(:)
        end do
        
    end subroutine trans_xyz_lqt
    
    pure subroutine normalize( vec )
        real, intent(inout), dimension(3) :: vec
        real :: r 
        r = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
        vec(1:3) = vec(1:3) / r
    end subroutine normalize
    
    pure function sp( a, b ) result( c )
        real, intent(in), dimension(3) :: a,b
        real                           :: c
    ! scalar product
        c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    end function sp
    
    pure function cp( a, b ) result( c )
        real, intent(in), dimension(3) :: a,b
        real, dimension(3)             :: c
    ! cross product
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function
    
end program
