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

module gfdb

    use util
    use sparse_trace
    use better_varying_string
    use hdf5
    use interpolation
    use omp_lib
    use gfdb_io_hdf

    implicit none

    public
    
    integer, parameter :: nblockx_default = 128            ! number of traces in one interpolation block in distance (power of 2)
    integer, parameter :: nblockx_overlap_default = 32     ! overlap between blocks (must be even)
    integer, parameter :: nblockx_payload_default = 128-32 ! this is the width of the usable part of the blocks
    
    integer, parameter :: nblockz_default = 32             ! number of traces in one interpolation block in depth (power of 2)
    integer, parameter :: nblockz_overlap_default = 8      ! overlap between blocks (must be even)
    integer, parameter :: nblockz_payload_default = 32-8   ! this is the width of the usable part of the blocks
    
    type :: t_trace_p 
    
      ! holds a pointer to a trace
      ! this is needed to make a array of pointers in fortran
      
        type(t_trace), pointer :: p => null()
        
    end type    
    
    
    type :: t_chunk
        
      ! this is an entity corresponding to a single file on disk
      ! it holds nxc*nz*ng traces
      
        integer :: ichunk ! this chunks number
        integer :: nxc = 0 ! number of receiver-distances in this chunk
        integer :: nz = 0  ! number of source depths
        integer :: ng = 0  ! number of greens functions (or elementary seismograms) (8/10 ff/nf)
        
        integer :: nipx = 1       ! total traces / real traces
                                  ! (increment between 2 non-interpolated traces)
        integer :: nipz = 1       ! total traces / real traces
                                  ! (increment between 2 non-interpolated traces)
        
        logical :: writemode = .false.
        logical :: readmode = .false.
        type(varying_string) :: filename
        integer(hid_t) :: file = 0
        integer(hid_t) :: dataset_index = 0
        
      ! bytes allocated inside traces cache        
        integer(kind=8) :: nbytes = 0
        integer(kind=8) :: last_access = 0

      ! references to the datasets are stored in this array for quick access:
        type(hobj_ref_t_f), dimension(:,:,:), allocatable :: references    ! (ng,nz/nipz,nxc/nipx)
      
      ! pointers to traces, where allready used traces live
        type(t_trace_p), dimension(:,:,:), allocatable :: traces           ! (ng,nz,nxc)
  
      ! these are buffers, needed by chunk_get_trace,
      ! they have been put here, so that their allocation can be reused
        integer, dimension(:), allocatable :: ofs, pofs
        real, dimension(:), allocatable :: packed

      ! access stack 
      ! used to keep track of which chunks have been used most recently
        integer :: accessed_earlier = 0
        integer :: accessed_later = 0
        
    end type
    

    type :: t_gfdb
        
      ! the t_gfdb corresponds to one green function database for a specific sampling rate.
        real    :: dt = 0.
        real    :: dx = 0.
        real    :: dz = 0.
        real    :: firstx = 0.
        real    :: firstz = 0.
        
      ! filenames of cache files are made by appending .index or .<ichunk>.chunk
      ! to filenamebase:
        
        type(varying_string) :: filenamebase

      ! the following are the wanted chunk sizes,
      ! the nx on the last chunk may differ from the others.
      
        integer :: nx  = 0     ! number of receiver-distances
        integer :: nxc = 0     ! wanted number of receiver-distances per chunk
                               ! (nxc of last chunk differs!)
        integer :: nz  = 0     ! number of source depths 
        integer :: ng  = 0     ! number of greens functions (= 8 or 10)
        
        integer :: nipx = 1    ! total traces / real traces
                               ! (increment between 2 non-interpolated traces)
        integer :: nipz = 1    ! total traces / real traces
                               ! (increment between 2 non-interpolated traces)
                               
      ! always a block of nblockx traces will be interpolated at once. 
      ! interpolated traces near to the edge of the block contain artifacts,
      ! so some overlap is needed.
        
        integer :: nblockx = nblockx_default                   ! number of traces in one iterpolation block in distance (power of 2)
        integer :: nblockx_overlap = nblockx_overlap_default   ! overlap between blocks (must be even)
        integer :: nblockx_payload =  nblockx_payload_default  ! this is the width of the usable part of the blocks
        
        integer :: nblockz = nblockz_default                   ! number of traces in one iterpolation block in depth (power of 2)
        integer :: nblockz_overlap = nblockz_overlap_default   ! overlap between blocks (must be even)
        integer :: nblockz_payload = nblockz_payload_default   ! this is the width of the usable part of the blocks
        
        integer(kind=8) :: nbytes = 0        ! bytes allocated inside traces cache
        integer(kind=8) :: nbytes_limit = 0
        
      ! cache contains nchunks chunks:
        integer :: nchunks = 0
        
        type(t_chunk), dimension(:), allocatable :: chunks
        
        type(t_trace_p), dimension(:), allocatable :: interpolated_traces

        integer :: accessed_latest = 0
        integer :: accessed_earliest = 0

    end type
    
    
    
    public gfdb_init, gfdb_save_trace, gfdb_close, gfdb_get_trace, gfdb_uncache_trace
    public gfdb_get_indices, gfdb_infos, gfdb_dump_infomap, gfdb_dump_contents
    public gfdb_dump_missing, gfdb_cached_traces_memory, gfdb_set_cached_traces_memory_limit
    public gfdb_housekeeping
    
  contains
  
    subroutine gfdb_cleanup()
        
        call gfdb_io_deinit()
        
    end subroutine
  
    subroutine gfdb_init( c, fnbase, nchunks, nx, nz, ng, dt, dx, dz, firstx, firstz, nipx, nipz )
    
      ! initialize a new cache object
      ! either by given parameters or by opening file
    
        type(t_gfdb), intent(inout) :: c
        type(varying_string), intent(in) :: fnbase
        integer, intent(in), optional :: nchunks, nx,nz,ng
        real, intent(in), optional :: dt, dx, dz
        real, intent(in), optional :: firstx, firstz
        integer, intent(in), optional :: nipx
        integer, intent(in), optional :: nipz

        integer :: ichunk, nxcthis
        type(varying_string) :: filename
        logical :: ok
        
        call gfdb_io_init(ok)
        if (.not. ok) call die()
        
        if (allocated(c%chunks)) call gfdb_destroy( c )

        if (present(nchunks) .and. present(nx) .and. present(nz) .and. &
            present(ng) .and. present(dt) .and. present(dx) .and. present(dz)) then
            
            c%filenamebase = fnbase
            c%nx = nx
            if (nx >= nchunks) then
                c%nchunks = nchunks
            else 
                c%nchunks = nx
            end if
            c%nxc = nx/c%nchunks + 1
            if (c%nxc > nx) c%nxc = nx
            do while (nx - c%nxc*(c%nchunks-1) <= 0)
                c%nxc = c%nxc - 1
            end do
            c%nz = nz
            c%ng = ng
            c%dt = dt
            c%dx = dx
            c%dz = dz
            if (present(firstx)) c%firstx = firstx
            if (present(firstz)) c%firstz = firstz
            
        else   ! read cache metainformation from file
                
            c%filenamebase = fnbase
            filename = fnbase // ".index"
            c%nipx = 1
            c%nipz = 1
            
            if (present(nipx)) c%nipx = nipx
            if (present(nipz)) c%nipz = nipz
            
            call gfdb_io_read_index(filename, c%dt,c%dx,c%dz, c%firstx, c%firstz, &
                                                c%nchunks, c%nx, c%nxc, c%nz, c%ng, ok)
            if (.not. ok) call die()
            
          ! if interpolation is wanted add disired amount of oversampling
            if (c%nipx .ne. 1) then
                c%nx  = c%nx * c%nipx
                c%nxc = c%nxc * c%nipx
                c%dx  = c%dx / c%nipx
                c%nblockx = nblockx_default
                c%nblockx_payload = nblockx_payload_default
                c%nblockx_overlap = nblockx_overlap_default
            else 
                c%nblockx = 1
                c%nblockx_payload = 1
                c%nblockx_overlap = 0
            end if
            if (c%nipz .ne. 1) then
                c%nz  = c%nz * c%nipz
                c%dz  = c%dz / c%nipz
                c%nblockz = nblockz_default
                c%nblockz_payload = nblockz_payload_default
                c%nblockz_overlap = nblockz_overlap_default
            else 
                c%nblockz = 1
                c%nblockz_payload = 1
                c%nblockz_overlap = 0
            end if
        end if
        
        allocate( c%chunks(c%nchunks) )

      ! initialize associated chunks
        do ichunk=1,c%nchunks
            nxcthis = c%nxc
            if (ichunk == c%nchunks) nxcthis = c%nx-(ichunk-1)*c%nxc
            call chunk_init( c%chunks(ichunk), &
                             c%filenamebase // "." // ichunk // ".chunk", &
                             ichunk, nxcthis,c%nz,c%ng, c%nipx, c%nipz)
        end do
        
        
        if (present(nchunks) .and. present(nx) .and. present(nz) .and. &
            present(ng) .and. present(dt) .and. present(dx) .and. present(dz)) then
            call gfdb_create( c )
        end if
    end subroutine
    
    subroutine gfdb_infos( c, indexmemory_per_chunk, ntraces_in_chunks )
        
        type(t_gfdb), intent(inout) :: c
        integer, intent(out) :: indexmemory_per_chunk
        integer, intent(inout), dimension(:,:), allocatable :: ntraces_in_chunks
        
        integer :: ntraces, ntraces_used, ichunk, ntraces_per_chunk
        type(hobj_ref_t_f) :: dummyref
        
        dummyref%ref = 0
        
        ntraces_per_chunk = c%nxc*c%nz*c%ng;
        indexmemory_per_chunk = kind(dummyref%ref)*ntraces_per_chunk;
        
        if (allocated( ntraces_in_chunks )) deallocate( ntraces_in_chunks )
        allocate( ntraces_in_chunks(2,c%nchunks) )
        do ichunk=1,c%nchunks
            call chunk_infos( c, c%chunks(ichunk), ntraces, ntraces_used )
            ntraces_in_chunks(1,ichunk) = ntraces_used
            ntraces_in_chunks(2,ichunk) = ntraces
        end do
        
    end subroutine
    
    subroutine gfdb_dump_infomap( c, iunit )
        type(t_gfdb), intent(inout) :: c
        integer, intent(in) :: iunit
        integer :: ichunk
        
        do ichunk=1,c%nchunks
            call chunk_dump_infomap( c, c%chunks(ichunk), iunit,c%nxc/c%nipx  )
        end do
        
    end subroutine

    pure function gfdb_cached_traces_memory( db )

        type(t_gfdb), intent(in) :: db
        integer :: gfdb_cached_traces_memory

        gfdb_cached_traces_memory = db%nbytes
 
    end function

    pure subroutine gfdb_set_cached_traces_memory_limit( db, nbytes_limit )

        type(t_gfdb), intent(inout)        :: db
        integer(kind=8), intent(in)     :: nbytes_limit

        db%nbytes_limit = nbytes_limit 

    end subroutine

    pure subroutine gfdb_accesslist_remove(db, ichunk)

        type(t_gfdb), intent(inout)  :: db
        integer, intent(in) :: ichunk

        integer :: ilate, iearl

        if (ichunk == 0) return

        ilate =  db%chunks(ichunk)%accessed_later
        iearl = db%chunks(ichunk)%accessed_earlier

        if (ilate .eq. 0 .and. iearl .eq. 0) return   ! not in list

        db%chunks(ichunk)%accessed_later = 0
        db%chunks(ichunk)%accessed_earlier = 0
                
        if (ilate .ne. 0) then
            db%chunks(ilate)%accessed_earlier = iearl
        else
            db%accessed_latest = iearl
        end if
        
        if (iearl .ne. 0) then
            db%chunks(iearl)%accessed_later = ilate
        else
            db%accessed_earliest = ilate
        end if

    end subroutine

    subroutine gfdb_accesslist_push_latest(db, ichunk)

        type(t_gfdb), intent(inout)  :: db
        integer, intent(in) :: ichunk

        integer :: ioldlatest

        if (ichunk == 0) return
        if (ichunk == db%accessed_latest) return ! is already at latest position
        
        ioldlatest = db%accessed_latest
        if (ioldlatest .ne. 0) then
            db%chunks(ioldlatest)%accessed_later = ichunk
        else
            db%accessed_earliest = ichunk
        end if

        db%chunks(ichunk)%accessed_later = 0
        db%chunks(ichunk)%accessed_earlier = ioldlatest
        db%accessed_latest = ichunk

    end subroutine

    subroutine gfdb_housekeeping(db)

      ! call periodically at times, when it is safe to uncache traces.
      ! 
        
        type(t_gfdb), intent(inout)        :: db

        integer :: ichunk, ichunknext

        if (db%nbytes_limit > 0 .and. db%nbytes_limit < db%nbytes) then

            call warn( "gfdb: cache memory limit reached. emptying..." )

            ichunk = db%accessed_earliest
            do while (ichunk .ne. 0 .and. db%nbytes > db%nbytes_limit/2) 
                ichunknext = db%chunks(ichunk)%accessed_later
                call chunk_close( db, db%chunks(ichunk) )
                ichunk = ichunknext
            end do

            if (db%nbytes > db%nbytes_limit/2) then
                call die("gfdb: cache management failure; if this program had no errors, this could not have happened...")
            end if
            
        end if
        


    end subroutine

    subroutine chunk_infos( db, c, ntraces, ntraces_used )
    
        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c
        integer, intent(out) :: ntraces, ntraces_used
        integer :: ixc, iz, ig

        ntraces = c%nxc/c%nipx*c%nz/c%nipz*c%ng
        
        call chunk_open_read( db, c ) 
        
      ! count number of valid references to get number of traces
        ntraces_used = 0
        do ixc=1,c%nxc/c%nipx
            do iz=1,c%nz/c%nipz
                do ig=1,c%ng
                    if ( c%references(ig,iz,ixc)%ref /= 0 ) then
                        ntraces_used = ntraces_used + 1
                    end if
                end do
            end do
        end do
        call chunk_close(db, c )
        
    end subroutine
    
    
    subroutine chunk_dump_infomap( db, c, iunit, nxcfull )
        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: iunit, nxcfull
        integer :: ixc, iz, ig
        integer :: ialloc
        
        call chunk_open_read( db, c ) 
        
        do ixc=1,c%nxc/c%nipx
            do iz=1,c%nz/c%nipz
                ialloc = 0
                do ig=1,c%ng
                    if ( c%references(ig,iz,ixc)%ref /= 0 ) then
                        ialloc = ialloc+1
                    end if
                end do
                
                write (iunit,*) (c%ichunk-1)*nxcfull+ixc, iz, ialloc
                
            end do
        end do
        
        call chunk_close( db,c )
    
    end subroutine
    
    subroutine gfdb_dump_contents( db, iunit )
        type(t_gfdb), intent(inout) :: db
        integer, intent(in) :: iunit
        integer :: ichunk
        
        do ichunk=1,db%nchunks
            call chunk_dump_contents( db, db%chunks(ichunk), iunit,db%nxc/db%nipx  )
        end do
        
    end subroutine
    
    subroutine chunk_dump_contents( db, c, iunit, nxcfull )
        type(t_gfdb), intent(inout) :: db    
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: iunit, nxcfull
        integer :: ixc,ix, iz, ig
        real :: x,z
        
        call chunk_open_read( db, c ) 
        
        do ixc=1,c%nxc/c%nipx
            ix = (c%ichunk-1)*nxcfull+ixc
            do iz=1,c%nz/c%nipz
                call gfdb_get_position( db, ix, iz, x, z )
                do ig=1,c%ng
                    if ( c%references(ig,iz,ixc)%ref /= 0 ) then
                        write (iunit,*) x, z, ig
                    end if
                end do
            end do
        end do
        
        call chunk_close(db,c )
    
    end subroutine
    
    subroutine gfdb_dump_missing( db, iunit )
        type(t_gfdb), intent(inout) :: db
        integer, intent(in) :: iunit
        integer :: ichunk
        
        do ichunk=1,db%nchunks
            call chunk_dump_missing( db, db%chunks(ichunk), iunit,db%nxc/db%nipx  )
        end do
        
    end subroutine
    
    subroutine chunk_dump_missing( db, c, iunit, nxcfull )

        type(t_gfdb), intent(inout) :: db    
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: iunit, nxcfull
        integer :: ixc,ix, iz, ig
        real :: x,z
        
        call chunk_open_read( db, c ) 
        
        do ixc=1,c%nxc/c%nipx
            ix = (c%ichunk-1)*nxcfull+ixc
            do iz=1,c%nz/c%nipz
                call gfdb_get_position( db, ix, iz, x, z )
                do ig=1,c%ng
                    if ( c%references(ig,iz,ixc)%ref == 0 ) then
                        write (iunit,*) x, z, ig
                    end if
                end do
            end do
        end do
        
        call chunk_close(db,c )
    
    end subroutine

    subroutine gfdb_dump_stats( db, iunit )

        type(t_gfdb), intent(inout) :: db
        integer, intent(in) :: iunit
        integer :: ix, iz, ig
        real :: x,z
        type(t_trace), pointer   :: tracep
        type(t_strip) :: strip
        integer :: nflat
        real, dimension(2) :: data_range
        integer, dimension(2) :: span
        
        do ix = 1,db%nx/db%nipx
            do iz=1,db%nz/db%nipz
                call gfdb_get_position( db, ix, iz, x, z )
                do ig=1,db%ng
                    call gfdb_get_trace(db, ix, iz, ig, tracep)
                    if (associated(tracep)) then
                        call trace_unpack(tracep, strip)
                        call stats(strip, data_range, nflat)
                        span = strip_span(strip)
                        write (iunit,*) x, z, ig, &
                            span(1)*db%dt, span(2)*db%dt, &
                            data_range(1), data_range(2), nflat*db%dt
                    end if
                    call gfdb_uncache_trace(db, ix, iz, ig)
                end do
            end do
        end do

    end subroutine

    subroutine stats( strip, data_range, nflat )

        type(t_strip), intent(in) :: strip
        real, dimension(2), intent(out) :: data_range
        integer, intent(out) :: nflat
        
        integer :: n, i
        real, dimension(size(strip%data)) :: ddata
        real, dimension(2) :: ddata_range
        real, parameter :: flat_level = 0.01

        data_range(1) = minval(strip%data)
        data_range(2) = maxval(strip%data)

        n = size(strip%data)
        ddata(:n-1) = strip%data(:n-1)
        ddata_range(1) = minval(ddata(:n-1))
        ddata_range(2) = maxval(ddata(:n-1))
        nflat = 0
        do i=1,n-1
            if (abs(ddata(i)) .gt. flat_level*(ddata_range(2)-ddata_range(1))) then
                nflat = i-1
                exit
            end if
        end do

    end subroutine

    subroutine gfdb_destroy( db )
    
      ! destroy gfdb object db
      ! free all memory associated with db
      ! reset everything to 0
      ! it should be safe to destroy a not initialized cache
         
        type(t_gfdb), intent(inout) :: db
        integer :: ichunk
        
        call gfdb_close( db)

        if (allocated(db%chunks)) then
            do ichunk=1,db%nchunks
                call chunk_destroy( db,  db%chunks(ichunk) )
            end do
            deallocate( db%chunks )
        end if
        
        call delete(db%filenamebase)
        db%nchunks = 0
        db%nx = 0
        db%nxc = 0
        db%nz = 0
        db%ng = 0
        db%dt = 0.
        db%dx = 0.
        db%dz = 0.
        db%nipx = 1
        db%nipz = 1
        db%nbytes = 0
        
    end subroutine
    
    subroutine chunk_destroy(db, c)
        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c
        
        call chunk_close( db, c )
        c%nxc = 0
        c%nz = 0
        c%ng = 0
        call delete(c%filename)
        c%writemode = .false.
        c%readmode = .false.
        c%nipx = 1
        c%nipz = 1
        c%nbytes = 0

        if (allocated(c%ofs)) deallocate(c%ofs)
        if (allocated(c%pofs)) deallocate(c%pofs)
        if (allocated(c%packed)) deallocate(c%packed)
    
    end subroutine
    
    subroutine chunk_init( c, filename, ichunk, nxc,nz,ng, nipx, nipz )
    
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: ichunk,nxc,nz,ng
        type(varying_string), intent(in) :: filename
        integer, intent(in), optional :: nipx
        integer, intent(in), optional :: nipz
        
        c%filename = filename
        c%ichunk = ichunk
        c%nxc = nxc
        c%nz = nz
        c%ng = ng
        
        c%nipx = 1
        c%nbytes = 0
        if (present(nipx)) c%nipx = nipx
        if (present(nipz)) c%nipz = nipz
        
        c%file = 0
        c%dataset_index = 0
        
    end subroutine
    
    subroutine gfdb_create( c )
    
        type(t_gfdb), intent(inout) :: c
        
        integer :: ichunk
        integer(hid_t) :: file
        integer(hid_t) :: egal
        

        type(varying_string) :: filename
        logical :: ok        

        if (.not. allocated(c%chunks)) return

        filename = c%filenamebase // ".index"
        
        call gfdb_io_create_index(filename, &
                c%dt, c%dx*c%nipx, c%dz*c%nipz, c%firstx, c%firstz, c%nchunks, &
                c%nx/c%nipx, c%nxc/c%nipx, c%nz/c%nipz, c%ng, ok)
        
        if (.not. ok) call die()

        do ichunk=1,c%nchunks
            call chunk_create( c%chunks(ichunk))
        end do
    
    end subroutine
    
    subroutine chunk_create( c )
    
        type(t_chunk), intent(inout) :: c
        logical :: ok

        call gfdb_io_create_chunk(c%filename, c%ng, c%nxc/c%nipx, c%nz/c%nipz, &
            int(c%nxc/c%nipx*(floor(log10(real(c%nxc/c%nipx)))+2),SIZE_T), ok)

        if (.not. ok) call die()
             
    end subroutine
    
    subroutine gfdb_save_trace( db, ix, iz, ig, data )
    
      ! save a trace to the database
        
        type(t_gfdb), intent(inout) :: db
        integer, intent(in) :: ix,iz,ig
        type(t_trace), intent(in) :: data
        
        integer :: ichunk, ixc
        
        if (.not. allocated(db%chunks)) return
                        
        call gfdb_index_to_chunk( db, ix, ichunk, ixc )
        
        call chunk_save_trace( db, db%chunks(ichunk), ixc,iz,ig, data )
        
    end subroutine
    
    subroutine chunk_save_trace( db, c, ixc, iz, ig, data )
    
      ! save a trace to this chunks file
        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: ixc,iz,ig
        type(t_trace), intent(in) :: data
        
        type(hobj_ref_t_f), dimension(1)  :: reference
        integer :: packed_size, nstrips
        integer :: ixc_file, iz_file
        logical :: ok

        if (ixc > c%nxc .or. ixc < 1) return
        if (iz > c%nz .or. iz < 1) return
        if (ig > c%ng .or. ig < 1) return
        
      ! cannot save interpolated trace this way.
        if (mod(ixc-1,c%nipx) /= 0) return
        if (mod(iz-1,c%nipz) /= 0) return
        ixc_file = (ixc-1)/c%nipx + 1
        iz_file = (iz-1)/c%nipz + 1
        
        if (trace_is_empty(data)) return
        
        call chunk_open_write( db, c ) ! ensure chunk is open for write
        
      ! check that there is no trace already stored at that position
        if (c%references(ig,iz_file,ixc_file)%ref /= 0) then
            call warn( "gfdb: chunk_save_trace() (chunk "//c%ichunk// &
                       "): trace already exists for index "//&   
                       "("//ixc_file//","//iz_file//","//ig//")" )
            return
        end if
   
      ! pack data to a single, continuous array
      ! generate arrays with offsets in the packed array (pofs)
      ! and offsets of the data strips (ofs)
        
        call trace_to_storable( data, c%packed, packed_size, c%pofs, c%ofs, nstrips )

        call gfdb_io_save_trace(c%filename, c%file, c%dataset_index, &
                ixc_file, iz_file, ig, &
                c%packed, packed_size, c%pofs, c%ofs, nstrips, &
                int(c%nz/c%nipz*(floor(log10(real(c%nz/c%nipz)))+2),SIZE_T), &
                int(8*(floor(log10(8.))+2),SIZE_T), &
                reference, ok)

        if (.not. ok) call die()

        c%references(ig,iz_file,ixc_file) = reference(1)

    end subroutine
    
    subroutine gfdb_get_indices( c, x, z, ix, iz )
    
      ! make indices suitable for get_trace, based on spacing of gf's
      
        type(t_gfdb), intent(in) :: c
        real, intent(in) :: x, z
        integer, intent(out) :: ix, iz
        
        ix = nint((x-c%firstx)/c%dx)+1
        iz = nint((z-c%firstz)/c%dz)+1
        
    end subroutine
    
    subroutine gfdb_get_indices_bilin( c, x, z, nxundersample, nzundersample, ix, iz, dix, diz )
    
      ! make upper and lower indices suitable for get_trace, based on spacing of gf's
      ! additionally return fractional fading factor.
      ! if nxundersample or nzundersample are not 1, get indices into downsampled grid.
      
        type(t_gfdb), intent(in) :: c
        real, intent(in) :: x, z
        integer, intent(in) :: nxundersample, nzundersample
        integer, dimension(2), intent(out) :: ix, iz
        real, intent(out) :: dix, diz
        
        
        ix(1) = int(floor((x-c%firstx)/(c%dx*nxundersample)))*nxundersample+1
        iz(1) = int(floor((z-c%firstz)/(c%dz*nzundersample)))*nzundersample+1
        ix(2) = ix(1)+nxundersample
        iz(2) = iz(1)+nzundersample
    
        dix = (x-c%firstx-(ix(1)-1)*c%dx)/(c%dx*nxundersample)
        diz = (z-c%firstz-(iz(1)-1)*c%dz)/(c%dz*nzundersample)
        
    end subroutine
    
    subroutine gfdb_get_position( c, ix, iz, x, z )
    
      ! inverse of gfdb_get_indices()
      
        type(t_gfdb), intent(in) :: c
        integer, intent(in) :: ix, iz
        real, intent(out) :: x, z

        x = c%firstx + (ix-1)*c%dx
        z = c%firstz + (iz-1)*c%dz
        
    end subroutine
    
    subroutine gfdb_get_trace( db, ix, iz, ig, tracep )
    
      ! get a pointer to the trace indexed by (ix,iz,ig) from the database.
      ! if the indices are out of bounds, null is returned
      ! if there is no stored trace at this location (or if the stored 
      ! trace at this location is empty), an  empty trace is returned 
      ! (check with trace_is_empty()).
      
      ! this simply routes the request to the appropriate chunk
      
        type(t_gfdb), intent(inout)            :: db
        integer, intent(in)                   :: ix, iz, ig
        type(t_trace), pointer   :: tracep
        
        integer :: ichunk, ixc
        
        tracep => null()  
      
        if (.not. allocated(db%chunks) .or. &
            ix > db%nx .or. ix < 1 .or. &
            iz > db%nz .or. iz < 1 .or. &
            ig > db%ng .or. ig < 1) then
            call warn("gfdb: gfdb_get_trace(): invalid request: out of bounds: "//&    
                      "("//ix//","//iz//","//ig//")")
            return
        end if
        
        call gfdb_index_to_chunk( db, ix, ichunk, ixc )
        
        !$omp critical
        call chunk_get_trace( db, db%chunks(ichunk), ixc,iz,ig, tracep )
        !$omp end critical

    end subroutine
    
    subroutine gfdb_get_trace_bilin( db, ix, iz, ig, dix, diz, tracep )
    
        type(t_gfdb), intent(inout)           :: db
        integer, dimension(2), intent(in)     :: ix, iz
        integer, intent(in)                   :: ig
        real, intent(in)                      :: dix, diz
        type(t_trace), pointer                :: tracep
    
      ! use bilinear interpolation to get a trace between gfdb grid nodes.
      ! ix, iz, dix and diz are indices and offsets as produced by gfdb_get_indices_bilin()
      ! 
      ! if dix==0 and diz==0 
      
        
        type(t_trace), pointer  :: t00, t01, t10, t11
        integer, dimension(2)   :: span
        integer                 :: iipt, nipt

        tracep => null()
        t00 => null()
        t01 => null()
        t10 => null()
        t11 => null()

      ! if we are exactly at a grid node, no interpolation is needed...
        if (dix .eq. 0. .and. diz .eq. 0.) then
            call gfdb_get_trace( db, ix(1), iz(1), ig, tracep )
            return
        end if
        
      ! get the traces
        call gfdb_get_trace( db, ix(1), iz(1), ig, t00 )
        call gfdb_get_trace( db, ix(1), iz(2), ig, t01 )
        call gfdb_get_trace( db, ix(2), iz(1), ig, t10 )
        call gfdb_get_trace( db, ix(2), iz(2), ig, t11 )
        
        if (.not. (associated(t00) .and. associated(t01) .and. &
                   associated(t10) .and. associated(t11))) then
            tracep => null()
            return
        end if


      ! get union of spans
        span(1) = min( t00%span(1), t01%span(1), t10%span(1), t11%span(1) )
        span(2) = max( t00%span(2), t01%span(2), t10%span(2), t11%span(2) )
      
      ! prepare the global buffer traces which will be reused on first invocation
        if (.not. allocated( db%interpolated_traces )) then
            nipt = omp_get_max_threads()
            allocate( db%interpolated_traces(nipt) )
            do iipt=1,nipt
                db%interpolated_traces(iipt)%p => null()
            end do
        end if

        if (omp_get_num_threads() > size(db%interpolated_traces)) then
            call die("number of threads exceeds assumed maximum number of threads")
        end if

        iipt = omp_get_thread_num()+1
        
        if (.not. associated( db%interpolated_traces(iipt)%p ) ) then
            allocate( db%interpolated_traces(iipt)%p )
        end if

        tracep => db%interpolated_traces(iipt)%p

        if ( trace_is_empty( tracep ) ) then
            call trace_create_simple_nodata( tracep, span )
        else ! resize the contained buffer
            !!! relying on internals of trace_t here...
            call resize( tracep%strips(1)%data, span(1), span(2)-span(1)+1 )
            tracep%span(:) = span(:)
        end if
        

      ! summation
        tracep%strips(1)%data(:) = 0.
        
        call trace_multiply_add_nogrow( t00, tracep%strips(1)%data, span, (1.-dix)*(1.-diz) )
        call trace_multiply_add_nogrow( t01, tracep%strips(1)%data, span, (1.-dix)*diz )
        call trace_multiply_add_nogrow( t10, tracep%strips(1)%data, span, dix*(1.-diz) )
        call trace_multiply_add_nogrow( t11, tracep%strips(1)%data, span, dix*diz )
          
    end subroutine
    
    subroutine chunk_get_trace( db, c, ixc,iz,ig, tracep )
    
      ! get a pointer to the trace indexed by (ixc,iz,ig) from chunk c.
      ! if the indices are out of bounds, null is returned
      ! if there is no stored trace at this location (or if the stored 
      ! trace at this location is empty), an  empty trace is returned 
      ! (check with trace_is_empty()).
    
        type(t_gfdb), intent(inout)         :: db
        type(t_chunk), intent(inout)        :: c
        integer, intent(in)                 :: ixc, iz,ig
        type(t_trace), pointer :: tracep
        
        integer :: packed_size, nstrips
        integer :: ixc_file, ix
        integer :: iz_file
        integer :: nbytes
        logical :: ok
                
        tracep => null()  
      
        if (ixc > c%nxc .or. ixc < 1 .or. &
            iz > c%nz .or. iz < 1 .or. &
            ig > c%ng .or. ig < 1) then
            call warn("gfdb: chunk_get_trace(): invalid request: out of bounds : "//&   
                      "("//ixc//","//iz//","//ig//")")
            return
        end if
        
        ixc_file = (ixc-1)/c%nipx+1
        iz_file = (iz-1)/c%nipz+1
        
      ! ensure that chunk is open for read
        call chunk_open_read( db, c ) 
        
        call gfdb_accesslist_remove( db, c%ichunk )
        call gfdb_accesslist_push_latest(db, c%ichunk )
     
      ! directly return pointer if it is already in memory
        if (associated(c%traces(ig,iz,ixc)%p)) then
            tracep => c%traces(ig,iz,ixc)%p
            return
        end if
      
      ! trigger block interpolation if needed
        ix = (c%ichunk-1)*db%nxc+ixc
        if (mod(ix-1, db%nipx ) /= 0 .or. &
            mod(iz-1, db%nipz ) /= 0 ) then
            call gfdb_interpolate_block(db, ix, iz )
            tracep => c%traces(ig,iz,ixc)%p
            return
        end if

        if (c%references(ig,iz_file,ixc_file)%ref == 0) then
            call warn( "gfdb: chunk_get_trace() (chunk "//c%ichunk// &
                       "): no trace available for index "//&   
                       "("//ixc_file//","//iz_file//","//ig//")" )
            return
        end if

        call gfdb_io_get_trace(c%dataset_index, &
                c%references(ig,iz_file,ixc_file), &
                c%packed, packed_size, c%pofs, c%ofs, nstrips, ok)
               
        if (.not. ok) call die()

        allocate(c%traces(ig,iz,ixc)%p)
        tracep => c%traces(ig,iz,ixc)%p

        call trace_destroy( tracep )
        
      ! convert to sparse trace  
        call trace_from_storable( tracep, c%packed(1:packed_size), &
                                  c%pofs(1:nstrips), c%ofs(1:nstrips) )
        
        nbytes = trace_size_bytes( tracep )
        c%nbytes = c%nbytes + nbytes
        db%nbytes = db%nbytes + nbytes
    
    end subroutine
    
    subroutine gfdb_set_trace( db, ix, iz, ig, data, span )
    
      ! insert an interpolated trace into the gfdb.
      ! this does not save the trace as gfdb_save_trace does.
      
      ! this simply routes the request to the appropriate chunk
      
        type(t_gfdb), intent(inout)           :: db
        integer, intent(in)                   :: ix, iz, ig
        real, dimension(:), intent(in)        :: data
        integer, dimension(:), intent(in)     :: span
        
        integer :: ichunk, ixc
        
        if (.not. allocated(db%chunks) .or. &
            ix > db%nx .or. ix < 1 .or. &
            iz > db%nz .or. iz < 1 .or. &
            ig > db%ng .or. ig < 1) then
            call warn("gfdb: gfdb_set_trace(): invalid request: out of bounds: "//&    
                      "("//ix//","//iz//","//ig//")")
            return
        end if
        call gfdb_index_to_chunk( db, ix, ichunk, ixc )
        
        call chunk_set_trace( db, db%chunks(ichunk), ixc,iz,ig, data, span )
    
    end subroutine
    
    subroutine chunk_set_trace( db, c, ixc, iz, ig, data, span )
    
      ! insert an interpolated trace into this chunk
      ! this does not save the trace as gfdb_save_trace does.
        
        type(t_gfdb), intent(inout)              :: db
        type(t_chunk), intent(inout)          :: c
        integer, intent(in)                   :: ixc, iz, ig
        real, dimension(:), intent(in)        :: data
        integer, dimension(:), intent(in)     :: span

        type(t_trace), pointer                :: tracep
        integer :: nbytes
        
        tracep => null()  
      
        if (ixc > c%nxc .or. ixc < 1 .or. &
            iz > c%nz .or. iz < 1 .or. &
            ig > c%ng .or. ig < 1) then
            call warn("gfdb: chunk_set_trace(): invalid request: out of bounds : "//&   
                      "("//ixc//","//iz//","//ig//")")
            return
        end if
                
      ! ensure that chunk is open for read  
      ! ( here, nothing is read, but this ensures, that the pointer arrays are allocated...)
        call chunk_open_read( db, c )
     
      ! fail if trace already in
        if (associated(c%traces(ig,iz,ixc)%p)) then
            call warn( "gfdb: chunk_set_trace() (chunk "//c%ichunk// &
                       "): trace already exists for index "//&   
                       "("//ixc//","//iz//","//ig//")" )
            return
        end if
        
        
        allocate(c%traces(ig,iz,ixc)%p)
        tracep => c%traces(ig,iz,ixc)%p
                
        call trace_create_simple( tracep, data, span )

        nbytes = trace_size_bytes( tracep )
        c%nbytes = c%nbytes + nbytes
        db%nbytes = db%nbytes + nbytes

    end subroutine
    
    subroutine gfdb_interpolate_block( db, ix_in, iz_in )
    
        type(t_gfdb), intent(inout)         :: db
        integer, intent(in)                 :: ix_in, iz_in
        
      ! in interpolating gfdb interpolate all traces in block containing (ix,iz)
        
        integer :: iblockx, ixfirst, ixlast, ix, nblockt, ig, ilastrealx, ix_get
        integer :: iblockz, izfirst, izlast, iz, ilastrealz, iz_get
        integer :: inextrealx, inextrealz
        type(t_trace), pointer :: tracep
        integer, dimension(2) :: span, dataspan
        real, dimension(:,:,:), allocatable :: field_orig, field_interpol
        integer, dimension(:,:,:), allocatable :: spans
      
        allocate( spans(2,db%nblockz,db%nblockx) )
        tracep => null()
        
      ! in which block are we?
        iblockx = (ix_in-1)/db%nblockx_payload + 1
        iblockz = (iz_in-1)/db%nblockz_payload + 1
        
      ! index of first and last trace in block
        ixfirst = (iblockx-1)*db%nblockx_payload + 1 - db%nblockx_overlap/2
        izfirst = (iblockz-1)*db%nblockz_payload + 1 - db%nblockz_overlap/2
        ixlast = ixfirst+db%nblockx-1
        izlast = izfirst+db%nblockz-1
        
        call warn("interpolating: blockx=" // iblockx // ", blockz=" // iblockz  )
        
      ! load all real traces in block into cache & det. required extent 
        span = (/ huge(span(1)), -huge(span(2)) /)
        do ix=ixfirst,ixlast,db%nipx
            do iz=izfirst,izlast,db%nipz
                do ig=1,db%ng
                  ! repeat end points 
                    ix_get = (min(max(ix,1),db%nx)-1)/db%nipx *db%nipx+1
                    iz_get = (min(max(iz,1),db%nz)-1)/db%nipz *db%nipz+1
                    call gfdb_get_trace( db, ix_get, iz_get, ig, tracep )
                    if (.not. associated(tracep)) then
                        call warn("gfdb: gfdb_interpolate_block(): Missing trace in interpolation block " &
                                //iblockx//". Interpolation result may be inaccurate.")
                        spans(:,iz-izfirst+1,ix-ixfirst+1) = (/0,0/) 
                        cycle
                    end if
                    span(1) = min(tracep%span(1), span(1))
                    span(2) = max(tracep%span(2), span(2))
                    spans(:,iz-izfirst+1,ix-ixfirst+1) = tracep%span(:)
                end do
            end do
        end do
        span = allowed_span( span, min(64,int((span(2)-span(1))*1.2)) )
        
        nblockt = span(2)-span(1)+1
        if (nblockt <= 1) return
       
      ! 2d in and out arrays for the interpolator
        allocate( field_orig(nblockt, db%nblockz/db%nipz, db%nblockx/db%nipx) )
        allocate( field_interpol(nblockt, db%nblockz, db%nblockx) )
        
        do ig=1,db%ng
            field_orig(:,:,:) = 0.
            field_interpol(:,:,:) = 0.
          
          ! copy data into the work array
            do iz=izfirst,izlast,db%nipz
                iblockz = (iz-izfirst) + 1
                do ix=ixfirst,ixlast,db%nipx
                    iblockx = (ix-ixfirst) + 1
                    ix_get = (min(max(ix,1),db%nx)-1)/db%nipx *db%nipx+1
                    iz_get = (min(max(iz,1),db%nz)-1)/db%nipz *db%nipz+1
                    call gfdb_get_trace( db, ix_get, iz_get, ig, tracep )
                    call trace_multiply_add_nogrow( tracep, field_orig(:,(iblockz-1)/db%nipz+1,(iblockx-1)/db%nipx+1), span )
                end do
            end do
            
          ! interpolate
            call interpolate3d( field_orig, field_interpol, int(0.1*(span(2)-span(1))), db%nblockx_overlap/2, db%nblockz_overlap/2 )
                
          ! put interpolated data into gfdb traces
            do iz=izfirst+db%nblockz_overlap/2,izlast-db%nblockz_overlap/2
                iblockz = (iz-izfirst) + 1
                do ix=ixfirst+db%nblockx_overlap/2,ixlast-db%nblockx_overlap/2
                    iblockx = (ix-ixfirst) + 1
                    if (mod(ix-1, db%nipx)==0 .and. &
                        mod(iz-1, db%nipz)==0) cycle ! insert only the interpolated traces... 
                        
                    if (ix < 1 .or. db%nx < ix .or. &
                        iz < 1 .or. db%nz < iz) cycle
                    
                  ! set span of interpolated trace to union of neighboring real traces
                    ilastrealx = ((iblockx-1)/db%nipx)*db%nipx+1
                    ilastrealz = ((iblockz-1)/db%nipz)*db%nipz+1
                    
                    inextrealx = ilastrealx + db%nipx
                    inextrealz = ilastrealz + db%nipz
                    
                    dataspan(1) = spans(1,ilastrealz,ilastrealx)
                    if (inextrealx <= size(spans,3)) &
                        dataspan(1) = min(spans(1,ilastrealz,inextrealx), dataspan(1) )
                    if (inextrealz <= size(spans,2)) &
                        dataspan(1) = min(spans(1,inextrealz,ilastrealx), dataspan(1) )
                    if (inextrealx <= size(spans,3) .and. inextrealz <= size(spans,2)) &
                        dataspan(1) = min( spans(1,inextrealz,inextrealx), dataspan(1) )
                    
                    dataspan(2) = spans(2,ilastrealz,ilastrealx)
                    if (inextrealx <= size(spans,3)) &
                        dataspan(2) = max(spans(2,ilastrealz,inextrealx), dataspan(2) )
                    if (inextrealz <= size(spans,2)) &
                        dataspan(2) = max(spans(2,inextrealz,ilastrealx), dataspan(2) )
                    if (inextrealx <= size(spans,3) .and. inextrealz <= size(spans,2)) &
                        dataspan(2) = max( spans(2,inextrealz,inextrealx), dataspan(2) )
                    call gfdb_set_trace( db, ix, iz, ig, &
                        field_interpol(dataspan(1)-span(1)+1:dataspan(2)-span(1)+1,iblockz,iblockx), &
                        dataspan  )
                end do
            end do
            
        end do
        
        deallocate(field_orig)
        deallocate(field_interpol)
       
        deallocate(spans)   
       
    end subroutine
    
    subroutine interpolate3d( fin, fout, ntmargin, nxmargin, nzmargin )
    
        real, dimension(:,:,:), intent(inout)  :: fin  ! input field,  size: (nt, nz_in, nx_in)
        real, dimension(:,:,:), intent(out) :: fout ! output field, size: (nt, nz_out, nx_out)
        integer, intent(in) :: ntmargin, nxmargin, nzmargin
        
      ! Interpolate 3D array fin into 3D array fout
      
      ! Sizes of the arrays are powers of 2.
      ! Margins of size ntmargin in the first dimension and nxmargin in the second, 
      ! may be safely cluttered with junk.
        
        integer :: nx_in, nx_out, ix_in, ix_out, nipx
        integer :: nz_in, nz_out, iz_in, iz_out, nipz
        real, dimension(size(fin,1), size(fin,3))   :: inslice_x
        real, dimension(size(fout,1),size(fout,3))  :: outslice_x
        real, dimension(size(fin,1), size(fin,2))   :: inslice_z
        real, dimension(size(fout,1),size(fout,2))  :: outslice_z
        real, allocatable, dimension(:,:,:) :: finter
        
        nz_in = size(fin,2)
        nz_out = size(fout,2)
        nx_in = size(fin,3)
        nx_out = size(fout,3)
                
        nipx = nx_out/nx_in 
        nipz = nz_out/nz_in
       
        if ( nipz == 1 ) then
            call gulunay2d(fin(:,1,:), fout(:,1,:), ntmargin, nxmargin )
            return
        end if
        
        if (nipx == 1 ) then
            call gulunay2d(fin(:,:,1), fout(:,:,1), ntmargin, nzmargin)
            return

        end if

        if (nipx == 4 .and. nipz == 4) then
            allocate(finter(size(fout,1),size(fout,2)/2,size(fout,3)/2))
            call gulunay3d(fin,finter, ntmargin, nxmargin/2, nzmargin/2)
            call gulunay3d(finter,fout, ntmargin, nxmargin, nzmargin)
            deallocate( finter )
            return
        end if
      
     ! 3D interpolation if possible
       if (nipx /= 1 .and. nipz /= 1 .and. nipx == nipz ) then
            call gulunay3d(fin, fout, ntmargin, nxmargin, nzmargin)
            return
       end if        
      
      ! until then: pseudo 3D interpolation:
      ! horizontal interpolation
        do iz_in=1,nz_in
            iz_out=(iz_in-1)*nipz+1
            inslice_x(:,:) = fin(:,iz_in,:)
            call gulunay2d(inslice_x, outslice_x, ntmargin, nxmargin )
            fout(:,iz_out,:) = outslice_x(:,:)
        end do
        
      ! vertical interpolation
        do ix_out=1,nx_out
            ix_in=(ix_out-1)/nipx+1
            if (mod(ix_in-1,nipx)==0) then ! on rows available in input, use original input
                inslice_z(:,:) = fin(:,:,ix_in)
            else ! use horizontally interpolated
                inslice_z(:,:) = fout(:,::nipz,ix_out)
            end if
            call gulunay2d(inslice_z, outslice_z, ntmargin, nxmargin )
            fout(:,:,ix_out) = outslice_z(:,:)
        end do
        
    end subroutine
    
   
    pure function allowed_span( span, minlength ) result(newspan)
    
        integer, dimension(2), intent(in) :: span
        integer, intent(in)               :: minlength
        integer, dimension(2)             :: newspan
        
        integer :: length, lengthp
        
        newspan(:) = span(:)
        length = slen(newspan)
        if (length < minlength) then
            length = minlength
        end if
        lengthp = next_power_of_two( length )
        newspan(1) = newspan(1) - floor((lengthp-length)/2.)
        newspan(2) = newspan(1) + lengthp - 1
        
    end function
    
    pure function next_power_of_two(n) result(m)
    
        integer, intent(in) :: n
        integer :: m
        
        m = 2**ceiling(log(real(n))/log(2.))
    
    end function
    
    subroutine gfdb_uncache_trace( db, ix,iz,ig )
    
      ! remove a cached trace from the trace cache
        
        type(t_gfdb), intent(inout)            :: db
        integer, intent(in)                   :: ix, iz, ig
        
        integer :: ichunk, ixc
        
        if (.not. allocated(db%chunks) .or. &
            ix > db%nx .or. ix < 1 .or. &
            iz > db%nz .or. iz < 1 .or. &
            ig > db%ng .or. ig < 1) then
            call warn("gfdb: gfdb_uncache_trace(): invalid request: out of bounds: " // &
                      "("//ix//","//iz//","//ig//")")
            return
        end if
        call gfdb_index_to_chunk( db, ix, ichunk, ixc )
        call chunk_uncache_trace( db, db%chunks(ichunk), ixc,iz,ig )
    
    end subroutine
    
    pure function slen( span )
        integer, dimension(2), intent(in) :: span
        integer :: slen
        slen = span(2) - span(1) + 1
    end function
    
    subroutine chunk_uncache_trace(db, c, ixc,iz,ig )
    
        type(t_gfdb), intent(inout)           :: db
        type(t_chunk), intent(inout)          :: c
        integer, intent(in)                   :: ixc, iz,ig
    
        integer :: nbytes

        if (ixc > c%nxc .or. ixc < 1 .or. &
            iz > c%nz .or. iz < 1 .or. &
            ig > c%ng .or. ig < 1) then
            call warn("gfdb: chunk_uncache_trace(): invalid request: out of bounds : "//&   
                      "("//ixc//","//iz//","//ig//")")
            return
        end if
                
        call chunk_open_read( db, c ) ! ensure that chunk is open for read
        
        if (associated(c%traces(ig,iz,ixc)%p)) then
             nbytes = trace_size_bytes( c%traces(ig,iz,ixc)%p )
             c%nbytes = c%nbytes - nbytes
             db%nbytes = db%nbytes - nbytes

             call trace_destroy( c%traces(ig,iz,ixc)%p )
             deallocate( c%traces(ig,iz,ixc)%p )
             c%traces(ig,iz,ixc)%p => null()
        end if
        
        

    end subroutine
    
    subroutine gfdb_index_to_chunk( c, ix, ichunk, ixc )
    
      ! traslate ix to ichunk and ixc
     
        type(t_gfdb), intent(in) :: c
        integer, intent(in)         :: ix
        integer, intent(out)        :: ichunk, ixc
        
        ichunk = (ix-1) / c%nxc + 1
        if (ichunk > c%nchunks) ichunk = c%nchunks
        ixc = ix - (ichunk-1)*c%nxc
    
    end subroutine
    
    subroutine gfdb_close( db )
        
        type(t_gfdb), intent(inout) :: db
        integer :: ichunk, iipt

        if (.not. allocated(db%chunks)) return
        
        do ichunk=1,db%nchunks     
            call chunk_close(db, db%chunks(ichunk))
        end do
        if (allocated(db%interpolated_traces)) then
            do iipt=1,size(db%interpolated_traces)
                if (associated( db%interpolated_traces(iipt)%p )) then
                    call trace_destroy( db%interpolated_traces(iipt)%p )
                    deallocate( db%interpolated_traces(iipt)%p )
                end if
            end do
            deallocate( db%interpolated_traces )
        end if
       
    end subroutine
    
    subroutine chunk_open_read( db, c )
    
      ! open chunk for read
      ! this opens the associated database file
      ! and reads the index dataset.
      ! it also initializes the traces cache
        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c
        
        integer :: ig,iz,ixc
        logical :: ok
        
        if ( c%readmode ) return
        if ( c%writemode ) call chunk_close( db, c )

        call gfdb_io_chunk_open_read( c%filename, c%file, ok )
        if (.not. ok) call die()
        
        allocate( c%traces(c%ng, c%nz, c%nxc) )        
        forall (ig=1:c%ng, iz=1:c%nz, ixc=1:c%nxc) 
            c%traces(ig,iz,ixc)%p => null()
        end forall

        allocate( c%references(c%ng, c%nz/c%nipz, c%nxc/c%nipx) )

        call gfdb_io_chunk_read_index( c%file, c%dataset_index, c%references, ok )
        if (.not. ok) call die()
        
        c%readmode = .true.

    end subroutine 
    
    subroutine chunk_open_write(db,c)
    
      ! open chunk for write
      ! this opens the associated database file 
      ! and the index dataset.

        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c

        logical :: ok

        
        if ( c%writemode ) return
        if ( c%readmode ) call chunk_close(db,c)
     
        call gfdb_io_chunk_open_write( c%filename, c%file, ok )
        if (.not. ok) call die()
                                    
      ! read references to detect duplicate inserts
        allocate( c%references(c%ng, c%nz/c%nipz, c%nxc/c%nipx) )
   
        call gfdb_io_chunk_read_index( c%file, c%dataset_index, c%references, ok )
        if (.not. ok) call die()
        
        c%writemode = .true.

    end subroutine
    
    subroutine chunk_close( db, c )
        type(t_gfdb), intent(inout) :: db
        type(t_chunk), intent(inout) :: c   
        integer(hid_t) :: err
        integer(hid_t), dimension(2) :: e
        integer :: ixc,iz,ig
        integer :: nbytes
        logical :: ok
        
        if (c%file /= 0) then
            call gfdb_io_chunk_close(c%file, c%dataset_index, ok)
            if (.not. ok) call die()
        end if
        
        c%file = 0
        c%dataset_index = 0
        c%writemode = .false.
        c%readmode = .false.
        
        if (allocated(c%references)) deallocate( c%references )
        if (allocated(c%traces)) then
            do ixc=1,c%nxc
                do iz=1,c%nz
                    do ig=1,c%ng
                        if (associated( c%traces(ig,iz,ixc)%p )) then
                            nbytes = trace_size_bytes( c%traces(ig,iz,ixc)%p )
                            c%nbytes = c%nbytes - nbytes
                            db%nbytes = db%nbytes - nbytes

                            call trace_destroy( c%traces(ig,iz,ixc)%p )
                            deallocate( c%traces(ig,iz,ixc)%p )
                        end if
                    end do
                end do
            end do
            deallocate( c%traces )
        end if
        
        call h5garbage_collect_f(err)

        call gfdb_accesslist_remove( db, c%ichunk )

    end subroutine
    


end module




