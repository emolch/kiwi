! $Id: gfdb.f90 702 2008-03-31 07:07:40Z sebastian $ 
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

module gfdb

    use util
    use sparse_trace
    use better_varying_string
    use hdf5
    use interpolation
    
    implicit none

    public
    
    logical :: h5inited = .false.
    
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
        integer :: ng = 0  ! number of greens functions (or elementary seismograms) (=8)
        
        integer :: nipx = 1       ! total traces / real traces
                                  ! (increment between 2 non-interpolated traces)
        integer :: nipz = 1       ! total traces / real traces
                                  ! (increment between 2 non-interpolated traces)
        
        logical :: writemode = .false.
        logical :: readmode = .false.
        type(varying_string) :: filename
        integer(hid_t) :: file = 0
        integer(hid_t) :: dataset_index = 0
        
      ! references to the datasets are stored in this array for quick access:
        type(hobj_ref_t_f), dimension(:,:,:), allocatable :: references    ! (ng,nz/nipz,nxc/nipx)
      
      ! pointers to traces, where allready used traces live
        type(t_trace_p), dimension(:,:,:), allocatable :: traces           ! (ng,nz,nxc)
  
      ! these are buffers, needed by chunk_get_trace,
      ! they have been put here, so that their allocation can be reused
        integer, dimension(:), allocatable :: ofs, pofs
        real, dimension(:), allocatable :: packed
        
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
        integer :: ng  = 0     ! number of greens functions (= 8)
        
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
        
        
      ! cache contains nchunks chunks:
        integer :: nchunks = 0
        
        type(t_chunk), dimension(:), allocatable :: chunks
        
        type(t_trace), pointer :: interpolated_trace => null()

    end type
    
    interface h5_save_scalar
        
        subroutine h5_save_scalar_integer( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            integer,intent(in) :: value
            integer(hid_t), intent(out) :: error

        end subroutine
        
        subroutine h5_save_scalar_real( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            real,intent(in) :: value
            integer(hid_t), intent(out) :: error

        end subroutine
        
    end interface
    
    interface h5_open_scalar
        
        subroutine h5_open_scalar_integer( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            integer,intent(out) :: value
            integer(hid_t), intent(out) :: error

        end subroutine
        
        subroutine h5_open_scalar_real( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            real,intent(out) :: value
            integer(hid_t), intent(out) :: error

        end subroutine
        
    end interface
    
    public gfdb_init, gfdb_save_trace, gfdb_close, gfdb_get_trace, gfdb_uncache_trace
    public gfdb_get_indices, gfdb_infos, gfdb_dump_infomap, gfdb_dump_contents
    public gfdb_dump_missing
    
  contains
  
    subroutine gfdb_cleanup()
        
        integer(hid_t) :: egal
    
        if (h5inited) then
            call h5close_f(egal)
            !call h5close_types(egal)
            h5inited = .false.
        end if
        
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
        integer(hid_t) :: file, e(10)
        integer(hid_t) :: egal
        
        if (.not. h5inited) then
            !call h5init_types_( e(1) )
            call h5open_f( e(1) )
            h5inited = .true.
        end if
        
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
        
            e = 0
        
            c%filenamebase = fnbase
            filename = fnbase // ".index"
            c%nipx = 1
            c%nipz = 1
            
            if (present(nipx)) c%nipx = nipx
            if (present(nipz)) c%nipz = nipz
            
            call h5eset_auto_f(0,egal)
            call h5fopen_f(char(filename), H5F_ACC_RDONLY_F, file, e(1))
            call h5eset_auto_f(1,egal)
            if (e(1) /= 0) call die( "gfdb: failed to open file: " //filename )
            
          ! read dt,dx,dz
            call h5_open_scalar( file, "dt", c%dt, e(1) )
            call h5_open_scalar( file, "dx", c%dx, e(2) )
            call h5_open_scalar( file, "dz", c%dz, e(3) )
            if (any(e/=0)) &
                call die( "gfdb: failed to read dataset from file: "//filename )
        
          ! try to read firstx, firstz; leave these at zero when not found 
          ! (for backward compatibility)            
            call h5eset_auto_f(0,egal)
            call h5_open_scalar( file, "firstx", c%firstx, e(4) )
            call h5_open_scalar( file, "firstz", c%firstz, e(5) )
            call h5eset_auto_f(1,egal)
            call h5eclear_f(egal)            
            e=0
            
          ! save nx, nxc, nz, ng
            call h5_open_scalar( file, "nchunks", c%nchunks, e(1) )
            call h5_open_scalar( file, "nx", c%nx, e(2) )
            call h5_open_scalar( file, "nxc", c%nxc, e(3) )
            call h5_open_scalar( file, "nz", c%nz, e(4) )
            call h5_open_scalar( file, "ng", c%ng, e(5) )
            if (any(e/=0)) &
                call die( "gfdb: failed to read dataset from file: "//filename )

            call h5fclose_f( file, e(1) )
            if (any(e/=0)) &
                call die( "gfdb: problems closing file: "//filename )
            
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
            call chunk_infos( c%chunks(ichunk), ntraces, ntraces_used )
            ntraces_in_chunks(1,ichunk) = ntraces_used
            ntraces_in_chunks(2,ichunk) = ntraces
        end do
        
    end subroutine
    
    subroutine gfdb_dump_infomap( c, iunit )
        type(t_gfdb), intent(inout) :: c
        integer, intent(in) :: iunit
        integer :: ichunk
        
        do ichunk=1,c%nchunks
            call chunk_dump_infomap( c%chunks(ichunk), iunit,c%nxc/c%nipx  )
        end do
        
    end subroutine

    subroutine chunk_infos( c, ntraces, ntraces_used )
    
        type(t_chunk), intent(inout) :: c
        integer, intent(out) :: ntraces, ntraces_used
        integer :: ixc, iz, ig

        ntraces = c%nxc/c%nipx*c%nz/c%nipz*c%ng
        
        call chunk_open_read( c ) 
        
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
        call chunk_close( c )
        
    end subroutine
    
    
    subroutine chunk_dump_infomap( c, iunit, nxcfull )
    
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: iunit, nxcfull
        integer :: ixc, iz, ig
        integer :: ialloc
        
        call chunk_open_read( c ) 
        
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
        
        call chunk_close( c )
    
    end subroutine
    
    subroutine gfdb_dump_contents( c, iunit )
        type(t_gfdb), intent(inout) :: c
        integer, intent(in) :: iunit
        integer :: ichunk
        
        do ichunk=1,c%nchunks
            call chunk_dump_contents( c%chunks(ichunk), c, iunit,c%nxc/c%nipx  )
        end do
        
    end subroutine
    
    subroutine chunk_dump_contents( c, db, iunit, nxcfull )
    
        type(t_chunk), intent(inout) :: c
        type(t_gfdb), intent(inout) :: db
        integer, intent(in) :: iunit, nxcfull
        integer :: ixc,ix, iz, ig
        real :: x,z
        
        call chunk_open_read( c ) 
        
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
        
        call chunk_close( c )
    
    end subroutine
    
    subroutine gfdb_dump_missing( c, iunit )
        type(t_gfdb), intent(inout) :: c
        integer, intent(in) :: iunit
        integer :: ichunk
        
        do ichunk=1,c%nchunks
            call chunk_dump_missing( c%chunks(ichunk), c, iunit,c%nxc/c%nipx  )
        end do
        
    end subroutine
    
    subroutine chunk_dump_missing( c, db, iunit, nxcfull )
    
        type(t_chunk), intent(inout) :: c
        type(t_gfdb), intent(inout) :: db
        integer, intent(in) :: iunit, nxcfull
        integer :: ixc,ix, iz, ig
        real :: x,z
        
        call chunk_open_read( c ) 
        
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
        
        call chunk_close( c )
    
    end subroutine


    subroutine gfdb_destroy( c )
    
      ! destroy cache object c 
      ! free all memory associated with c
      ! reset everything to 0
      ! it should be safe to destroy a not initialized cache
         
        type(t_gfdb), intent(inout) :: c
        integer :: ichunk
        
        call gfdb_close( c )

        if (allocated(c%chunks)) then
            do ichunk=1,c%nchunks
                call chunk_destroy( c%chunks(ichunk) )
            end do
            deallocate( c%chunks )
        end if
        
        call delete(c%filenamebase)
        c%nchunks = 0
        c%nx = 0
        c%nxc = 0
        c%nz = 0
        c%ng = 0
        c%dt = 0.
        c%dx = 0.
        c%dz = 0.
        c%nipx = 1
        c%nipz = 1
        
    end subroutine
    
    subroutine chunk_destroy( c )
    
        type(t_chunk), intent(inout) :: c
        
        call chunk_close( c )
        c%nxc = 0
        c%nz = 0
        c%ng = 0
        call delete(c%filename)
        c%writemode = .false.
        c%readmode = .false.
        c%nipx = 1
        c%nipz = 1
        
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
        if (present(nipx)) c%nipx = nipx
        if (present(nipz)) c%nipz = nipz
        
        c%file = 0
        c%dataset_index = 0
        
    end subroutine
    
    subroutine gfdb_create( c )
    
        type(t_gfdb), intent(inout) :: c
        
        integer :: ichunk
        integer(hsize_t),dimension(1) :: dims
        integer(hid_t) :: file
        integer(hid_t) :: egal
        integer(hid_t),dimension(7) :: e
        type(varying_string) :: filename
        
        if (.not. allocated(c%chunks)) return

        e = 0
        dims = (/1/)
        
        filename = c%filenamebase // ".index"
        call h5eset_auto_f(0,egal)
        call h5fcreate_f(char(filename),H5F_ACC_TRUNC_F,  file,e(1))
        call h5eset_auto_f(1,egal)
        if (e(1) /= 0) call die( "gfdb: failed to create file: "//filename )
        
      ! save dt,dx,dz
        call h5_save_scalar( file, "dt", c%dt, e(1) )
        call h5_save_scalar( file, "dx", c%dx*c%nipx, e(2) )
        call h5_save_scalar( file, "dz", c%dz*c%nipz, e(3) )
        call h5_save_scalar( file, "firstx", c%firstx, e(4) )
        call h5_save_scalar( file, "firstz", c%firstz, e(5) )
        if (any(e/=0)) &
            call die( "gfdb: failed to write dataset to file: "//filename )

        
      ! save nx, nxc, nz, ng
        call h5_save_scalar( file, "nchunks", c%nchunks, e(1) )
        call h5_save_scalar( file, "nx", c%nx/c%nipx, e(2) )
        call h5_save_scalar( file, "nxc", c%nxc/c%nipx, e(3) )
        call h5_save_scalar( file, "nz", c%nz/c%nipz, e(4) )
        call h5_save_scalar( file, "ng", c%ng, e(5) )
        if (any(e/=0)) &
            call die( "gfdb: failed to write dataset to file: "//filename )
        
        call h5fclose_f(file,  e(1))
        if (any(e/=0)) &
                call die( "gfdb: problems closing file: "//filename )
                
        do ichunk=1,c%nchunks
            call chunk_create( c%chunks(ichunk))
        end do
    
    end subroutine
    
    subroutine chunk_create( c )
    
        type(t_chunk), intent(inout) :: c
        
        type(hobj_ref_t_f), dimension(:), allocatable :: refs
        integer(hid_t) :: error, file, dataspace, dataset, group, egal
        integer(hsize_t), dimension(3) :: dims
        integer(hsize_t), dimension(1) :: dims1
        integer :: i, length
        
        dims(1) = c%ng
        dims(2) = c%nz/c%nipz
        dims(3) = c%nxc/c%nipx
        call h5eset_auto_f(0,egal)
        call h5fcreate_f(char(c%filename),H5F_ACC_TRUNC_F,  file,error)
        call h5eset_auto_f(1,egal)
        if (error /= 0) call die( "gfdb: failed to create file: "//c%filename )

      ! create index dataset for faster access to individual traces
        call h5screate_simple_f(3,dims,  dataspace,error)
        if (error /= 0) &
            call die( "gfdb: failed to write dataset to file: "//c%filename )
      
        call h5dcreate_f(file,"index",H5T_STD_REF_OBJ,dataspace,  dataset,error)
        if (error /= 0) &
            call die( "gfdb: failed to create index dataset in file: "//c%filename )

        dims1(1) = dims(1)*dims(2)*dims(3)
        length = int(dims1(1)) ! remove a compiler warning 
        allocate(refs(length))
        
        do i=1,dims1(1)
            refs(i)%ref = 0
        end do
        
        call h5dwrite_f(dataset, H5T_STD_REF_OBJ, refs, dims1, error )
        if (error /= 0) &
            call die( "gfdb: failed to write index dataset to file: "//c%filename )

        deallocate(refs)
            
        call h5gcreate_f(file,"gf",group,error,int(c%nxc/c%nipx*(floor(log10(real(c%nxc/c%nipx)))+2)))
        if (error /= 0) &
            call die( "gfdb: failed to create group gf in file: "//c%filename )
        
        call h5gclose_f(group, error)    
        
        call h5dclose_f(dataset,  error)
        call h5sclose_f(dataspace,  error)        
        call h5fclose_f(file,  error)
        if (error/=0) &
             call die( "gfdb: problems closing file: "//c%filename )
             
    end subroutine
    
    subroutine gfdb_save_trace( c, ix, iz, ig, data )
    
      ! save a trace to the database
        
        type(t_gfdb), intent(inout) :: c
        integer, intent(in) :: ix,iz,ig
        type(t_trace), intent(in) :: data
        
        integer :: ichunk, ixc
        
        if (.not. allocated(c%chunks)) return
                        
        call gfdb_index_to_chunk( c, ix, ichunk, ixc )
        
        call chunk_save_trace( c%chunks(ichunk), ixc,iz,ig, data )
        
    end subroutine
    
    subroutine chunk_save_trace( c, ixc, iz, ig, data )
    
      ! save a trace to this chunks file
    
        type(t_chunk), intent(inout) :: c
        integer, intent(in) :: ixc,iz,ig
        type(t_trace), intent(in) :: data
        
        integer(hid_t), dimension(8) :: e
        integer(hid_t) :: dataset, dataspace, dataspace1, attribute, group_dist
        integer(hid_t) :: dataspace_for_ref, memspace, group
        integer(hsize_t),dimension(1) :: dims
        
        integer :: totalsize, nstrips
        type(hobj_ref_t_f), dimension(1) :: reference
        type(varying_string) :: datasetname, gloc
        integer(hsize_t), dimension(3,1) :: coord
        integer :: ixc_file, iz_file
        
        e=0
        
        if (ixc > c%nxc .or. ixc < 1) return
        if (iz > c%nz .or. iz < 1) return
        if (ig > c%ng .or. ig < 1) return
        
      ! cannot save interpolated trace this way.
        if (mod(ixc-1,c%nipx) /= 0) return
        if (mod(iz-1,c%nipz) /= 0) return
        ixc_file = (ixc-1)/c%nipx + 1
        iz_file = (iz-1)/c%nipz + 1
        
        if (trace_is_empty(data)) return
        
        call chunk_open_write( c ) ! ensure chunk is open for write
        
      ! check that there is no trace already stored at that position
        if (c%references(ig,iz_file,ixc_file)%ref /= 0) then
            call warn( "gfdb: chunk_save_trace() (chunk "//c%ichunk// &
                       "): trace already exists for index "//&   
                       "("//ixc_file//","//iz_file//","//ig//")" )
            return
        end if
        e=0
        
      ! open/create group for the dataset
                
        gloc = var_str("/gf/") // ixc_file
        call h5_opencreategroup( c%file, char(gloc), group_dist, e(1), int(c%nz/c%nipz*(floor(log10(real(c%nz/c%nipz)))+2)) )
        gloc = iz_file
        call h5_opencreategroup( group_dist, char(gloc), group, e(2), int(8*(floor(log10(8.))+2)) )
        if (any(e /= 0)) call die( "gfdb: failed to open/create group "//gloc//" in file: " &
                                      // c%filename )
                
      ! pack data to a single, continuous array
      ! generate arrays with offsets in the packed array (pofs)
      ! and offsets of the data strips (ofs)
      
        call trace_to_storable( data, c%packed, totalsize, c%pofs, c%ofs, nstrips )
        
      ! create dataset for the array
        
        dims(1) = totalsize
        call h5screate_simple_f(1,dims,  dataspace,e(1))
        datasetname = ig
        call h5dcreate_f(group,char(datasetname), &
                         H5T_NATIVE_REAL,dataspace,  dataset,e(2))
        if (any(e /= 0)) call die( "gfdb: failed to create dataset "// &
                                    datasetname // " in file: " // c%filename )
       
      ! save offsets as attributes to the dataset  
                         
        dims(1) = nstrips
        call h5screate_simple_f(1,dims, dataspace1, e(1))
        
        call h5acreate_f(dataset,"pofs",H5T_NATIVE_INTEGER,dataspace1, &
                         attribute,e(2))
        call h5awrite_f(attribute,H5T_NATIVE_INTEGER,c%pofs,dims,e(3))
        call h5aclose_f(attribute, e(4))
        
        call h5acreate_f(dataset,"ofs",H5T_NATIVE_INTEGER,dataspace1, &
                         attribute,e(5))
        call h5awrite_f(attribute,H5T_NATIVE_INTEGER,c%ofs,dims,e(6))
        call h5aclose_f(attribute, e(7))
        call h5sclose_f(dataspace1, e(8))
        if (any(e /= 0)) call die( "gfdb: failed to save attributes to dataset "// &
                                    datasetname // " in file: " // c%filename )
       
      ! save the array 
      
        dims(1) = totalsize
        call h5dwrite_f(dataset, H5T_NATIVE_REAL, c%packed, dims, e(1))
        if (any(e /= 0)) call die( "gfdb: failed to write dataset "// &
                                    datasetname // " in file: " // c%filename )
      
      ! save a reference to the dataset in the index block  
        
        call h5rcreate_f(group, char(datasetname), reference(1), e(1))
        coord(:,1) = (/ ig, iz_file, ixc_file /)
        
        dims(1) = 1
        call h5dget_space_f(c%dataset_index,dataspace_for_ref,e(2))
        call h5sselect_elements_f(dataspace_for_ref,H5S_SELECT_SET_F,3, 1, coord, e(3))
        call h5screate_simple_f(1, dims, memspace, e(4))
        
        call h5dwrite_f(c%dataset_index, H5T_STD_REF_OBJ, reference, dims, e(5), &
                         mem_space_id=memspace, file_space_id=dataspace_for_ref)
                         
        call h5sclose_f(memspace, e(6))
        call h5sclose_f(dataspace_for_ref, e(7))
        if (any(e /= 0)) call die( "gfdb: failed save reference to "// &
                                    datasetname // " in index dataset of file: " &
                                    // c%filename )
                                    
        c%references(ig,iz_file,ixc_file) = reference(1)
        
        call h5dclose_f(dataset, e(1))
        call h5sclose_f(dataspace, e(1))
        call h5gclose_f(group_dist, e(1))
        call h5gclose_f(group, e(1))

        
    end subroutine
    
    subroutine gfdb_get_indices( c, x, z, ix, iz )
    
      ! make indices suitable for get_trace, based on spacing of gf's
      
        type(t_gfdb), intent(in) :: c
        real, intent(in) :: x, z
        integer, intent(out) :: ix, iz
        
        ix = nint((x-c%firstx)/c%dx)+1
        iz = nint((z-c%firstz)/c%dz)+1
        
    end subroutine
    
    subroutine gfdb_get_indices_bilin( c, x, z, ix, iz, dix, diz )
    
      ! make indices suitable for get_trace, based on spacing of gf's
      ! additionally return fractional indices for use with bilinear interpolation
      ! these are in the range of [-0.5,0.5]
      
        type(t_gfdb), intent(in) :: c
        real, intent(in) :: x, z
        integer, intent(out) :: ix, iz
        real, intent(out) :: dix, diz
        
        ix = nint((x-c%firstx)/c%dx)+1
        iz = nint((z-c%firstz)/c%dz)+1
    
        dix = (x-c%firstx)/c%dx - (ix-1)
        diz = (z-c%firstz)/c%dz - (iz-1)
        
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
        
        call chunk_get_trace( db, db%chunks(ichunk), ixc,iz,ig, tracep )
    
    end subroutine
    
    subroutine gfdb_get_trace_bilin( db, ix, iz, ig, dix, diz, tracep )
    
        type(t_gfdb), intent(inout)           :: db
        integer, intent(in)                   :: ix, iz, ig
        real, intent(in)                      :: dix, diz
        type(t_trace), pointer   :: tracep
    
      ! use bilinear interpolation to get a trace between gfdb grid nodes.
      ! ix, iz, dix and diz are indices and offsets as produced by gfdb_get_indices_bilin()
      ! 
      ! if dix==0 and diz==0 
      
        
        type(t_trace), pointer  :: t00, t01, t10, t11
        integer                 :: ixl, izl
        real                    :: dixl, dizl
        integer, dimension(2)   :: span
        tracep => null()
        t00 => null()
        t01 => null()
        t10 => null()
        t11 => null()
        
      ! if we are exactly at a grid node, no interpolation is needed...
        if (dix .eq. 0. .and. diz .eq. 0.) then
            call gfdb_get_trace( db, ix, iz, ig, tracep )
        end if
        
      ! select appropriate grid square
        ixl = ix
        dixl = dix
        if (dixl < 0) then
            ixl = ixl - 1
            dixl = 1 + dixl
        end if
        
        izl = iz
        dizl = diz
        if (dizl < 0) then
            izl = izl - 1
            dizl = 1 + dizl
        end if
        
      ! get the traces
        call gfdb_get_trace( db, ixl,   izl,   ig, t00 )
        call gfdb_get_trace( db, ixl,   izl+1, ig, t01 )
        call gfdb_get_trace( db, ixl+1, izl,   ig, t10 )
        call gfdb_get_trace( db, ixl+1, izl+1, ig, t11 )
        
        if (.not. (associated(t00) .and. associated(t01) .and. &
                   associated(t10) .and. associated(t11))) then
            tracep => null()
            return
        end if


      ! get union of spans
        span(1) = min( t00%span(1), t01%span(1), t10%span(1), t11%span(1) )
        span(2) = max( t00%span(2), t01%span(2), t10%span(2), t11%span(2) )
      
      ! prepare the global buffer trace which will be reused on first invocation
        if (.not. associated( db%interpolated_trace ) ) then
            allocate( db%interpolated_trace )
        end if

        if ( trace_is_empty( db%interpolated_trace ) ) then
            call trace_create_simple_nodata( db%interpolated_trace, span )
        else ! resize the contained buffer
            !!! relying on internals of trace_t here...
            call resize( db%interpolated_trace%strips(1)%data, span(1), span(2)-span(1)+1 ) !!! relying on internals of trace_t ....
            db%interpolated_trace%span(:) = span(:)
        end if
        
      ! summation
        db%interpolated_trace%strips(1)%data(:) = 0
        
        call trace_multiply_add_nogrow( t00, db%interpolated_trace%strips(1)%data, span, (1-dixl)*(1-dizl) )
        call trace_multiply_add_nogrow( t01, db%interpolated_trace%strips(1)%data, span, (1-dixl)*dizl )
        call trace_multiply_add_nogrow( t10, db%interpolated_trace%strips(1)%data, span, dixl*(1-dizl) )
        call trace_multiply_add_nogrow( t11, db%interpolated_trace%strips(1)%data, span, dixl*dizl )
        
        tracep => db%interpolated_trace
        
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
        
        integer(hsize_t),dimension(1) :: adims
        integer(hid_t) :: dataset, attribute, space
        integer :: compactsize, nstrips
        integer(hid_t), dimension(8) :: e
        integer(hid_t) :: egal
        integer(hsize_t) :: length
        integer :: ixc_file, ix
        integer :: iz_file
        
        
        e=0
        
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
        call chunk_open_read( c ) 
        
     
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
        
        allocate(c%traces(ig,iz,ixc)%p)
        tracep => c%traces(ig,iz,ixc)%p
        
        call trace_destroy( tracep )
        
      ! lookup reference to the dataset for a quick jump
        call h5eset_auto_f(0,egal)
        call h5rdereference_f(c%dataset_index, c%references(ig,iz_file,ixc_file),&
                              dataset, e(1) )
        call h5eset_auto_f(1,egal)
        call h5eclear_f(egal)
        if (e(1) /= 0) then
            call warn( "gfdb: chunk_get_trace() (chunk "//c%ichunk// &
                       "): no trace available for index "//&   
                       "("//ixc_file//","//iz_file//","//ig//")" )
            return
        end if
        
      ! get offsets
                
        call h5aopen_idx_f(dataset, 0, attribute, e(1))
        call h5aget_space_f( attribute, space, e(2) )
        call h5sget_simple_extent_npoints_f(space, length, e(3) )
        call h5sclose_f( space, e(4) ) 
        nstrips = int(length)
        
        if (any(e /= 0)) call die( "gfdb: failed to get size of attributes of dataset with" // &
                                   "ixc="//ixc // ", iz="//iz //" and ig="//ig )
       
                                   
        if (.not. allocated(c%pofs)) &
            call resize( c%pofs, 1, nstrips )
       
        if (nstrips > size(c%pofs)) &
           call resize( c%pofs, 1, nstrips )

        if (.not. allocated(c%ofs)) &
            call resize( c%ofs, 1, nstrips )
            
        if (nstrips > size(c%ofs)) &
            call resize( c%ofs, 1, nstrips )

        adims(1) = nstrips
        call h5aread_f(attribute, H5T_NATIVE_INTEGER, c%pofs, adims, e(1))
        call h5aclose_f( attribute, e(2) )
        
        call h5aopen_idx_f(dataset, 1, attribute, e(3))
        call h5aread_f(attribute, H5T_NATIVE_INTEGER, c%ofs, adims, e(4))
        call h5aclose_f( attribute, e(5) )
        
        if (any(e /= 0)) call die( "gfdb: failed to get attributes of dataset with" // &
                                   "ixc="//ixc // ", iz="//iz //" and ig="//ig )
        
      ! get size of data
      
        call h5dget_space_f(dataset,space,e(1))
        call h5sget_simple_extent_npoints_f(space, length, e(2))
        call h5sclose_f( space, e(3) )
        
      ! get the data
  
        compactsize = int(length)
        if (.not. allocated(c%packed)) &
            call resize(c%packed, 1, compactsize)
            
        if (compactsize > size(c%packed)) &
            call resize(c%packed, 1, compactsize)

        adims(1) = compactsize
        call h5dread_f(dataset, H5T_NATIVE_REAL, c%packed, adims, e(4))
        call h5dclose_f(dataset, e(5))
        if (any(e /= 0)) call die( "gfdb: failed to get dataset with" // &
                                   "ixc="//ixc // ", iz="//iz //" and ig="//ig )
    
      ! convert to sparse trace
        call trace_from_storable( tracep, c%packed(1:compactsize), &
                                  c%pofs(1:nstrips), c%ofs(1:nstrips) )
    
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
        
        call chunk_set_trace( db%chunks(ichunk), ixc,iz,ig, data, span )
    
    end subroutine
    
    subroutine chunk_set_trace( c, ixc, iz, ig, data, span )
    
      ! insert an interpolated trace into this chunk
      ! this does not save the trace as gfdb_save_trace does.
            
        type(t_chunk), intent(inout)          :: c
        integer, intent(in)                   :: ixc, iz, ig
        real, dimension(:), intent(in)        :: data
        integer, dimension(:), intent(in)     :: span

        type(t_trace), pointer                :: tracep
        
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
        call chunk_open_read( c )
     
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
    
    subroutine gfdb_uncache_trace( c, ix,iz,ig )
    
      ! remove a cached trace from the trace cache
        
        type(t_gfdb), intent(inout)            :: c
        integer, intent(in)                   :: ix, iz, ig
        
        integer :: ichunk, ixc
        
        if (.not. allocated(c%chunks) .or. &
            ix > c%nx .or. ix < 1 .or. &
            iz > c%nz .or. iz < 1 .or. &
            ig > c%ng .or. ig < 1) then
            call warn("gfdb: gfdb_uncache_trace(): invalid request: out of bounds: " // &
                      "("//ix//","//iz//","//ig//")")
            return
        end if
        call gfdb_index_to_chunk( c, ix, ichunk, ixc )
        call chunk_uncache_trace( c%chunks(ichunk), ixc,iz,ig )
    
    end subroutine
    
    pure function slen( span )
        integer, dimension(2), intent(in) :: span
        integer :: slen
        slen = span(2) - span(1) + 1
    end function
    
    subroutine chunk_uncache_trace( c, ixc,iz,ig )
    
        type(t_chunk), intent(inout)          :: c
        integer, intent(in)                 :: ixc, iz,ig
    
        if (ixc > c%nxc .or. ixc < 1 .or. &
            iz > c%nz .or. iz < 1 .or. &
            ig > c%ng .or. ig < 1) then
            call warn("gfdb: chunk_uncache_trace(): invalid request: out of bounds : "//&   
                      "("//ixc//","//iz//","//ig//")")
            return
        end if
                
        call chunk_open_read( c ) ! ensure that chunk is open for read
        
        if (associated(c%traces(ig,iz,ixc)%p)) then       
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
    
    subroutine gfdb_close( c )
        
        type(t_gfdb), intent(inout) :: c
        integer :: ichunk

        if (.not. allocated(c%chunks)) return
        
        do ichunk=1,c%nchunks     
            call chunk_close(c%chunks(ichunk))
        end do

        if (associated( c%interpolated_trace )) then
            call trace_destroy( c%interpolated_trace )
            deallocate( c%interpolated_trace )
        end if

    end subroutine
    
    subroutine chunk_open_read( c )
    
      ! open chunk for read
      ! this opens the associated database file
      ! and reads the index dataset.
      ! it also initializes the traces cache
        
        type(t_chunk), intent(inout) :: c
        
        integer(hid_t) :: error, egal
        integer(hsize_t), dimension(3) :: dims
        integer :: ig,iz,ixc
        
        if ( c%readmode ) return
        if ( c%writemode ) call chunk_close( c )

      ! the file will stay open
        call h5eset_auto_f(0,egal)
        call h5fopen_f(char(c%filename), H5F_ACC_RDONLY_F,  c%file, error)
        call h5eset_auto_f(1,egal)
        if (error /= 0) call die( "gfdb: failed to open file: "//c%filename )
        
      ! dataset must be kept open, so that dereferencing works
        call h5dopen_f(c%file, "index", c%dataset_index, error)
        if (error /= 0) call die( "gfdb: failed to open index dataset in file: " &
                                    // c%filename )

        c%readmode = .true.
        
        allocate( c%references(c%ng, c%nz/c%nipz, c%nxc/c%nipx) )
        allocate( c%traces(c%ng, c%nz, c%nxc) )
        
        forall (ig=1:c%ng, iz=1:c%nz, ixc=1:c%nxc) 
            c%traces(ig,iz,ixc)%p => null()
        end forall
        
        dims(1) =  c%ng
        dims(2) =  c%nz/c%nipz
        dims(3) =  c%nxc/c%nipx
        call workaround(c%dataset_index, c%references, dims, error )
        if (error /= 0) call die( "gfdb: failed to read index dataset from file: " &
                                    // c%filename )
        c%readmode = .true.

    end subroutine 
    
    subroutine chunk_open_write( c )
    
      ! open chunk for write
      ! this opens the associated database file 
      ! and the index dataset.
        
        type(t_chunk), intent(inout) :: c
        integer(hsize_t), dimension(3) :: dims
        integer(hid_t) :: error, egal
        
        if ( c%writemode ) return
        if ( c%readmode ) call chunk_close( c )
        
     
      ! file will is kept open
        call h5eset_auto_f(0,egal)
        call h5fopen_f(char(c%filename), H5F_ACC_RDWR_F,  c%file, error)
        call h5eset_auto_f(1,egal)
        if (error /= 0) call die( "gfdb: failed to open file: "//c%filename )
        
      ! the index dataset is also be kept open
        
        call h5dopen_f(c%file, "index", c%dataset_index, error)
        if (error /= 0) call die( "gfdb: failed to open index dataset in file: " &
                                    // c%filename )
                                    
      ! read references to detect duplicate inserts
        
        dims(1) =  c%ng
        dims(2) =  c%nz/c%nipz
        dims(3) =  c%nxc/c%nipx
        allocate( c%references(c%ng, c%nz/c%nipz, c%nxc/c%nipx) )
        call workaround(c%dataset_index, c%references, dims, error )
        if (error /= 0) call die( "gfdb: failed to read index dataset from file: " &
                                    // c%filename )
        
        c%writemode = .true.

    end subroutine
    
    subroutine chunk_close( c )
        
        type(t_chunk), intent(inout) :: c   
        integer(hid_t) :: err
        integer(hid_t), dimension(2) :: e
        integer :: ixc,iz,ig
        
        if (c%file /= 0) then
            call h5dclose_f(c%dataset_index,e(1))
            call h5fclose_f( c%file, e(2) )
            if (any (e /= 0)) then 
                call warn("gfdb: problems closing file for chunk "//c%ichunk )
            end if
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
                            call trace_destroy( c%traces(ig,iz,ixc)%p )
                            deallocate( c%traces(ig,iz,ixc)%p )
                        end if
                    end do
                end do
            end do
            deallocate( c%traces )
        end if
        
        call h5garbage_collect_f(err)
        
    end subroutine
    
    subroutine h5_opencreategroup( file, name, group, error, sizehint )


        integer(hid_t), intent(in) :: file
        character(len=*), intent(in) :: name
        integer(hid_t), intent(out) :: group, error
        integer(size_t), intent(in) :: sizehint
        
        integer(hid_t) :: egal

        call h5eset_auto_f(0,egal)
        call h5gopen_f( file, name, group, error )
        if (error /= 0) then
            call h5gcreate_f( file, name, group, error, sizehint )
        end if
        
        call h5eset_auto_f(1,egal)
    
    end subroutine
    
    subroutine workaround(dataset, refs, dims, error )
      
      ! reading directly into a 3D-array does not seem to 
      ! work for references, as it works for other datatypes...
      
        integer(hid_t), intent(in) :: dataset
        integer(hsize_t), dimension(3), intent(in) :: dims
        type(hobj_ref_t_f), dimension(:,:,:), intent(out) :: refs
        integer(hid_t), intent(out) :: error
        
        type(hobj_ref_t_f), dimension(:), allocatable :: refs2
        integer :: i,ig,iz,ix  
        integer(hsize_t), dimension(1) :: dims1
       
        dims1(1) = dims(1)*dims(2)*dims(3)
        
        allocate(refs2(int(dims1(1))))
        call h5dread_f(dataset, H5T_STD_REF_OBJ, refs2, dims1, error )
        i= 1
        do ix=1,dims(3)
            do iz=1,dims(2)
                do ig=1,dims(1)
                    refs(ig,iz,ix) = refs2(i)                    
                    i=i+1
                end do
            end do
        end do
        deallocate(refs2) 
        
    end subroutine


end module



! the rest below here is boring...

subroutine h5_save_scalar_integer( file, datasetname, value, error )
        use hdf5

    integer(hid_t), intent(in) :: file
    character(len=*), intent(in) :: datasetname
    integer,intent(in) :: value
    integer(hid_t), intent(out) :: error
    
    integer(hid_t) :: dataspace, dataset
    integer(hsize_t),dimension(1) :: dims
    
    dims = (/1/)

    call h5screate_f( H5S_SCALAR_F,  dataspace,error)
    if (error /= 0) return
    
    call h5dcreate_f(file,datasetname,H5T_NATIVE_INTEGER,dataspace, &
                        dataset,error)
    if (error /= 0) return
    
    call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,value,dims,error)
    if (error /= 0) return
    
    call h5dclose_f(dataset,  error)
    if (error /= 0) return
    
    call h5sclose_f(dataspace,  error)
    
end subroutine

subroutine h5_save_scalar_real( file, datasetname, value, error )
        use hdf5

    integer(hid_t), intent(in) :: file
    character(len=*), intent(in) :: datasetname
    real,intent(in) :: value
    integer(hid_t), intent(out) :: error
    
    integer(hid_t) :: dataspace, dataset
    integer(hsize_t),dimension(1) :: dims
    
    dims = (/1/)

    call h5screate_f( H5S_SCALAR_F,  dataspace,error)
    if (error /= 0) return
    
    call h5dcreate_f(file,datasetname,H5T_NATIVE_REAL,dataspace, &
                        dataset,error)
    if (error /= 0) return
    
    call h5dwrite_f(dataset,H5T_NATIVE_REAL,value,dims,error)
    if (error /= 0) return
    
    call h5dclose_f(dataset,  error)
    if (error /= 0) return
    
    call h5sclose_f(dataspace,  error)
    
end subroutine

subroutine h5_open_scalar_integer( file, datasetname, value, error )
    use hdf5

    integer(hid_t), intent(in) :: file
    character(len=*), intent(in) :: datasetname
    integer,intent(out) :: value
    integer(hid_t), intent(out) :: error
    
    integer(hid_t) :: dataset
    integer(hsize_t),dimension(1) :: dims
    
    dims = (/1/)
    
    call h5dopen_f(file, datasetname, dataset, error)
    if (error /= 0) return
    
    call h5dread_f(dataset,H5T_NATIVE_INTEGER,value,dims,error)
    if (error /= 0) return
    
    call h5dclose_f(dataset,  error)
            
end subroutine

subroutine h5_open_scalar_real( file, datasetname, value, error )
    use hdf5

    integer(hid_t), intent(in) :: file
    character(len=*), intent(in) :: datasetname
    real,intent(out) :: value
    integer(hid_t), intent(out) :: error
    
    integer(hid_t) :: dataset
    integer(hsize_t),dimension(1) :: dims
    
    dims = (/1/)
    
    call h5dopen_f(file, datasetname, dataset, error)
    if (error /= 0) return
    
    call h5dread_f(dataset,H5T_NATIVE_REAL,value,dims,error)
    if (error /= 0) return
    
    call h5dclose_f(dataset,  error)
            
end subroutine
