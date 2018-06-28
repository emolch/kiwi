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

module gfdb_io_hdf

    use hdf5
    use util 
    use better_varying_string

    implicit none

    logical :: h5inited = .false.
    
    public gfdb_io_init
    public gfdb_io_deinit
    public gfdb_io_read_index
    public gfdb_io_create_index
    public gfdb_io_create_chunk
    public gfdb_io_save_trace
    public gfdb_io_get_trace
    public gfdb_io_chunk_open_read
    public gfdb_io_chunk_open_write
    public gfdb_io_chunk_close
    public gfdb_io_chunk_read_index


    interface h5_save_scalar

        
        subroutine h5_save_scalar_integer( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            integer,intent(in) :: value
            integer, intent(out) :: error

        end subroutine
        
        subroutine h5_save_scalar_real( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            real,intent(in) :: value
            integer, intent(out) :: error

        end subroutine
        
    end interface
    
    interface h5_open_scalar
        
        subroutine h5_open_scalar_integer( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            integer,intent(out) :: value
            integer, intent(out) :: error

        end subroutine
        
        subroutine h5_open_scalar_real( file, datasetname, value, error )
            use hdf5
            integer(hid_t), intent(in) :: file
            character(len=*), intent(in) :: datasetname
            real,intent(out) :: value
            integer, intent(out) :: error

        end subroutine
        
    end interface
    

  contains
  
    subroutine gfdb_io_init( ok )
   
        logical, intent(out) :: ok
        integer :: e
        
        ok = .false.
     
        if (.not. h5inited) then
            call h5open_f( e )
            if (e /= 0) then
                call error( "cannot initialize hdf5 library" )
                return
            end if
            h5inited = .true.
        end if

        ok = .true.

    end subroutine
  
    subroutine gfdb_io_deinit( )

        integer :: egal

        if (h5inited) then
            call h5close_f(egal)
            h5inited = .false.
        end if

    end subroutine

    subroutine gfdb_io_read_index(filename, dt,dx,dz, firstx, firstz, &
                                                nchunks, nx, nxc, nz, ng, ok )
            
        type(varying_string), intent(in) :: filename
        real, intent(out) :: dt, dx, dz, firstx, firstz
        integer, intent(out) :: nchunks, nx, nxc, nz, ng
        logical, intent(out) :: ok
        
        integer(hid_t) :: file
        integer :: e(10), egal
        
        ok = .false.
        
        e = 0

        firstx = 0
        firstz = 0
        
        call h5eset_auto_f(0,egal)
        call h5fopen_f(char(filename), H5F_ACC_RDONLY_F, file, e(1))
        call h5eset_auto_f(1,egal)
        if (e(1) /= 0) then
            call error( "gfdb: failed to open file: " //filename )
            return
        end if
      
        call h5_open_scalar( file, "dt", dt, e(1) )
        call h5_open_scalar( file, "dx", dx, e(2) )
        call h5_open_scalar( file, "dz", dz, e(3) )
        if (any(e/=0)) then
            call error( "gfdb: failed to read dataset from file: "//filename )
            return
        end if
    
      ! try to read firstx, firstz; leave these at zero when not found 
      ! (for backward compatibility)            
        call h5eset_auto_f(0,egal)
        call h5_open_scalar( file, "firstx", firstx, e(4) )
        call h5_open_scalar( file, "firstz", firstz, e(5) )
        call h5eset_auto_f(1,egal)
        call h5eclear_f(egal)            
        e=0
        
        call h5_open_scalar( file, "nchunks", nchunks, e(1) )
        call h5_open_scalar( file, "nx", nx, e(2) )
        call h5_open_scalar( file, "nxc", nxc, e(3) )
        call h5_open_scalar( file, "nz", nz, e(4) )
        call h5_open_scalar( file, "ng", ng, e(5) )
        if (any(e/=0)) then
            call error( "gfdb: failed to read dataset from file: "//filename )
            return
        end if

        call h5fclose_f( file, e(1) )
        if (any(e/=0)) then
            call error( "gfdb: problems closing file: "//filename )
            return
        end if

        ok = .true.
            
    end subroutine
  
    subroutine gfdb_io_create_index(filename, dt,dx,dz, firstx, firstz, &
                                                nchunks, nx, nxc, nz, ng, ok)
                                                
        type(varying_string), intent(in) :: filename
        real, intent(in) :: dt, dx, dz, firstx, firstz
        integer, intent(in) :: nchunks, nx, nxc, nz, ng
        logical, intent(out) :: ok
            
        integer(hid_t) :: file
        integer :: e(7), egal
            
        ok = .false.
        e = 0
        
        call h5eset_auto_f(0,egal)
        call h5fcreate_f(char(filename),H5F_ACC_TRUNC_F,  file,e(1))
        call h5eset_auto_f(1,egal)
        if (e(1) /= 0) then
            call error( "gfdb: failed to create file: "//filename )
            return
        end if
        
      ! save dt,dx,dz
        call h5_save_scalar( file, "dt", dt, e(1) )
        call h5_save_scalar( file, "dx",dx, e(2) )
        call h5_save_scalar( file, "dz", dz, e(3) )
        call h5_save_scalar( file, "firstx", firstx, e(4) )
        call h5_save_scalar( file, "firstz", firstz, e(5) )
        if (any(e/=0)) then
            call error( "gfdb: failed to write dataset to file: "//filename )
            return
        end if
        
      ! save nx, nxc, nz, ng
        call h5_save_scalar( file, "nchunks", nchunks, e(1) )
        call h5_save_scalar( file, "nx", nx, e(2) )
        call h5_save_scalar( file, "nxc", nxc, e(3) )
        call h5_save_scalar( file, "nz", nz, e(4) )
        call h5_save_scalar( file, "ng", ng, e(5) )
        if (any(e/=0)) then
            call error( "gfdb: failed to write dataset to file: "//filename )
            return
        end if
        
        call h5fclose_f(file,  e(1))
        if (any(e/=0)) then
            call error( "gfdb: problems closing file: "//filename )
            return
        end if
        
        ok = .true.
  
    end subroutine
  
    subroutine gfdb_io_create_chunk(filename, nxc, nz, ng, size_hint, ok)
        
        type(varying_string), intent(in) :: filename
        integer, intent(in)              :: ng, nz, nxc
        integer(size_t), intent(in)      :: size_hint
        logical, intent(out)             :: ok        

        type(hobj_ref_t_f), dimension(:), allocatable :: refs
        integer(hid_t) :: file, dataspace, dataset, group
        integer :: e, egal
        integer(hsize_t), dimension(3) :: dims
        integer(hsize_t), dimension(1) :: dims1
        integer :: i, length

        ok = .false.
        
        dims(1) = ng
        dims(2) = nz
        dims(3) = nxc

        call h5eset_auto_f(0,egal)
        call h5fcreate_f(char(filename),H5F_ACC_TRUNC_F,  file, e)
        call h5eset_auto_f(1,egal)
        if (e /= 0) then
            call error( "gfdb: failed to create file: "//filename )
            return
        end if

      ! create index dataset for faster access to individual traces
        call h5screate_simple_f(3,dims,  dataspace, e)
        if (e /= 0) then
            call error( "gfdb: failed to write dataset to file: "//filename )
            return
        end if
      
        call h5dcreate_f(file,"index",H5T_STD_REF_OBJ,dataspace,  dataset, e)
        if (e /= 0) then
            call error( "gfdb: failed to create index dataset in file: "//filename )
            return
        end if

        dims1(1) = dims(1)*dims(2)*dims(3)
        length = int(dims1(1)) ! remove a compiler warning 
        allocate(refs(length))
        
        do i=1,dims1(1)
            refs(i)%ref = 0
        end do
        
        call h5dwrite_f(dataset, H5T_STD_REF_OBJ, refs, dims1, e )
        if (e /= 0) then
            call error( "gfdb: failed to write index dataset to file: "//filename )
            return
        end if

        deallocate(refs)
            
        call h5gcreate_f(file,"gf",group, e, size_hint)
        if (e /= 0) then
            call error( "gfdb: failed to create group gf in file: "//filename )
            return
        end if
        
        call h5gclose_f(group, e)    
        
        call h5dclose_f(dataset,  e)
        call h5sclose_f(dataspace,  e)        
        call h5fclose_f(file,  e)
        if (e /= 0) then
            call error( "gfdb: problems closing file: "//filename )
            return
        end if
        
        ok = .true.

    end subroutine


    subroutine gfdb_io_save_trace(filename, file, dataset_index, &
            ixc, iz, ig, &
            packed, packed_size, pofs, ofs, nstrips, &
            size_hint1, size_hint2, &
            reference, ok )

        type(varying_string), intent(in)                :: filename
        integer(hid_t), intent(in)                      :: file,  dataset_index
        integer, intent(in)                             :: ixc, iz, ig
        integer, dimension(:), intent(in)               :: ofs, pofs
        real, dimension(:), intent(in)                  :: packed
        integer, intent(in)                             :: packed_size, nstrips
        integer(size_t), intent(in)                     :: size_hint1, size_hint2
        type(hobj_ref_t_f), dimension(1), intent(out)   :: reference
        logical, intent(out)                            :: ok

        integer, dimension(8) :: e
        integer(hid_t) :: dataset, dataspace, dataspace1, attribute, group_dist
        integer(hid_t) :: dataspace_for_ref, memspace, group
        integer(hsize_t),dimension(1) :: dims
        
        type(varying_string) :: datasetname, gloc
        integer(hsize_t), dimension(3,1) :: coord

        e=0
        reference(1)%ref = 0
        ok = .false.
        
      ! open/create group for the dataset
                
        gloc = var_str("/gf/") // ixc
        call h5_opencreategroup( file, char(gloc), group_dist, e(1), size_hint1 )
        gloc = iz
        call h5_opencreategroup( group_dist, char(gloc), group, e(2), size_hint2 )
        if (any(e /= 0)) then
            call error( "gfdb: failed to open/create group "//gloc//" in file: " &
                                      // filename )
            return
        end if
      
      ! create dataset for the array
        
        dims(1) = packed_size
        call h5screate_simple_f(1,dims,  dataspace,e(1))
        datasetname = ig
        call h5dcreate_f(group,char(datasetname), &
                         H5T_NATIVE_REAL,dataspace,  dataset,e(2))
        if (any(e /= 0)) then
            call error( "gfdb: failed to create dataset "// &
                                    datasetname // " in file: " // filename )
            return
        end if
       
      ! save offsets as attributes to the dataset  
                         
        dims(1) = nstrips
        call h5screate_simple_f(1,dims, dataspace1, e(1))
        
        call h5acreate_f(dataset,"pofs",H5T_NATIVE_INTEGER,dataspace1, &
                         attribute,e(2))
        call h5awrite_f(attribute,H5T_NATIVE_INTEGER,pofs,dims,e(3))
        call h5aclose_f(attribute, e(4))
        
        call h5acreate_f(dataset,"ofs",H5T_NATIVE_INTEGER,dataspace1, &
                         attribute,e(5))
        call h5awrite_f(attribute,H5T_NATIVE_INTEGER,ofs,dims,e(6))
        call h5aclose_f(attribute, e(7))
        call h5sclose_f(dataspace1, e(8))
        if (any(e /= 0)) then
            call error( "gfdb: failed to save attributes to dataset "// &
                                    datasetname // " in file: " // filename )
            return
        end if
       
      ! save the array 
      
        dims(1) = packed_size
        call h5dwrite_f(dataset, H5T_NATIVE_REAL, packed, dims, e(1))
        if (any(e /= 0)) then
            call error( "gfdb: failed to write dataset "// &
                                    datasetname // " in file: " // filename )
            return
        end if
      
      ! save a reference to the dataset in the index block  
        
        call h5rcreate_f(group, char(datasetname), reference(1), e(1))
        coord(:,1) = (/ ig, iz, ixc /)
        dims(1) = 1
        call h5dget_space_f(dataset_index,dataspace_for_ref,e(2))
        call h5sselect_elements_f(dataspace_for_ref,H5S_SELECT_SET_F,3, int(1,SIZE_T), coord, e(3))
        call h5screate_simple_f(1, dims, memspace, e(4))
        
        call h5dwrite_f(dataset_index, H5T_STD_REF_OBJ, reference, dims, e(5), &
                         mem_space_id=memspace, file_space_id=dataspace_for_ref)
                         
        call h5sclose_f(memspace, e(6))
        call h5sclose_f(dataspace_for_ref, e(7))
        if (any(e /= 0)) then
            call error( "gfdb: failed save reference to "// &
                                    datasetname // " in index dataset of file: " &
                                    // filename )
            return
        end if
        
        call h5dclose_f(dataset, e(1))
        call h5sclose_f(dataspace, e(1))
        call h5gclose_f(group_dist, e(1))
        call h5gclose_f(group, e(1))

        ok = .true.

    end subroutine


    subroutine gfdb_io_get_trace( dataset_index, reference, &
            packed, packed_size, pofs, ofs, nstrips, &
            ok)
            
        integer(hid_t), intent(in)                        :: dataset_index
        type(hobj_ref_t_f), intent(in)                    :: reference
        integer, dimension(:), allocatable, intent(inout) :: ofs, pofs
        real, dimension(:),  allocatable, intent(inout)   :: packed
        integer, intent(out)                              :: packed_size
        integer, intent(out)                              :: nstrips
        logical, intent(out)                              :: ok

        integer(hsize_t),dimension(1)   :: adims
        integer(hid_t)                  :: dataset, attribute, space
        integer, dimension(8)    :: e
        integer                  :: egal
        integer(hsize_t)                :: length
        
        e = 0
        ok = .false.

      ! lookup reference to the dataset for a quick jump

        call h5eset_auto_f(0,egal)
        call h5rdereference_f(dataset_index, reference, dataset, e(1) )
        call h5eset_auto_f(1,egal)
        call h5eclear_f(egal)
        if (e(1) /= 0) then
            call error( "gfdb: failed to dereference a dataset reference" )
            return
        end if
        
      ! get offsets
                
        call h5aopen_idx_f(dataset, 0, attribute, e(1))
        call h5aget_space_f( attribute, space, e(2) )
        call h5sget_simple_extent_npoints_f(space, length, e(3) )
        call h5sclose_f( space, e(4) ) 
        nstrips = int(length)
        
        if (any(e /= 0)) then
            call error( "gfdb: failed to get size of attributes of a dataset" )
            return
        end if

        if (.not. allocated(pofs)) &
            call resize( pofs, 1, nstrips )
       
        if (nstrips > size(pofs)) &
           call resize( pofs, 1, nstrips )

        if (.not. allocated(ofs)) &
            call resize( ofs, 1, nstrips )
            
        if (nstrips > size(ofs)) &
            call resize( ofs, 1, nstrips )

        adims(1) = nstrips
        call h5aread_f(attribute, H5T_NATIVE_INTEGER, pofs, adims, e(1))
        call h5aclose_f( attribute, e(2) )
        
        call h5aopen_idx_f(dataset, 1, attribute, e(3))
        call h5aread_f(attribute, H5T_NATIVE_INTEGER, ofs, adims, e(4))
        call h5aclose_f( attribute, e(5) )
        
        if (any(e /= 0)) then
            call error( "gfdb: failed to get attributes of a dataset" )
            return
        end if

      ! get size of data
      
        call h5dget_space_f(dataset,space,e(1))
        call h5sget_simple_extent_npoints_f(space, length, e(2))
        call h5sclose_f( space, e(3) )
        
      ! get the data
  
        packed_size = int(length)
        if (.not. allocated(packed)) &
            call resize(packed, 1, packed_size)
            
        if (packed_size > size(packed)) &
            call resize(packed, 1, packed_size)

        adims(1) = packed_size
        call h5dread_f(dataset, H5T_NATIVE_REAL, packed, adims, e(4))
        call h5dclose_f(dataset, e(5))
        if (any(e /= 0)) then
            call error( "gfdb: failed to get a dataset" )
            return
        end if

        ok = .true.

    end subroutine

    subroutine gfdb_io_chunk_open_read( filename, file, ok )

        type(varying_string), intent(in)    :: filename
        integer(hid_t), intent(out)         :: file
        logical, intent(out)                :: ok
    
        integer :: e, egal

        ok = .false.

      ! the file will stay open
        call h5eset_auto_f(0,egal)
        call h5fopen_f(char(filename), H5F_ACC_RDONLY_F,  file, e)
        call h5eset_auto_f(1,egal)
        if (e /= 0) then
            call error( "gfdb: failed to open file: "//filename )
            return
        end if

        ok = .true.

    end subroutine

    subroutine gfdb_io_chunk_open_write( filename, file, ok )

        type(varying_string), intent(in)    :: filename
        integer(hid_t), intent(out)         :: file
        logical, intent(out)                :: ok
    
        integer :: e, egal

        ok = .false.

      ! file will is kept open
        call h5eset_auto_f(0,egal)
        call h5fopen_f(char(filename), H5F_ACC_RDWR_F,  file, e)
        call h5eset_auto_f(1,egal)
        if (e /= 0) then
            call error( "gfdb: failed to open file: "//filename )
            return 
        end if
        
        ok = .true.

    end subroutine


    subroutine gfdb_io_chunk_read_index( file, dataset_index, references, ok )

        integer(hid_t), intent(in)                          :: file
        integer(hid_t), intent(out)                         :: dataset_index
        type(hobj_ref_t_f), dimension(:,:,:), intent(out)   :: references
        logical, intent(out)                                :: ok

        integer :: e
        integer(hsize_t), dimension(3) :: dims


        ok = .false.

        ! dataset must be kept open, so that dereferencing works
        call h5dopen_f(file, "index", dataset_index, e)
        if (e /= 0) then
            call error( "gfdb: failed to open index dataset" )
            return
        end if

        dims(1) =  size(references,1)
        dims(2) =  size(references,2)
        dims(3) =  size(references,3)

        call workaround(dataset_index, references, dims, e )
        if (e /= 0) then
            call error( "gfdb: failed to read index dataset" )
            return
        end if
        
        ok = .true.
        
    end subroutine

    subroutine gfdb_io_chunk_close( file, dataset_index, ok )

        integer(hid_t), intent(in)                          :: file
        integer(hid_t), intent(in)                          :: dataset_index
        logical, intent(out)                                :: ok
        integer, dimension(2) :: e

        ok = .false.
        e = 0

        call h5dclose_f(dataset_index,e(1))
        call h5fclose_f( file, e(2) )
        if (any (e /= 0)) then 
            call error("gfdb: error closing chunk file" )
            return
        end if

        ok = .true.

    end subroutine

    subroutine h5_opencreategroup( file, name, group, error, sizehint )


        integer(hid_t), intent(in) :: file
        character(len=*), intent(in) :: name
        integer(hid_t), intent(out) :: group
        integer, intent(out) :: error
        integer(size_t), intent(in) :: sizehint
        
        integer :: egal

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
        integer, intent(out) :: error
        
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

subroutine h5_save_scalar_integer( file, datasetname, value, error )
    use hdf5

    integer(hid_t), intent(in) :: file
    character(len=*), intent(in) :: datasetname
    integer,intent(in) :: value
    integer, intent(out) :: error
    
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
    integer, intent(out) :: error
    
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
    integer, intent(out) :: error
    
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
    integer, intent(out) :: error
    
    integer(hid_t) :: dataset
    integer(hsize_t),dimension(1) :: dims
    
    dims = (/1/)
    
    call h5dopen_f(file, datasetname, dataset, error)
    if (error /= 0) return
    
    call h5dread_f(dataset,H5T_NATIVE_REAL,value,dims,error)
    if (error /= 0) return
    
    call h5dclose_f(dataset,  error)
            
end subroutine
