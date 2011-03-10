module gfdb_io_hdf

    use hdf5
    
    implicit none

    logical :: h5inited = .false.
    
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
    
  contains
  
    subroutine gfdb_io_init( ok )
   
        logical, intent(out) :: ok
        integer(hid_t) :: e
        
        if (.not. h5inited) then
            call h5open_f( e )
            if (e /= 0) then
                ok = .false.
                call error( "cannot initialize hdf5 library" )
                return
            end if
            h5inited = .true.
        end if

    end subroutine
  
  
    subroutine gfdb_io_read_index(filename, dt,dx,dz, firstx, firstz, &
                                                nchunks, nx, nxc, nz, ng, ok )
            
            type(varying_string), intent(in) :: filename
            real, intent(out) :: dt, dx, dz, firstx, firstz
            integer, intent(out) :: nchunks, nx, nxc, nz, ng
            logical, intent(out) :: ok
            
            integer(hid_t) :: file, e(10), egal
            
            ok = .true.
            
            firstx = 0
            firstz = 0
            
            call h5eset_auto_f(0,egal)
            call h5fopen_f(char(filename), H5F_ACC_RDONLY_F, file, e(1))
            call h5eset_auto_f(1,egal)
            if (e(1) /= 0) then
                call error( "gfdb: failed to open file: " //filename )
                ok = .false.
                return
            end if
          ! read dt,dx,dz
            call h5_open_scalar( file, "dt", dt, e(1) )
            call h5_open_scalar( file, "dx", dx, e(2) )
            call h5_open_scalar( file, "dz", dz, e(3) )
            if (any(e/=0)) then
                call error( "gfdb: failed to read dataset from file: "//filename )
                ok = .false.
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
            
          ! save nx, nxc, nz, ng
            call h5_open_scalar( file, "nchunks", nchunks, e(1) )
            call h5_open_scalar( file, "nx", nx, e(2) )
            call h5_open_scalar( file, "nxc", nxc, e(3) )
            call h5_open_scalar( file, "nz", nz, e(4) )
            call h5_open_scalar( file, "ng", ng, e(5) )
            if (any(e/=0)) then
                call error( "gfdb: failed to read dataset from file: "//filename )
                ok = .false.
                return
            end if

            call h5fclose_f( file, e(1) )
            if (any(e/=0)) then
                call error( "gfdb: problems closing file: "//filename )
                ok = .false.
                return
            end if
            
    end subroutine
  
    subroutine gfdb_io_create_index(filename, dt,dx,dz, firstx, firstz, &
                                                nchunks, nx, nxc, nz, ng, ok
                                                
        type(varying_string), intent(in) :: filename
        real, intent(in) :: dt, dx, dz, firstx, firstz
        integer, intent(in) :: nchunks, nx, nxc, nz, ng
        logical, intent(out) :: ok
            
        integer(hid_t) :: file, e(7), egal
            
        ok = .true.
        e = 0
        

        call h5eset_auto_f(0,egal)
        call h5fcreate_f(char(filename),H5F_ACC_TRUNC_F,  file,e(1))
        call h5eset_auto_f(1,egal)
        if (e(1) /= 0) then
            call error( "gfdb: failed to create file: "//filename )
            ok = .false.
            return
        end if
        
      ! save dt,dx,dz
        call h5_save_scalar( file, "dt", c%dt, e(1) )
        call h5_save_scalar( file, "dx", c%dx*c%nipx, e(2) )
        call h5_save_scalar( file, "dz", c%dz*c%nipz, e(3) )
        call h5_save_scalar( file, "firstx", c%firstx, e(4) )
        call h5_save_scalar( file, "firstz", c%firstz, e(5) )
        if (any(e/=0)) then
            call error( "gfdb: failed to write dataset to file: "//filename )
            ok = .false.
            return
        end if
        
      ! save nx, nxc, nz, ng
        call h5_save_scalar( file, "nchunks", c%nchunks, e(1) )
        call h5_save_scalar( file, "nx", c%nx/c%nipx, e(2) )
        call h5_save_scalar( file, "nxc", c%nxc/c%nipx, e(3) )
        call h5_save_scalar( file, "nz", c%nz/c%nipz, e(4) )
        call h5_save_scalar( file, "ng", c%ng, e(5) )
        if (any(e/=0)) then
            call error( "gfdb: failed to write dataset to file: "//filename )
            ok = .false.
            return
        end if
        
        call h5fclose_f(file,  e(1))
        if (any(e/=0)) then
            call error( "gfdb: problems closing file: "//filename )
            ok = .false.
            return
        end if
  
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
