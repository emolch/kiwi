! $Id: gfdb_build.f90 623 2007-04-25 14:49:07Z sebastian $ 
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

module gfdb_build_ahfull_

    use util
    use unit
    use gfdb
    use constants
    use better_varying_string
    use sparse_trace
    use seismogram_io   
    use elseis_oo
    
    implicit none
      
    private  
    public addentry, set_database_material_stf, close_database
    
    type(t_gfdb), save   :: db
    real, dimension(3,3) :: source_a = reshape((/1,1,0,1,0,0,0,0,0/),(/3,3/))
    real, dimension(3,3) :: source_b = reshape((/0,0,1,0,0,1,1,1,0/),(/3,3/))
    real, dimension(3,3) :: source_c = reshape((/0,0,0,0,0,0,0,0,1/),(/3,3/))
    real                 :: rho, alpha, beta
    type(elseis_t), save :: es
    integer :: oversample
    integer :: traces_added = 0

    
  contains
  
    subroutine set_database_material_stf( basefn, material_tab, stf_tab, oversample_ )
        
        type(varying_string), intent(in) :: basefn
        real, dimension(:,:), intent(in) :: material_tab, stf_tab
        integer, intent(in) :: oversample_
                
        call gfdb_init(db, basefn)
        
        ! prepare material
        rho = material_tab(1,1)
        alpha = material_tab(2,1)
        beta = material_tab(3,1)        
        call set_material( es, rho, alpha, beta )
      
        oversample = oversample_
        
      ! prepare stf
        call set_stf( es, stf_tab(2,:), db%dt/oversample )
    
    end subroutine
     
    subroutine close_database()
        
        call es_destroy( es )
        call gfdb_destroy( db )
        
    end subroutine
    
    subroutine addentry( buffer, ok )
    
        character(len=*), intent(in) :: buffer
        logical, intent(out)         :: ok
                
        integer :: iostat, ix, iz
        real :: x, z
       
        real, dimension(3) :: r_location, s_location, rel_location
        real*8 :: epiangle
        real :: tbegin_total, tend_total, d, tstf
        real, dimension(2)  :: tbegin, tend
        logical :: nfflag, ffflag
        real, dimension(:,:), allocatable :: seismograms, seismograms2
        real, dimension(3,3) :: rotmat
        real :: firstarrival_p, lastarrival_p, firstarrival_s, lastarrival_s 
        integer :: n,p,q, nsamples, itbegin, itend, nsamples2, jm,jp, i
        integer :: nwindows, iwindow
                  
        read (unit=buffer,fmt=*,iostat=iostat) x, z,  nfflag, ffflag
        
        if (iostat == 0) then
            ok = .true.
            
          ! det locations of source, receiver and relative location
            s_location(1) = 0
            s_location(2) = 0
            r_location(2) = 0
            s_location(3) = z
            epiangle = x/earthradius
            r_location(1) = real( sin(epiangle)*earthradius )
            r_location(3) = real( (1.-cos(epiangle))*earthradius )
            rel_location = r_location - s_location
            d = dist(s_location, r_location)
            
          ! time window  
            
            tstf = (size(es%stf)-1)*es%dt
            
            firstarrival_p = snapdown(d/alpha, db%dt)
            lastarrival_p  = snapup(d/alpha + tstf,db%dt)
            firstarrival_s = snapdown(d/beta ,db%dt)
            lastarrival_s  = snapup(d/beta + tstf,db%dt) + db%dt*2 ! add 2 samples of zero/static at the end
            
            tbegin_total = firstarrival_p
            tend_total = lastarrival_s
            
            if (lastarrival_p >= firstarrival_s .or. nfflag) then   ! s and p time windows overlap or in near field
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
            

            nsamples = nint((tend_total - tbegin_total) / es%dt + 1) 
            
            allocate( seismograms(9,nsamples) )
            seismograms(:,:) = 0.
            
            call set_coords( es, rel_location,  nfflag, ffflag )
            
            do n=1,3
                do p=1,3
                    do q=1,3
                        call set_npq( es, n, p, q )
                        do iwindow=1,nwindows
                            itbegin = nint((tbegin(iwindow)-tbegin_total) / es%dt)+1
                            itend = nint((tend(iwindow)-tbegin_total) / es%dt)+1
                            call add_elseis_mt( es, &
                                                tbegin(iwindow), &
                                                source_a(p,q), &
                                                seismograms(n,itbegin:itend) )
                            call add_elseis_mt( es, &
                                                tbegin(iwindow), &
                                                source_b(p,q), &
                                                seismograms(n+3,itbegin:itend) )
                            call add_elseis_mt( es, &
                                                tbegin(iwindow), &
                                                source_c(p,q), &
                                                seismograms(n+6,itbegin:itend) )
                        end do
                    end do
                end do
            end do
            
            
          ! downsample 
            nsamples2 = (nsamples-1)/oversample + 2
            allocate( seismograms2(9,nsamples2) )
            seismograms2(:,:) = 0.
            
            do i=1,nsamples
                jm = (i-oversample/2-1)/oversample+1
                jp = (i-oversample/2)/oversample+1
                if (jm==jp) then
                    seismograms2(:,jm) = seismograms2(:,jm) + seismograms(:,i)/oversample
                else
                    seismograms2(:,jm) = seismograms2(:,jm) + seismograms(:,i)/oversample*0.5
                    seismograms2(:,jp) = seismograms2(:,jp) + seismograms(:,i)/oversample*0.5
                end if
            end do
            ! to get the last sample(s) right... (repeating end point)
            do i=nsamples+1,nsamples+1+oversample*3
                jm = (i-oversample/2-1)/oversample+1
                jp = (i-oversample/2)/oversample+1
                if (jm==jp) then
                    if (jm<=nsamples2) seismograms2(:,jm) = seismograms2(:,jm) + seismograms(:,nsamples)/oversample
                else
                    if (jm<=nsamples2) seismograms2(:,jm) = seismograms2(:,jm) + seismograms(:,nsamples)/oversample*0.5
                    if (jp<=nsamples2) seismograms2(:,jp) = seismograms2(:,jp) + seismograms(:,nsamples)/oversample*0.5
                end if
            end do
            
            
            
            deallocate(seismograms)
            
          ! rotate to local receiver coordinates. here: (South,Right,Down) 
            rotmat(:,:) = 0.
            rotmat(1,1) = real(cos(epiangle))
            rotmat(2,2) = 1.
            rotmat(3,3) = real(cos(epiangle))
            rotmat(1,3) = real(sin(epiangle))
            rotmat(3,1) = real(-sin(epiangle))
            
            call rotate_seismogram( rotmat, seismograms2(1:3,:) )
            call rotate_seismogram( rotmat, seismograms2(4:6,:) )
            call rotate_seismogram( rotmat, seismograms2(7:9,:) )
            
          ! put in database
            call gfdb_get_indices( db, x, z, ix, iz )
            
            call gfdb_save_array( db, ix,iz ,1, tbegin_total, seismograms2(1,:) ) 
            call gfdb_save_array( db, ix,iz ,2, tbegin_total, seismograms2(4,:) ) 
            call gfdb_save_array( db, ix,iz ,3, tbegin_total, seismograms2(7,:) ) 
            call gfdb_save_array( db, ix,iz ,4, tbegin_total, seismograms2(2,:) ) 
            call gfdb_save_array( db, ix,iz ,5, tbegin_total, seismograms2(5,:) ) 
            call gfdb_save_array( db, ix,iz ,6, tbegin_total, seismograms2(3,:) ) 
            call gfdb_save_array( db, ix,iz ,7, tbegin_total, seismograms2(6,:) ) 
            call gfdb_save_array( db, ix,iz ,8, tbegin_total, seismograms2(9,:) ) 
            traces_added = traces_added + 8
            
            
          ! periodically close the gfdb, so that hdf is forced to deallocate
          ! all it's memory
            if (traces_added > 30000) then
                call gfdb_close( db )
                traces_added = 0
            end if
            
            deallocate(seismograms2)
            
        end if
              
    end subroutine 
    
    subroutine gfdb_save_array( db, ix,iz,ig, tbegin, array ) 
        
        type(t_gfdb), intent(inout) ::db
        integer, intent(in) :: ix,iz,ig
        real, intent(in) :: tbegin
        real, dimension(:), intent(in) :: array
      
      ! wrapper to gfdb_save_trace to store an array instead of a trace
        
        type(t_strip) :: conti
        type(t_trace) :: tr
        integer, dimension(2) :: span
        
        span(1) = nint(tbegin/db%dt)
        span(2) = span(1) + size(array) - 1
        
        call strip_init( span, array, conti )
        call trace_pack( conti, tr )
        call gfdb_save_trace( db, ix, iz, ig, tr )
    
        call strip_destroy( conti )
        call trace_destroy( tr )
        
    end subroutine
    
    pure subroutine rotate_seismogram( rotmat, xyz )
        
        real, dimension(3,3), intent(in)     :: rotmat
        real, dimension(:,:), intent(inout)  :: xyz
        integer                              :: i
        
    ! inplace apply rotation to vectors in xyz
   
        do i=1,size(xyz,2)
            xyz(:,i) = matmul( rotmat, xyz(:,i) )
        end do
        
    end subroutine rotate_seismogram
    
    pure function dist(a,b) result(r)
    
        real, intent(in), dimension(3) :: a,b
        real :: r
        
        r = sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)

    end function dist
    
    pure real function snapdown(t,dt) 
        real, intent(in) :: t, dt
        snapdown = floor(t/dt)*dt
    end function
    
    pure real function snapup(t,dt)
        real, intent(in) :: t, dt
        snapup = ceiling(t/dt)*dt
    end function
    
end module

program gfdb_build_ahfull
  
  ! This program is used to create and fill a Greens function database to be used with the kiwi programs
  ! It calculates Greens functions for an analytical fullspace
  !
  ! usage: gfdb_build_ahfull database material stf <<EOF
  ! x z 
  ! ...
  ! EOF
  !
  ! Complete documentation is available on
  ! 
  !   http://kinherd.org/power/trac/wiki/GFDBBuildAhfullTool
  !
  
    use util
    use unit
    use better_varying_string
    use varying_string_getarg
    use read_line
    use read_table
    use gfdb_build_ahfull_
    
    ! use f90_unix_env
    
    implicit none
    
    type(varying_string) :: basefn, string
    integer              :: iostat, iline
    logical              :: ok
    character, parameter :: eol = char(10)
    real, allocatable, dimension(:,:) :: material_tab, stf_tab
    integer :: oversample
    
    g_pn = 'gfdb_build_ahfull'
    g_usage = 'usage: ' // g_pn // ' database material stf <<EOF' // eol // & 
              "x z nfflag ffflag" // eol // &
              "..." // eol // &
              "EOF" // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/GFDBBuildAhfullTool"
    
    if (iargc() /= 4) call usage()
    
    call vs_getarg( 1, basefn )
    call readtable_arg( 2, 3, 1, material_tab )
    call readtable_arg( 3, 2, 2, stf_tab )
    call vs_getarg( 4, string )
    oversample = string
    
    call set_database_material_stf( basefn, material_tab, stf_tab, oversample )
    
    iline = 1
    line_loop : do
        call readline( addentry, iostat, ok )
        if (iostat == IOSTAT_EOF) exit line_loop
        if (.not. ok) call die( "reading line "// iline &
                                 //" from standard input failed." )        
        iline = iline+1
    end do line_loop    
    
    call close_database()

    call cleanup()
    
  contains
    
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
        if (size(field,2) <  min_rows) call die( 'expected at least ' // min_rows // ' row(s) in file "' // fn // '"' )

    end subroutine
    
end program
