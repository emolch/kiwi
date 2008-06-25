! $Id: source_bilat.f90 687 2008-02-13 14:39:42Z sebastian $ 
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

module source_bilat

    use orthodrome
    use parameterized_source
    use discrete_source
    use piecewise_linear_function
    use euler
    use better_varying_string
    use unit
    use util
    use constants
    
    implicit none 

    private
    
    integer, public, parameter :: n_source_params_bilat = 14
    integer, public, parameter :: n_grid_dims_bilat = 3
    
    public psm_to_tdsm_bilat, psm_write_info_file_bilat, psm_set_bilat
    public psm_get_param_name_bilat, psm_get_param_unit_bilat, psm_get_param_id_bilat
    public psm_cleanup_bilat
    
    private psm_to_tdsm_size_bilat, psm_to_tdsm_table_bilat
    
    type(varying_string), private, dimension(:), allocatable :: psm_param_names_bilat, psm_param_units_bilat
    
  ! these values may be used to hint step sizes and the like
    real, dimension(n_source_params_bilat), public, parameter :: psm_params_norm_bilat = &
    (/ 1., 10000., 10000., 10000., 7e18, 360., 90., 360., 360., 10000., 10000., 10000., 3000.,1. /)
    
    real, private, parameter :: big = huge(psm_params_norm_bilat(1))
    
  ! logical, physical or computational accuracy limits for params
    real, dimension(n_source_params_bilat), public, parameter :: psm_params_min_hard_bilat = &
    (/ -big, -100000., -100000., 0., 1., -big, -big, -big, -big, 0., 0., 0., 100., 0./)
    real, dimension(n_source_params_bilat), public, parameter :: psm_params_max_hard_bilat = &
    (/ big, 100000., 100000., 1000000., 7e25, big, big, big, big, 10000000., 10000000., 10000000., 100000., 10./)
    
  ! these limits give range of realistic, non-redundant or practical params
    real, dimension(n_source_params_bilat), public, parameter :: psm_params_min_soft_bilat = &
    (/ -20., -10000., -10000., 0., 1., -180., 0., -180., -180., 0., 0., 0., 1000., 0./)
    real, dimension(n_source_params_bilat), public, parameter :: psm_params_max_soft_bilat = &
    (/ 20., 10000., 10000., 150000., 7e25, 180., 90., 180., 180., 100000., 100000., 100000., 10000., 5./)
    
  ! defaults for external use
    real, dimension(n_source_params_bilat), public, parameter :: psm_params_default_bilat = &
    (/ 0.,0., 0., 10000., 7e18, 0., 80., 0., 0., 10000., 0., 7000., 3500., 1./)
    
   ! M0 = 7e18 Nm <=> Mw = 6.5
   
  contains
  
    subroutine psm_cleanup_bilat()
    
        integer :: i
        
        if (allocated( psm_param_names_bilat )) then
            do i=1,n_source_params_bilat
                call delete( psm_param_names_bilat(i) )
            end do
            deallocate( psm_param_names_bilat )
        end if
        if (allocated( psm_param_units_bilat )) then
            do i=1,n_source_params_bilat
                call delete( psm_param_units_bilat(i) )
            end do
            deallocate( psm_param_units_bilat )
        end if
        
    end subroutine
  
    subroutine psm_init_param_names_bilat()
    
        if (.not. allocated( psm_param_names_bilat )) then
            allocate( psm_param_names_bilat( n_source_params_bilat ) )
            psm_param_names_bilat(1) = "time"
            psm_param_names_bilat(2) = "north-shift"
            psm_param_names_bilat(3) = "east-shift"
            psm_param_names_bilat(4) = "depth"
            psm_param_names_bilat(5) = "moment"
            psm_param_names_bilat(6) = "strike"
            psm_param_names_bilat(7) = "dip"
            psm_param_names_bilat(8) = "slip-rake"
            psm_param_names_bilat(9) = "rupture-rake"
            psm_param_names_bilat(10) = "length-a"
            psm_param_names_bilat(11) = "length-b"
            psm_param_names_bilat(12) = "width"
            psm_param_names_bilat(13) = "rupture-velocity"
            psm_param_names_bilat(14) = "rise-time"
        end if
        
        if (.not. allocated( psm_param_units_bilat )) then
            allocate( psm_param_units_bilat( n_source_params_bilat ) )
            psm_param_units_bilat(1) = "s"
            psm_param_units_bilat(2) = "m"
            psm_param_units_bilat(3) = "m"
            psm_param_units_bilat(4) = "m"
            psm_param_units_bilat(5) = "Nm"
            psm_param_units_bilat(6) = "degrees"
            psm_param_units_bilat(7) = "degrees"
            psm_param_units_bilat(8) = "degrees"
            psm_param_units_bilat(9) = "degrees"
            psm_param_units_bilat(10) = "m"
            psm_param_units_bilat(11) = "m"
            psm_param_units_bilat(12) = "m"
            psm_param_units_bilat(13) = "m/s"
            psm_param_units_bilat(14) = "s"
        end if
        
    end subroutine
  
    subroutine psm_get_param_name_bilat( iparam, name )
    
        integer, intent(in)                 :: iparam
        type(varying_string), intent(inout) :: name
        
        call psm_init_param_names_bilat()
        
        name = ""
        if (iparam < 1 .or. n_source_params_bilat < iparam) return
        name = psm_param_names_bilat(iparam)
        
    end subroutine
    
    subroutine psm_get_param_unit_bilat( iparam, unit )
    
        integer, intent(in)                 :: iparam
        type(varying_string), intent(inout) :: unit
        
        call psm_init_param_names_bilat()
        
        unit = ""
        if (iparam < 1 .or. n_source_params_bilat < iparam) return
        unit = psm_param_units_bilat(iparam)
        
    end subroutine
   
    subroutine psm_get_param_id_bilat( name, iparam )
    
        type(varying_string), intent(in)     :: name
        integer, intent(out)                 :: iparam

        integer :: i
        
        call psm_init_param_names_bilat()
        iparam = 0
        do i=1, n_source_params_bilat
            if (name .eq. psm_param_names_bilat(i)) then
                iparam = i
                return
            end if
        end do
    
    end subroutine
  
    subroutine psm_set_bilat( psm,  params, normalized_ )
    
        type(t_psm), intent(inout)     :: psm
        real, dimension(:), intent(in) :: params
        logical, intent(in), optional  :: normalized_
        
        logical :: must_reset_grid
        logical :: normalized
        
        normalized = .false.
        if (present(normalized_)) normalized = normalized_
        
        must_reset_grid = .false.
        if (.not. allocated(psm%grid_size) .or. (psm%sourcetype .ne. psm_bilat)) then
            must_reset_grid = .true.
        end if
        
        call resize( psm%params, 1, n_source_params_bilat )
        call resize( psm%grid_size, 1, n_grid_dims_bilat )
        call resize( psm%params_norm, 1, n_source_params_bilat )
        
        psm%params_norm(:) = psm_params_norm_bilat(:)
        
        if (must_reset_grid) then
            psm%grid_size = 1
        end if
        
        if (size(params,1) .ne. size(psm%params)) call die("wrong number of source parameters in psm_set_bilat()")
                
        if (normalized) then
            psm%params = params * psm%params_norm
        else
            psm%params = params 
        end if
        
        call psm_update_dep_params_bilat( psm )
        
    end subroutine
        
    subroutine psm_update_dep_params_bilat( psm )
    
        type(t_psm), intent(inout)     :: psm
    
        real    :: dip, strike, rupdir, rake
        real, dimension(3) :: pol
   
      ! set up rotation matrices, using the eulerian angles
      ! strike=praezessionswinkel, dip=nutationswinkel, rupdir=drehwinkel
        strike = d2r(psm%params(6))
        dip = d2r(psm%params(7))
        rake = d2r(psm%params(8))
        rupdir = d2r(psm%params(9))
        
        call init_euler(dip,strike,-rupdir, psm%rotmat_rup)
        call init_euler(dip,strike,-rake, psm%rotmat_slip)
        
      ! p- and t- axis
        pol = r2d(domeshot(polar(matmul(psm%rotmat_slip,(/ sqrt(2.),0.,-sqrt(2.)/)))))
        psm%pax(1:2) = pol(2:3)
        pol = r2d(domeshot(polar(matmul(psm%rotmat_slip,(/-sqrt(2.),0.,-sqrt(2.)/)))))
        psm%tax(1:2) = pol(2:3)

    end subroutine
    
    subroutine psm_to_tdsm_bilat( psm, tdsm, shortest_doi, ok )
    
      ! translate a specific psm to tdsm
      ! shortest_duration is the shortest duration of interest
      ! automatically determines shape of output grid
      
        type(t_psm), intent(inout) :: psm
        type(t_tdsm), intent(out) :: tdsm
        real, intent(in) :: shortest_doi
        logical, intent(out) :: ok
        
        real :: rupvel
        real :: maxdx, maxdy, maxdt
        integer :: nx, ny, nt
        
        ok = .true.
        rupvel = psm%params(13)
                
        maxdt = shortest_doi
        maxdx = 0.5 * shortest_doi * rupvel
        maxdy = shortest_doi * rupvel
        call psm_to_tdsm_size_bilat( psm, maxdx, maxdy, maxdt, nx, ny, nt )
    
      ! calculate centroid moment tensor density
        call psm_to_tdsm_table_bilat( psm, tdsm, nx,ny,nt ) 
        
        psm%grid_size(1) = nx
        psm%grid_size(2) = ny
        psm%grid_size(3) = nt
    
    end subroutine
    
  
    subroutine psm_to_tdsm_size_bilat( in, maxdx, maxdy, maxdt, nx, ny, nt)
    
      ! given a psm, determine perfect number of grid points
      ! based on wanted maxdx, maxdy and maxdt
    
        type(t_psm), intent(inout) :: in
        real, intent(in) :: maxdx, maxdy, maxdt
        integer, intent(out) :: nx, ny, nt
  
        real dursf, durfull, length, width, rupvel, risetime, length_a, length_b
       
         
        length_a = in%params(10)
        length_b = in%params(11)
        width  = in%params(12)
        rupvel = in%params(13)
        risetime = in%params(14)

        length = length_a+length_b
        
        nx = floor(length/maxdx) + 1
      ! need at least 2 points to get partial deriv. depend on length
        if (nx .le. 1) nx = 2
      ! a line source may be made with width=0
        if (length .eq. 0.) nx = 1
        
        ny = floor(width/maxdy) + 1
      ! need at least 2 points to get partial deriv. depend on width
        if (ny .le. 1) ny = 2
      ! a point source (or simultaneous line source) may be made with length=0
        if (width .eq. 0.) ny = 1
        
        dursf = length/nx/rupvel ! duration of rupture front passing subfault
        durfull = risetime + dursf ! total duration of rupture on a subfault
        
      ! set temporal extension of source grid according to durfull
        nt = floor(durfull/maxdt) + 1
        
      ! need at least 2 points to get partial deriv. depend on duration
        if (nt .le. 1) nt = 2
  
    end subroutine
    
    
    subroutine psm_to_tdsm_table_bilat( psm, out, nx, ny, nt )
    
      ! translate parameterized source model psm psm
      ! to centroid table out on a grid of size nx * ny * nt
      ! nx, ny, nt should be determined by the psm_to_tdsm_size() sub
        
        type(t_psm), intent(inout)  :: psm
        type(t_tdsm), intent(out)   :: out
        integer, intent(in)             :: nx,ny,nt
         
      ! make some aliases to the parameters, so that the code get's readable
        real :: depth, moment, length, width, rupvel, risetime
        real :: north, east, length_a, length_b
        
        
        real, dimension(nx*ny)            :: tshift
        real, dimension(3,nx*ny)          :: grid
        real, dimension(nt)               :: wt, toff
        integer                           :: ix, iy, it
        integer                           :: ip, np, id, nd
        real                              :: dt
        real, dimension(3,3)              :: trotmat
        real, dimension(3)                :: p
        real, dimension(3,3)              :: m_rot
        real, dimension(3,3)              :: m_unrot = reshape((/0,0,-1,0,0,0,-1,0,0/),(/3,3/))
        type(t_plf)                       :: stf
        real                              :: dursf, durfull, tbeg, ta, tb
        
        north = psm%params(2)
        east = psm%params(3)
        depth = psm%params(4)
        moment = psm%params(5)
        
        length_a = psm%params(10)
        length_b = psm%params(11)
        width  = psm%params(12)
        rupvel = psm%params(13)
        risetime = psm%params(14)
        
        length = length_a + length_b
        
      ! set up grid.
      ! it is first layed out in the x-y plane; rupture direction is y
      ! later it is rotated to the desired orientation
              
      ! set up spacial grid, allocate temporary array to hold the coordinates
        
        np = nx*ny
        
      ! lay out the grid, centered in the x-y plane
      ! let it start at tshift=0 at its bottom border
      ! must work for width=0 or length=0 for point/line-sources
      ! grid points are at the center of the subfaults
      
        call eikonal_grid_destroy( psm%cgrid )
        psm%cgrid%ndims(1) = nx
        psm%cgrid%ndims(2) = ny
        allocate( psm%cgrid%points(3,nx,ny) )  ! keep a copy of the stuff in here, so that
        allocate( psm%cgrid%times(nx,ny) )     ! we can output consistent output files...
        
        ip = 1
        do ix=1,nx
            do iy=1,ny
                grid(1,ip) = (2.*(ix-1.)-nx+1.)/(2.*nx) * length
                grid(2,ip) = (2.*(iy-1.)-ny+1.)/(2.*ny) * width
                grid(3,ip) = 0.
                tshift(ip) = abs(length/2. - length_b + grid(1,ip))/rupvel + psm%params(1) &
                             - max(length_a,length_b)/2./rupvel
                
              ! rotate and translate to desired depth
                p = matmul(psm%rotmat_rup,grid(:,ip))
                grid(1,ip)=p(1)+north; grid(2,ip)=p(2)+east; grid(3,ip)=p(3)+depth
                
              ! keep copy for output
                psm%cgrid%points(:,ix,iy) = grid(:,ip)
                psm%cgrid%times(ix,iy) = tshift(ip) + max(length_a,length_b)/2./rupvel
              
                ip = ip+1
            end do
        end do
      
      ! source time function is a convolution of a box of width risetime with 
      ! another box with width of duration of the subfault
      ! the integral shall be normalized to 1
      ! setup piecewise linear function to represent this:
        
        dursf = length/nx/rupvel ! duration of subfault
        if (risetime < dursf) then
            call plf_make( stf, (-dursf-risetime)/2., 0., &
                                (-dursf+risetime)/2., 1./dursf, &
                                (dursf-risetime)/2.,  1./dursf, &
                                (dursf+risetime)/2.,  0. )
        else
            call plf_make( stf, (-risetime-dursf)/2., 0., &
                                (-risetime+dursf)/2., 1./risetime, &
                                (risetime-dursf)/2.,  1./risetime, &
                                (risetime+dursf)/2.,  0. )
        end if
        
        durfull = dursf+risetime ! total time a subfault is rupturing
        
        tbeg = stf%f(1,1)
        dt = durfull/nt
        
      ! weight on each time interval
      ! and 
        do it=1,nt
            ta = tbeg+dt*(it-1)
            tb = tbeg+dt*it
            call plf_integrate_and_centroid( stf, ta, tb, wt(it), toff(it) )
        end do
        
      ! all grid points should get the same moment
        nd = np*nt
        if (allocated( out%centroids )) deallocate(out%centroids)
        allocate( out%centroids(nd) )
        
      ! all grid points should get the same moment
      ! the proper mt is made by rotating a moment tensor
      ! according to strike, dip and rake
        trotmat = transpose(psm%rotmat_slip)
        m_rot = matmul( psm%rotmat_slip, matmul( m_unrot, trotmat ) )
        m_rot(:,:) = m_rot(:,:) * moment / np        
        
      ! fill output arrays
        id = 1
        do ip=1,np
            do it=1,nt
                out%centroids(id)%north = grid(1,ip)
                out%centroids(id)%east = grid(2,ip)
                out%centroids(id)%depth = grid(3,ip)
                out%centroids(id)%time = tshift(ip) + toff(it)
                out%centroids(id)%m(1) = m_rot(1,1)*wt(it)
                out%centroids(id)%m(2) = m_rot(2,2)*wt(it)
                out%centroids(id)%m(3) = m_rot(3,3)*wt(it)
                out%centroids(id)%m(4) = m_rot(1,2)*wt(it)
                out%centroids(id)%m(5) = m_rot(1,3)*wt(it)
                out%centroids(id)%m(6) = m_rot(2,3)*wt(it)
                id = id+1
            end do
        end do
        
    end subroutine
    
    
    subroutine psm_write_info_file_bilat( psm, fn )
        
        type(t_psm), intent(in)  :: psm
        type(varying_string), intent(in) :: fn

        integer :: unit, ix, iy
        real    :: north, east, depth, length, width, length_a, length_b
        real :: l1,l2,l3, w1,w2
        
        north = psm%params(2)
        east = psm%params(3)
        depth = psm%params(4)
        length_a = psm%params(10)
        length_b = psm%params(11)
        length = length_a + length_b
        width  = psm%params(12)
       
        call claim_unit( unit )
        open( unit=unit, file=char(fn), status='unknown' )
        
        l1 = -length/2.
        l2 = -length/2. + length_b
        l3 = length/2.
        
        w1 = -width/2.
        w2 = width/2.
        
        write (unit,"(a)") "origin"
        write (unit,*) psm%origin%lat, psm%origin%lon
        write (unit,*)
        
        write (unit,"(a)") "center"
        write (unit,*) north, east, depth
        write (unit,*)
        write (unit,"(a)") "outline"
        write (unit,*) matmul(psm%rotmat_rup,(/l2,w1,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l2,w2,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l3,w2,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l3,w1,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l1,w1,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l1,w2,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l2,w2,0./))+(/north,east,depth/)
        
        write (unit,*)
        write (unit,"(a)") "rupture"
        write (unit,*) matmul(psm%rotmat_rup,(/l2,0.,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l3,0.,0./)) & 
                       - matmul(psm%rotmat_rup,(/l2,0.,0./))
        write (unit,*) matmul(psm%rotmat_rup,(/l2,0.,0./))+(/north,east,depth/)
        write (unit,*) matmul(psm%rotmat_rup,(/l1,0.,0./)) & 
                       - matmul(psm%rotmat_rup,(/l2,0.,0./))
        write (unit,*)
        
        write (unit,"(a)") "slip"
      ! at position
        write (unit,*) matmul(psm%rotmat_slip,(/0.,0.,-200./))+(/north, east, depth/)
      ! slip vector
        write (unit,*) matmul(psm%rotmat_slip,(/1000.,0.,0./))
        ! at position
        write (unit,*) matmul(psm%rotmat_slip,(/0.,0.,200./))+(/north, east, depth/)
      ! slip vector
        write (unit,*) matmul(psm%rotmat_slip,(/-1000.,0.,0./))
        write (unit,*)
        
        write (unit,"(a)") "p-axis"
        write (unit,*) psm%pax(:)
        write (unit,*)
      
        write (unit,"(a)") "t-axis"
        write (unit,*) psm%tax(:)
        write (unit,*)
        
        write (unit,"(a)") "eikonal-grid"
        write (unit,*) psm%cgrid%ndims(:)
        do iy=1,size(psm%cgrid%times,2)
            do ix=1,size(psm%cgrid%times,1)
                write (unit,*) psm%cgrid%points(:,ix,iy), psm%cgrid%times(ix,iy)
            end do
        end do
        write (unit,*) 
  
        
        close( unit ) 
        call release_unit( unit )
    
    end subroutine
    
    function polar( xyz )
        real, dimension(3), intent(in)  :: xyz
        real, dimension(3) :: polar
        
        polar(1) = sqrt(dot_product(xyz,xyz))
        polar(2) = atan2(xyz(2),xyz(1))
        polar(3) = acos(xyz(3)/polar(1))
        
    end function
    
    function domeshot( polar )
    
        real, dimension(3), intent(in)  :: polar
        real, dimension(3) :: domeshot

        domeshot(1) = polar(1)
        domeshot(2:3) = wrap(polar(2:3),pi,-pi)
        if (domeshot(3) > pi/2.) then
            domeshot(2) = wrap(domeshot(2)+pi,-pi,pi)
            domeshot(3) = pi-domeshot(3)
        end if
        
    end function
    
    elemental function wrap( x, mi, ma )
        real, intent(in) :: x, mi, ma
        real :: wrap
        wrap = x - floor((x-mi)/(ma-mi)) * (ma-mi)
    end function
    
end module
