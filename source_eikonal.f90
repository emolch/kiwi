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

module source_eikonal

    use orthodrome
    use geometry
    use parameterized_source
    use discrete_source
    use piecewise_linear_function
    use euler
    use better_varying_string
    use unit
    use util
    use constants
    use crust2x2
    use eikonal
    
    implicit none 

    private
    
    integer, public, parameter :: n_source_params_eikonal = 15
    integer, public, parameter :: n_grid_dims_eikonal = 2
    
    public psm_to_tdsm_eikonal, psm_write_info_file_eikonal, psm_set_eikonal
    public psm_get_param_name_eikonal, psm_get_param_unit_eikonal, psm_get_param_id_eikonal
    public psm_cleanup_eikonal
    
    private psm_to_tdsm_size_eikonal, psm_to_tdsm_table_eikonal
    
    type(varying_string), private, dimension(:), allocatable :: psm_param_names_eikonal, psm_param_units_eikonal
    
  ! these values may be used to hint step sizes and the like
    real, dimension(n_source_params_eikonal), public, parameter :: psm_params_norm_eikonal = &
    (/ 1., 10000., 10000., 10000., 7e18, 360., 90., 360.,   10000., 10000., 10000.,  360., 10000.,  1.,1. /)
    
    real, private, parameter :: big = huge(psm_params_norm_eikonal(1))
    
  ! logical, physical or computational accuracy limits for params
    real, dimension(n_source_params_eikonal), public, parameter :: psm_params_min_hard_eikonal = &
    (/ -big, -100000., -100000.,        0., 1., -big, -big, -big, -1e7, -1e7, 0., -1e7, -1e7, 0.1, 0./)
    real, dimension(n_source_params_eikonal), public, parameter :: psm_params_max_hard_eikonal = &
    (/ big, 100000.,    100000., 1000000., 7e25, big, big, big,    1e7, 1e7, 1e7, 1e7, 1e7, 10., 10./)
    
  ! these limits give range of realistic, non-redundant or practical params
    real, dimension(n_source_params_eikonal), public, parameter :: psm_params_min_soft_eikonal = &
    (/ -20., -10000., -10000., 0., 1., -180., 0., -180., -100000., -100000., 0., -100000., -100000., 0.5, 0./)
    real, dimension(n_source_params_eikonal), public, parameter :: psm_params_max_soft_eikonal = &
    (/ 20., 10000., 10000., 150000., 7e25, 180., 90., 180., 100000., 100000., 100000., 100000., 100000., 1.5, 5./)
    
  ! defaults for external use
    real, dimension(n_source_params_eikonal), public, parameter :: psm_params_default_eikonal = &
    (/ 0.,0., 0., 3000.,    7e18,    0., 80., 0.,    0., 0., 5000., 0., 0., 0.9, 1./)
    
   ! M0 = 7e18 Nm <=> Mw = 6.5
   
        
  contains
  
    subroutine psm_cleanup_eikonal()
    
        integer :: i
        
        if (allocated( psm_param_names_eikonal )) then
            do i=1,n_source_params_eikonal
                call delete( psm_param_names_eikonal(i) )
            end do
            deallocate( psm_param_names_eikonal )
        end if
        if (allocated( psm_param_units_eikonal )) then
            do i=1,n_source_params_eikonal
                call delete( psm_param_units_eikonal(i) )
            end do
            deallocate( psm_param_units_eikonal )
        end if
        
    end subroutine
  
    subroutine psm_init_param_names_eikonal()
    
        if (.not. allocated( psm_param_names_eikonal )) then
            allocate( psm_param_names_eikonal( n_source_params_eikonal ) )
            psm_param_names_eikonal(1) = "time"
            psm_param_names_eikonal(2) = "north-shift"
            psm_param_names_eikonal(3) = "east-shift"
            psm_param_names_eikonal(4) = "depth"
            psm_param_names_eikonal(5) = "moment"
            psm_param_names_eikonal(6) = "strike"
            psm_param_names_eikonal(7) = "dip"
            psm_param_names_eikonal(8) = "slip-rake"
            
            psm_param_names_eikonal(9) = "bord-shift-x"    ! (rightward)
            psm_param_names_eikonal(10) = "bord-shift-y"   ! (downward)
            psm_param_names_eikonal(11) = "bord-radius"
            
            psm_param_names_eikonal(12) = "nukl-shift-x"   ! (rightward)
            psm_param_names_eikonal(13) = "nukl-shift-y"   ! (downward)
            
            psm_param_names_eikonal(14) = "rel-rupture-velocity"
            psm_param_names_eikonal(15) = "rise-time"
        end if
        
        if (.not. allocated( psm_param_units_eikonal )) then
            allocate( psm_param_units_eikonal( n_source_params_eikonal ) )
            psm_param_units_eikonal(1) = "s"
            psm_param_units_eikonal(2) = "m"
            psm_param_units_eikonal(3) = "m"
            psm_param_units_eikonal(4) = "m"
            psm_param_units_eikonal(5) = "Nm"
            psm_param_units_eikonal(6) = "degrees"
            psm_param_units_eikonal(7) = "degrees"
            psm_param_units_eikonal(8) = "degrees"
            
            psm_param_units_eikonal(9) = "m"
            psm_param_units_eikonal(10) = "m"
            psm_param_units_eikonal(11) = "m"
            
            psm_param_units_eikonal(12) = "m"
            psm_param_units_eikonal(13) = "m"
            
            psm_param_units_eikonal(14) = "1"
            psm_param_units_eikonal(15) = "s"
        end if
        
    end subroutine
  
    subroutine psm_get_param_name_eikonal( iparam, name )
    
        integer, intent(in)                 :: iparam
        type(varying_string), intent(inout) :: name
        
        call psm_init_param_names_eikonal()
        
        name = ""
        if (iparam < 1 .or. n_source_params_eikonal < iparam) return
        name = psm_param_names_eikonal(iparam)
        
    end subroutine
    
    subroutine psm_get_param_unit_eikonal( iparam, unit )
    
        integer, intent(in)                 :: iparam
        type(varying_string), intent(inout) :: unit
        
        call psm_init_param_names_eikonal()
        
        unit = ""
        if (iparam < 1 .or. n_source_params_eikonal < iparam) return
        unit = psm_param_units_eikonal(iparam)
        
    end subroutine
   
    subroutine psm_get_param_id_eikonal( name, iparam )
    
        type(varying_string), intent(in)     :: name
        integer, intent(out)                 :: iparam

        integer :: i
        
        call psm_init_param_names_eikonal()
        iparam = 0
        do i=1, n_source_params_eikonal
            if (name .eq. psm_param_names_eikonal(i)) then
                iparam = i
                return
            end if
        end do
    
    end subroutine
  
    subroutine psm_set_eikonal( psm,  params, normalized, only_moment_changed)
    
        type(t_psm), intent(inout)     :: psm
        real, dimension(:), intent(in) :: params
        logical, intent(in)            :: normalized
        logical, intent(out)           :: only_moment_changed

        real, dimension(size(params)) :: new_params
        logical :: must_reset_grid
        integer :: ip
        
        must_reset_grid = .false.
        if (.not. allocated(psm%grid_size) .or. (psm%sourcetype .ne. psm_eikonal)) then
            must_reset_grid = .true.
        end if
        
        call resize( psm%params, 1, n_source_params_eikonal )
        call resize( psm%grid_size, 1, n_grid_dims_eikonal )
        call resize( psm%params_norm, 1, n_source_params_eikonal )
        
        psm%params_norm(:) = psm_params_norm_eikonal(:)
        
        if (must_reset_grid) then
            psm%grid_size = 1
        end if
        
        if (size(params,1) .ne. size(psm%params)) call die("wrong number of source parameters in psm_set_eikonal()")
                
        if (normalized) then
            new_params = params * psm%params_norm
        else
            new_params = params 
        end if
        
        only_moment_changed = .true.
        do ip=1,size(psm%params)
            if (ip .ne. 5 .and. ip .ne. 15 .and. psm%params(ip) .ne. new_params(ip)) then
                only_moment_changed = .false.
            end if
        end do

        psm%params = new_params
        
        psm%moment = psm%params(5)
        psm%risetime = psm%params(15)
        
        call psm_update_dep_params_eikonal( psm )
        
    end subroutine
        
    subroutine psm_update_dep_params_eikonal( psm )
    
        type(t_psm), intent(inout)     :: psm
    
        real    :: dip, strike, rake
        real, dimension(3) :: pol
   
      ! set up rotation matrices, using the eulerian angles
      ! strike=praezessionswinkel, dip=nutationswinkel, rupdir=drehwinkel
        strike = d2r(psm%params(6))
        dip = d2r(psm%params(7))
        rake = d2r(psm%params(8))
        
        call init_euler(dip, strike, -rake, psm%rotmat_slip)
        call init_euler(dip, strike, 0., psm%rotmat_rup)
        
      ! p- and t- axis
        pol = r2d(domeshot(polar(matmul(psm%rotmat_slip,(/ sqrt(2.),0.,-sqrt(2.)/)))))
        psm%pax(1:2) = pol(2:3)
        pol = r2d(domeshot(polar(matmul(psm%rotmat_slip,(/-sqrt(2.),0.,-sqrt(2.)/)))))
        psm%tax(1:2) = pol(2:3)

    end subroutine
    
    subroutine psm_to_tdsm_eikonal( psm, tdsm, shortest_doi, ok )
    
      ! translate a specific psm to tdsm
      ! shortest_duration is the shortest duration of interest
      ! automatically determines shape of output grid
            
        type(t_psm), intent(inout) :: psm
        type(t_tdsm), intent(out) :: tdsm
        real, intent(in) :: shortest_doi
        logical, intent(out) :: ok
        
        real :: maxdx, maxdy, maxdt, bord_shift_x, bord_shift_y, bord_radius
        integer :: nx, ny
        type(t_polygon) :: rupture_poly, rupture_poly_rc
        real, dimension(3) :: min_rc, max_rc
        real :: deltagrid

        bord_shift_x = psm%params(9)
        bord_shift_y = psm%params(10)
        bord_radius = psm%params(11)
      
        ok = .true.
      ! determine region of interest
        call psm_borderline_eikonal( psm, bord_shift_x, bord_shift_y, bord_radius, rupture_poly, rupture_poly_rc )
        
      ! fail if rupture area is completely eaten up by constraints
        if (size(rupture_poly%points,2) == 0) then
            call error("Empty rupture area")
            ok = .false.
            return
        end if

      ! solve eikonal equation on a fine and rectangular grid
        call polygon_box( rupture_poly_rc, min_rc, max_rc )
        deltagrid = min(100.*shortest_doi/2.,4000.)
        call psm_make_eikonal_grid( psm, min_rc, max_rc, deltagrid, psm%egrid, ok )
        if (.not. ok) return
        
      ! determine optimal grid size  
        maxdt = shortest_doi
        maxdx = 0.5 * shortest_doi * psm%egrid%minspeed
        maxdy = 0.5 * shortest_doi * psm%egrid%minspeed
        
        call psm_to_tdsm_size_eikonal( psm%egrid%last(1)-psm%egrid%first(1), &
                                       psm%egrid%last(2)-psm%egrid%first(2), &
                                           maxdx, maxdy, nx, ny)

      ! downsample the fine grid to the optimal grid
        call psm_downsample_grid( psm, psm%egrid, nx, ny, psm%cgrid )
        
      ! calculate table of centroids approximating moment tensor density
        call psm_to_tdsm_table_eikonal( psm, psm%cgrid, maxdt, tdsm ) 
        
        
        psm%grid_size(1) = nx
        psm%grid_size(2) = ny
    
    end subroutine
    
    subroutine psm_borderline_eikonal( psm, bord_shift_x, bord_shift_y, bord_radius, &
                                       rupture_poly, rupture_poly_rc )
    
      ! determine polygon which borders the source area

        type(t_psm), intent(in) :: psm
        real, intent(in) :: bord_radius, bord_shift_x, bord_shift_y
        
        type(t_polygon), intent(inout) :: rupture_poly, rupture_poly_rc
        type(t_polygon) :: circle_poly
        type(t_circle) :: circle
        
        integer :: n_initial_points, ipoint
        
      ! setup initial border circle
        circle%center = psm_circle_center( psm, bord_shift_x, bord_shift_y )
        circle%transform = -psm%rotmat_rup * bord_radius
        
      ! cut out parts of circle which lie outside constraint halfspaces
        n_initial_points = 180
        if (bord_radius == 0.) n_initial_points = 1
        call circle_to_polygon( circle, n_initial_points, circle_poly )
        call trim_polygon( circle_poly, psm%constraints, rupture_poly )
        
      ! get the border polygon in rupture coordinates
        call copy_polygon( rupture_poly, rupture_poly_rc )
        do ipoint=1,size(rupture_poly_rc%points,2)
            rupture_poly_rc%points(:,ipoint) = psm_ned_to_rc( psm, rupture_poly%points(:,ipoint) )
        end do
        
    end subroutine
    
    function psm_circle_center( psm, bord_shift_x, bord_shift_y )
    
      ! get center of bounding circle in ned coordinates

        type(t_psm), intent(in) :: psm
        real, intent(in) :: bord_shift_x, bord_shift_y
        real, dimension(3) :: psm_circle_center
        
        psm_circle_center = psm_rc_to_ned(psm, (/bord_shift_x, bord_shift_y ,0./))
    
    end function
    
    function psm_initial_point_rc( psm, borderline )

      ! Get a possible nucleation point which lies inside ruputure area.
      ! Coordinates returned are in rupture coordinates.

      ! The nucleation point which is set via psm%params is first projected
      ! to the surrounding circle, and then from there to the closest point on
      ! the polygon.
    
        type(t_psm), intent(in)     :: psm
        type(t_polygon), intent(in) :: borderline
        real, dimension(3)          :: psm_initial_point_rc
        
        real, dimension(3)          :: initial_point_ned
        real                        :: bord_radius
        real                        :: nukl_shift_x, nukl_shift_y
        real                        :: nukl_shift
        
      
        bord_radius = psm%params(11)
        nukl_shift_x = psm%params(12)
        nukl_shift_y = psm%params(13)
        
        nukl_shift = sqrt(nukl_shift_x**2+nukl_shift_y**2)
        
        if (nukl_shift > bord_radius) then
            nukl_shift_x = nukl_shift_x * bord_radius/nukl_shift
            nukl_shift_y = nukl_shift_y * bord_radius/nukl_shift
        end if
        
        psm_initial_point_rc = (/nukl_shift_x, nukl_shift_y, 0./)
        
        initial_point_ned = psm_rc_to_ned( psm, psm_initial_point_rc )
        if (.not. psm_point_in_constraints(psm, initial_point_ned)) then
            psm_initial_point_rc = psm_ned_to_rc &
                ( psm, nearest_point_on_polygon(borderline, initial_point_ned) )
        end if
        
    end function

    function psm_initial_point_intolerant_rc( psm, ok )

      ! Get nucleation point.
      ! Indicates point outside of rupture area with ok=false.
      ! Coordinates returned are in rupture coordinates.
    
        type(t_psm), intent(in)     :: psm
        logical, intent(out)        :: ok
        real, dimension(3)          :: psm_initial_point_intolerant_rc
        
        real, dimension(3)          :: initial_point_ned
        real                        :: bord_radius
        real                        :: nukl_shift_x, nukl_shift_y
        real                        :: nukl_shift
        
      
        bord_radius = psm%params(11)
        nukl_shift_x = psm%params(12)
        nukl_shift_y = psm%params(13)
        
        nukl_shift = sqrt(nukl_shift_x**2+nukl_shift_y**2)
        
        ok = .true.
        psm_initial_point_intolerant_rc = (/nukl_shift_x, nukl_shift_y, 0./)
        initial_point_ned = psm_rc_to_ned( psm, psm_initial_point_intolerant_rc )
        if (.not. psm_point_in_constraints(psm, initial_point_ned) .or. nukl_shift > bord_radius) then
            call error("position of nucleation point is outside of rupture region")
            ok = .false.
        end if

    end function
    
    
    subroutine psm_make_eikonal_grid( psm, min_rc, max_rc, approx_delta, grid, ok )
    
        type(t_psm), intent(in) :: psm
        real, dimension(3), intent(in) :: min_rc, max_rc
        real, intent(in) :: approx_delta
        type(t_eikonal_grid), intent(inout) :: grid
        logical, intent(out) :: ok

        integer :: iy,ix
        real :: rel_rupture_velocity, bord_radius, vp, vs, rho, minspeed, invalid_speed
        real, dimension(3) :: circle_center, point, point_rc, initial_point_rc
        real, dimension(2) ::  dims
        type(t_crust2x2_1d_profile) :: profile
        real :: bord_shift_x, bord_shift_y

        ok = .true.
        
        call eikonal_grid_destroy( grid )
        
        bord_shift_x = psm%params(9)
        bord_shift_y = psm%params(10)
        bord_radius = psm%params(11)
        
        rel_rupture_velocity = psm%params(14)
        
        grid%first = min_rc(1:2)
        grid%last = max_rc(1:2)
        dims = grid%last - grid%first
        grid%ndims = int(ceiling(dims/approx_delta))
        if (grid%ndims(1) == 0) grid%ndims(1) = 1
        if (grid%ndims(2) == 0) grid%ndims(2) = 1
        grid%delta = dims/(grid%ndims)
        
        allocate(grid%speed(grid%ndims(1),grid%ndims(2)))
        allocate(grid%times(grid%ndims(1),grid%ndims(2)))
        allocate(grid%points(3,grid%ndims(1),grid%ndims(2)))
        
        call crust2x2_get_profile(psm%origin, profile)
        circle_center = psm_circle_center( psm, bord_shift_x, bord_shift_y )
        
        !initial_point_rc = psm_initial_point_rc( psm, borderline )
        ! XXX
        initial_point_rc = psm_initial_point_intolerant_rc( psm, ok )
        if (.not. ok) return
        ! XXX

        grid%initialpoint = initial_point_rc(1:2)
      
      ! setup speeds
        
        minspeed = huge(vs)
        do iy=1,grid%ndims(2)
            do ix=1,grid%ndims(1)
                point_rc(1) = grid%first(1) + (ix-0.5)*grid%delta(1)
                point_rc(2) = grid%first(2) + (iy-0.5)*grid%delta(2)
                point_rc(3) = 0.
                point = psm_rc_to_ned( psm, point_rc )
                grid%points(:,ix,iy) = point
                
                if (veclen(point-circle_center) > bord_radius .or. &
                    .not. psm_point_in_constraints( psm, point )) then 
                    grid%speed(ix,iy) = 0.
                else 
                    call crust2x2_get_at_depth( profile, point(3), vp, vs, rho )
                    grid%speed(ix,iy) = vs * rel_rupture_velocity
                    minspeed = min(grid%speed(ix,iy), minspeed)
                end if
            end do
        end do
        grid%minspeed = minspeed
        
      ! solve eikonal problem
        
      ! replace zero speed with a small value, to prevent stuff coming in from outside
        invalid_speed = minspeed*0.5
        where (grid%speed == 0.) grid%speed = invalid_speed
        
        call eikonal_solver_fmm( grid%speed, grid%first, grid%delta, grid%initialpoint, grid%times )
        
      ! flag points outside of valid region with an invalid time
        where (grid%speed == invalid_speed) grid%times = -1.
        
    end subroutine
    
    subroutine psm_downsample_grid( psm, fgrid, nxc,nyc, cgrid )

        type(t_psm), intent(in)              :: psm
        type(t_eikonal_grid), intent(in)     :: fgrid
        integer, intent(in)                  :: nxc, nyc
        type(t_eikonal_grid), intent(inout)  :: cgrid
 
        integer :: npf
        integer :: ixf, iyf
        integer :: ixc, iyc
        real, dimension(nxc,nyc) :: ntimes
        real, dimension(3) :: point_rc

        call eikonal_grid_destroy( cgrid )
        
        ntimes = 0
        cgrid%first = fgrid%first
        cgrid%last  = fgrid%last
        cgrid%ndims = (/nxc,nyc/)
        cgrid%delta = (cgrid%last-cgrid%first)/cgrid%ndims
        where (cgrid%delta == 0. .or. cgrid%ndims == 0) cgrid%delta = 1.
        
        allocate(cgrid%speed(nxc,nyc))
        allocate(cgrid%times(nxc,nyc))
        allocate(cgrid%durations(nxc,nyc))
        allocate(cgrid%weights(nxc,nyc))
        allocate(cgrid%points(3,nxc,nyc))

        cgrid%times = -1.
        cgrid%speed = 0.
        cgrid%points = 0.
        npf = 0
        do iyf=1,fgrid%ndims(2)
            do ixf=1,fgrid%ndims(1)
                if (fgrid%times(ixf,iyf) < 0.) cycle
                point_rc = psm_ned_to_rc( psm, fgrid%points(:,ixf,iyf) )
                ixc = floor((point_rc(1)-cgrid%first(1))/cgrid%delta(1))+1
                iyc = floor((point_rc(2)-cgrid%first(2))/cgrid%delta(2))+1
                if (any((/ixc,iyc/) < 1) .or. any((/ixc,iyc/) > (/nxc,nyc/))) then
                    call warn("orphaned point in fine grid: "//ixc //" "//iyc)
                    cycle
                end if
                ntimes(ixc,iyc) = ntimes(ixc,iyc) + 1
                if (cgrid%times(ixc,iyc) == -1.) cgrid%times(ixc,iyc) = 0.
                cgrid%times(ixc,iyc) = cgrid%times(ixc,iyc) + fgrid%times(ixf,iyf)
                cgrid%speed(ixc,iyc) = cgrid%speed(ixc,iyc) + 1./fgrid%speed(ixf,iyf)
    
    
                cgrid%points(:,ixc,iyc) = cgrid%points(:,ixc,iyc) + fgrid%points(:,ixf,iyf)
                npf = npf + 1
            end do
        end do
        where (ntimes > 0)
            cgrid%times = 1./ntimes * cgrid%times
            cgrid%speed = 1./(1./ntimes * cgrid%speed)
            cgrid%points(1,:,:) = 1./ntimes * cgrid%points(1,:,:)
            cgrid%points(2,:,:) = 1./ntimes * cgrid%points(2,:,:)
            cgrid%points(3,:,:) = 1./ntimes * cgrid%points(3,:,:)
        end where
        cgrid%weights = ntimes/float(npf)
        
      ! second pass to get the approximated durations right:
        cgrid%durations = 0.
        do iyf=1,fgrid%ndims(2)
            do ixf=1,fgrid%ndims(1)
                if (fgrid%times(ixf,iyf) < 0.) cycle
                point_rc = psm_ned_to_rc( psm, fgrid%points(:,ixf,iyf) )
                ixc = floor((point_rc(1)-cgrid%first(1))/cgrid%delta(1))+1
                iyc = floor((point_rc(2)-cgrid%first(2))/cgrid%delta(2))+1
                if (any((/ixc,iyc/) < 1) .or. any((/ixc,iyc/) > (/nxc,nyc/))) then
                    call warn("orphaned point in fine grid")
                    cycle
                end if
                
                cgrid%durations(ixc,iyc) = cgrid%durations(ixc,iyc) + &
                    abs(fgrid%times(ixf,iyf)-cgrid%times(ixc,iyc))
        
            end do
        end do
        where (ntimes > 0)
            cgrid%durations = 4./ntimes * cgrid%durations
        end where
    end subroutine
    
    pure function psm_ned_to_rc( psm, point ) result(point_rc)
        type(t_psm), intent(in)         :: psm
        real, dimension(3), intent(in)  :: point
        real, dimension(3)              :: point_rc
        point_rc = matmul(transpose(psm%rotmat_rup),point- psm%params(2:4) )
    end function
    
    pure function psm_rc_to_ned( psm, point_rc ) result(point)
        type(t_psm), intent(in)         :: psm
        real, dimension(3), intent(in)  :: point_rc
        real, dimension(3)              :: point
        point = matmul(psm%rotmat_rup, point_rc) + psm%params(2:4)
    end function
    
    subroutine psm_to_tdsm_size_eikonal( sizex, sizey, maxdx, maxdy, nx, ny)
    
      ! given a psm, determine perfect number of grid points
      ! based on wanted maxdx, maxdy and maxdt
    
        real, intent(in) :: sizex,sizey
        real, intent(in) :: maxdx, maxdy
        integer, intent(out) :: nx, ny
  
        nx = floor(sizex/maxdx) + 1
      ! need at least 2 points to get partial deriv. depend on length
        if (nx .le. 1) nx = 2
      ! a line source may be made with width=0
        if (sizex .eq. 0.) nx = 1
        
        ny = floor(sizey/maxdy) + 1
      ! need at least 2 points to get partial deriv. depend on width
        if (ny .le. 1) ny = 2
      ! a point source (or simultaneous line source) may be made with length=0
        if (sizey .eq. 0.) ny = 1
        
    end subroutine
    
    subroutine psm_to_tdsm_table_eikonal( psm, grid, maxdt, out )
    
      ! Translate parameterized source model to centroid table

        type(t_psm), intent(inout)        :: psm
        type(t_eikonal_grid), intent(in)  :: grid
        real, intent(in)                  :: maxdt
        type(t_tdsm), intent(inout)       :: out
         
        real :: origin_time
        
        real, dimension(:), allocatable   :: tweights, toffsets
        integer                           :: ix, iy
        integer                           :: it, nt
        integer                           :: id, nd
        real, dimension(3,3)              :: trotmat
        real, dimension(3,3)              :: m_rot
        real, dimension(3,3)              :: m_unrot = reshape((/0,0,-1,0,0,0,-1,0,0/),(/3,3/))
        real :: centertime
        
        origin_time = psm%params(1)
       
      ! determine needed size for the centroid table
        nd = 0
        centertime = 0.
        do iy=1,grid%ndims(2)
            do ix=1,grid%ndims(1)
                if (grid%times(ix,iy) >= 0.) then
                    nd = nd + int(floor((grid%durations(ix,iy)+0.0)/maxdt))+1
                    centertime = centertime + grid%times(ix,iy)*grid%weights(ix,iy)
                end if
            end do
        end do
        
        if (allocated( out%centroids )) deallocate(out%centroids)
        allocate( out%centroids(nd) )
        
      ! the proper mt is made by rotating moment tensor
      ! according to strike, dip and rake
        trotmat = transpose(psm%rotmat_slip)
        m_rot = matmul( psm%rotmat_slip, matmul( m_unrot, trotmat ) )
        m_rot(:,:) = m_rot(:,:) 

      ! fill output arrays
        id = 1
        do iy=1,grid%ndims(2)
            do ix=1,grid%ndims(1)
                if (grid%times(ix,iy) < 0.) cycle
                
              ! using a zero risetime here, effect of (constant) risetime applied after synthesis of seismogram
                call discretize_subfault_time &
                   ( grid%durations(ix,iy), 0., maxdt, tweights, toffsets, nt )

                do it=1,nt
                    out%centroids(id)%north = grid%points(1,ix,iy)
                    out%centroids(id)%east = grid%points(2,ix,iy)
                    out%centroids(id)%depth = grid%points(3,ix,iy)
                    out%centroids(id)%time = grid%times(ix,iy) + toffsets(it) + origin_time - centertime
                   
                    out%centroids(id)%m(1) = m_rot(1,1)*tweights(it) * grid%weights(ix,iy)
                    out%centroids(id)%m(2) = m_rot(2,2)*tweights(it) * grid%weights(ix,iy)
                    out%centroids(id)%m(3) = m_rot(3,3)*tweights(it) * grid%weights(ix,iy)
                    out%centroids(id)%m(4) = m_rot(1,2)*tweights(it) * grid%weights(ix,iy)
                    out%centroids(id)%m(5) = m_rot(1,3)*tweights(it) * grid%weights(ix,iy)
                    out%centroids(id)%m(6) = m_rot(2,3)*tweights(it) * grid%weights(ix,iy)
                    id = id+1
                end do
            end do
        end do
        if (allocated(tweights)) deallocate(tweights)
        if (allocated(toffsets)) deallocate(toffsets)        

    end subroutine
    
    subroutine discretize_subfault_time( duration_subfault, risetime, maxdt, tweights, toffsets, nt )
    
        real, intent(in) :: duration_subfault, risetime, maxdt
        real, dimension(:), allocatable, intent(inout) :: tweights, toffsets
        integer, intent(out) :: nt
        type(t_plf) :: stf
        
        real :: durfull, dursf, ta, tb, tbeg,dt
        integer :: it
 
        dursf = duration_subfault
        
        
        durfull = dursf+risetime ! total time a subfault is rupturing
        
        nt = int(floor(durfull/maxdt))+1

        if (.not. allocated(tweights)) allocate(tweights(nt))
        if (.not. allocated(toffsets)) allocate(toffsets(nt))
        if (nt > size(tweights)) call resize( tweights, 1, nt )
        if (nt > size(toffsets)) call resize( toffsets, 1, nt )
        
        if (nt == 1) then
            tweights(1) = 1.
            toffsets(1) = 0.
            return
        end if

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
        tbeg = stf%f(1,1)
        dt = durfull/nt
        
        do it=1,nt
            ta = tbeg+dt*(it-1)
            tb = tbeg+dt*it
            call plf_integrate_and_centroid( stf, ta, tb, tweights(it), toffsets(it) )
        end do
    
        call plf_destroy( stf )
        
    end subroutine
    
    
    subroutine psm_write_info_file_eikonal( psm, fn )
        
        type(t_psm), intent(in)  :: psm
        type(varying_string), intent(in) :: fn

        integer :: unit, i, ix, iy
        real    :: north, east, depth, bord_shift_x, bord_shift_y, bord_radius
        type(t_polygon) :: rupture_poly, rupture_poly_rc
        type(t_crust2x2_1d_profile) :: profile
        real    :: vp, vs, rho
        real, dimension(3)    :: point_rc, initial_point_rc, initial_point
        logical :: ok
        
        north = psm%params(2)
        east = psm%params(3)
        depth = psm%params(4)        
        
        bord_shift_x = psm%params(9)
        bord_shift_y = psm%params(10)
        bord_radius = psm%params(11)
        
        call psm_borderline_eikonal( psm, bord_shift_x, bord_shift_y, bord_radius, rupture_poly, rupture_poly_rc )
        
        call claim_unit( unit )
        open( unit=unit, file=char(fn), status='unknown' )
           
        write (unit,"(a)") "origin"
        write (unit,*) psm%origin%lat, psm%origin%lon
        write (unit,*)
                
        
        write (unit,"(a)") "center"
        write (unit,*) north, east, depth
        write (unit,*)
                
        initial_point_rc = psm_initial_point_intolerant_rc( psm, ok )
        initial_point = psm_rc_to_ned( psm, initial_point_rc )
        write (unit,"(a)") "nucleation-point"
        write (unit,*) initial_point, initial_point_rc(1:2)
        write (unit,*) 
        
        write (unit,"(a)") "outline"
        do i=1,size(rupture_poly%points,2)
            write (unit,*) rupture_poly%points(:,i), rupture_poly_rc%points(1:2,i)
        end do
        write (unit,*)
        
        write (unit,"(a)") "rupture"
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

        call crust2x2_get_profile(psm%origin, profile)

        write (unit,"(a)") "eikonal-grid"
        write (unit,*) psm%cgrid%ndims(:)
        do iy=1,size(psm%cgrid%times,2)
            do ix=1,size(psm%cgrid%times,1)
                call crust2x2_get_at_depth( profile, psm%cgrid%points(3,ix,iy), vp, vs, rho )
                point_rc(:) = psm_ned_to_rc( psm, psm%cgrid%points(:,ix,iy) )
                write (unit,*) psm%cgrid%points(:,ix,iy), point_rc(1:2), psm%cgrid%times(ix,iy), vp, vs, rho
            end do
        end do
        write (unit,*) 
        
        write (unit,"(a)") "area"
        write (unit,*) polygon_area(rupture_poly)
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
    
    pure function veclen( a )
        real, dimension(3), intent(in)  :: a
        real :: veclen
        veclen = sqrt( dot_product(a,a) )
    end function
    
end module
