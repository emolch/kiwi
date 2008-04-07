! ! 
!    Copyright 2007 Lars Krieger, Sebastian Heimann
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

! ------------------------------------------------------------------------------

! code for a spacial point source with a temporal oscillatory LP-conditions 
! 
! Main period of excitation and temporal extension are given as parameters "prd" and " dur_exc"


module source_point_lp

    use constants
    use orthodrome
    use parameterized_source
    use discrete_source
    use piecewise_linear_function
    use euler
    use better_varying_string
    use unit
    use util
    
    !---------------------------------------
    implicit none 
    !---------------------------------------

    private


    integer, public, parameter :: n_source_params_point_lp = 13
    integer, public, parameter :: n_grid_dims_point_lp = 3
    
    public psm_to_tdsm_point_lp, psm_write_info_file_point_lp, psm_set_point_lp
    public psm_get_param_name_point_lp, psm_get_param_unit_point_lp, psm_get_param_id_point_lp
    public psm_cleanup_point_lp
    
    private psm_to_tdsm_table_point_lp
    
    type(varying_string), private, dimension(:), allocatable :: psm_param_names_point_lp, psm_param_units_point_lp
    
    real, dimension(n_source_params_point_lp), public, parameter :: psm_params_norm_point_lp = &
    (/ 1., 10000., 10000., 10000., 7e18, 1., 0., -1., 1., 1.,1.,20.,1. /)
    
    real, private, parameter :: big = huge(psm_params_norm_point_lp(1))   
    
  ! logical, physical or computational accuracy limits for params
    real, dimension(n_source_params_point_lp), public, parameter :: psm_params_min_hard_point_lp = &
    (/ -big, -100000., -100000., 0., 1., -1000., -1000., -1000., -1000., -1000., -1000., 0., 0./)
    real, dimension(n_source_params_point_lp), public, parameter :: psm_params_max_hard_point_lp = &
    (/ big, 100000., 100000., 1000000., 7e25, 1000., 1000., 1000., 1000., 1000., 1000., 120., 120./)
    
  ! these limits give range of realistic, non-redundant or practical params
    real, dimension(n_source_params_point_lp), public, parameter :: psm_params_min_soft_point_lp = &
    (/ -big, -10000., -10000., 0., 1., -100., -100., -100., -100., -100., -100., 0.001, 0.001/)
    real, dimension(n_source_params_point_lp), public, parameter :: psm_params_max_soft_point_lp = &
    (/ big, 10000., 10000., 150000., 7e24, 100., 100., 100., 100., 100., 100., 90., 50./)
    
  ! defaults for external use
    real, dimension(n_source_params_point_lp), public, parameter :: psm_params_default_point_lp = &
    (/ 0.,0., 0., 10000., 7e18, 0., -2., 2., 9., 0., -1., 40., 1. /)
    
  ! M0 = 7e18 Nm <=> Mw = 6.5

    !---------------------------------------
    !---------------------------------------



  contains
    
    !---------------------------------------
    
    subroutine psm_cleanup_point_lp()
      
      integer :: i
      
      if (allocated( psm_param_names_point_lp )) then
         do i=1,n_source_params_point_lp
            call delete( psm_param_names_point_lp(i) )
         end do
         deallocate( psm_param_names_point_lp )
      end if
      if (allocated( psm_param_units_point_lp )) then
         do i=1,n_source_params_point_lp
            call delete( psm_param_units_point_lp(i) )
         end do
         deallocate( psm_param_units_point_lp )
      end if
      
    end subroutine psm_cleanup_point_lp
    !---------------------------------------
   
    subroutine psm_init_param_names_point_lp()
      
      if (.not. allocated( psm_param_names_point_lp )) then
         allocate( psm_param_names_point_lp( n_source_params_point_lp ) )
         psm_param_names_point_lp(1) = "time"
         psm_param_names_point_lp(2) = "north-shift"
         psm_param_names_point_lp(3) = "east-shift"
         psm_param_names_point_lp(4) = "depth"
         psm_param_names_point_lp(5) = "moment"
         psm_param_names_point_lp(6) = "m_xx"
         psm_param_names_point_lp(7) = "m_yy"
         psm_param_names_point_lp(8) = "m_zz"
         psm_param_names_point_lp(9) = "m_xy"
         psm_param_names_point_lp(10) = "m_xz"
         psm_param_names_point_lp(11) = "m_yz"
         psm_param_names_point_lp(12) = "excitation-time"
         psm_param_names_point_lp(13) = "main-period"
      end if
      
      if (.not. allocated( psm_param_units_point_lp )) then
         allocate( psm_param_units_point_lp( n_source_params_point_lp ) )
         psm_param_units_point_lp(1) = "s"
         psm_param_units_point_lp(2) = "m"
         psm_param_units_point_lp(3) = "m"
         psm_param_units_point_lp(4) = "m"
         psm_param_units_point_lp(5) = "Nm"
         psm_param_units_point_lp(6) = "Nm"
         psm_param_units_point_lp(7) = "Nm"
         psm_param_units_point_lp(8) = "Nm"
         psm_param_units_point_lp(9) = "Nm"
         psm_param_units_point_lp(10)= "Nm"
         psm_param_units_point_lp(11)= "Nm"
         psm_param_units_point_lp(12)= "s"
         psm_param_units_point_lp(13)= "s"
      end if
      
    end subroutine psm_init_param_names_point_lp
    !---------------------------------------
    
    subroutine psm_get_param_name_point_lp( iparam, name )
      
      integer, intent(in)                 :: iparam
      type(varying_string), intent(inout) :: name
      
      call psm_init_param_names_point_lp()
      
      name = ""
      if (iparam < 1 .or. n_source_params_point_lp < iparam) return
      name = psm_param_names_point_lp(iparam)
      
    end subroutine psm_get_param_name_point_lp
    !---------------------------------------
    
    subroutine psm_get_param_unit_point_lp( iparam, unit )
      
      integer, intent(in)                 :: iparam
      type(varying_string), intent(inout) :: unit
      
      call psm_init_param_names_point_lp()
      
      unit = ""
      if (iparam < 1 .or. n_source_params_point_lp < iparam) return
      unit = psm_param_units_point_lp(iparam)
      
    end subroutine psm_get_param_unit_point_lp
    !---------------------------------------
    
    subroutine psm_get_param_id_point_lp( name, iparam )
      
      type(varying_string), intent(in)     :: name
      integer, intent(out)                 :: iparam
      
      integer :: i
      
      call psm_init_param_names_point_lp()
      iparam = 0
      do i=1, n_source_params_point_lp
         if (name .eq. psm_param_names_point_lp(i)) then
            iparam = i
            return
         end if
      end do
      
    end subroutine psm_get_param_id_point_lp
    !---------------------------------------
    
    subroutine psm_set_point_lp( psm,  params, normalized_ )
      
      type(t_psm), intent(inout)     :: psm
      real, dimension(:), intent(in) :: params
      logical, intent(in), optional  :: normalized_
      
      logical :: must_reset_grid
      logical :: normalized
      
      normalized = .false.
      if (present(normalized_)) normalized = normalized_
        
      must_reset_grid = .false.
      if (.not. allocated(psm%grid_size) .or. (psm%sourcetype .ne. psm_point_lp)) then
         must_reset_grid = .true.
      end if
      
      call resize( psm%params, 1, n_source_params_point_lp )
      call resize( psm%grid_size, 1, n_grid_dims_point_lp )
      call resize( psm%params_norm, 1, n_source_params_point_lp )
      
      psm%params_norm(:) = psm_params_norm_point_lp(:)
      
      if (must_reset_grid) then
         psm%grid_size = 1
      end if
      
      if (size(params,1) .ne. size(psm%params)) call die("wrong number of source parameters in psm_set_point_lp()")
      
      if (normalized) then
         psm%params = params * psm%params_norm
      else
         psm%params = params 
      end if
      
      
    end subroutine psm_set_point_lp
    !---------------------------------------
    
      
    subroutine psm_to_tdsm_point_lp( psm, tdsm, shortest_doi )
  
  
      ! translate a specific psm to tdsm
      ! shortest_duration is minimum of the shortest duration of interest and the conditions of nyquist (period/2)
      ! automatically determines shape of output grid
      
      type(t_psm), intent(inout) :: psm
      type(t_tdsm), intent(out) :: tdsm
      real, intent(in) :: shortest_doi
      
      real :: maxdt
      integer :: nx, ny, nt
      real :: dur_exc
      
      maxdt = shortest_doi

      dur_exc = psm%params(12)
      
      nx = 1   ! no extention in space
      ny = nx  ! spacial symmetry
      
      ! set temporal extension of source grid according to total duration of excitation
      nt = floor(dur_exc/maxdt) + 1
      
      ! need at least 2 points to get partial deriv. depend on duration
      if (nt .le. 1) nt = 2
      
      ! calculate centroid moment tensor density
      call psm_to_tdsm_table_point_lp( psm, tdsm, maxdt ,nt ) 
      
      psm%grid_size(1) = nx
      psm%grid_size(2) = ny
      psm%grid_size(3) = nt
      
    end subroutine psm_to_tdsm_point_lp
    
!!! new.....alpha-version.....test !!!!!
    
    
    subroutine psm_to_tdsm_table_point_lp( psm, out, maxdt, nt )
      
      ! translate parameterized source model psm 
      ! to centroid table out on a grid of size nx * ny * nt
      ! nx, ny, nt should be determined by the psm_to_tdsm_size() subroutine
      
      type(t_psm), intent(inout)      :: psm
      type(t_tdsm), intent(out) :: out
      real, intent(in)                :: maxdt
      integer, intent(in)             :: nt
      
      ! make some aliases to the parameters, so that the code get's readable
      real :: time0, depth0 !, mnull
      real :: north0, east0
      real    :: tfactor, dur_exc, prd, timestepsize
      integer :: it
      real :: rel_time
            
      time0  = psm%params(1)
      north0 = psm%params(2)
      east0  = psm%params(3)
      depth0 = psm%params(4)
  !    mnull  = psm%params(5)
      
      
      
      dur_exc      = psm%params(12)
      prd          = psm%params(13)
      timestepsize = maxdt
      
      ! build array of the moment-tensor-components in the size (6 * exdur/maxdt)
      ! start-values (first column) are the given, time-independent ones        
      ! fill the array with values of the stf(t):
       
      if (allocated( out%centroids )) deallocate(out%centroids)
      allocate( out%centroids(nt) )
        
      
      ! final moment has to be normalised with mnull after calculation of total source-signal 
      
      ! fill output array:
       
      do it=1,nt
         rel_time=(it-1) * timestepsize
         tfactor = stf( rel_time, prd, dur_exc ) 
         
         out%centroids(it)%north = north0
         out%centroids(it)%east  = east0
         out%centroids(it)%depth = depth0
         out%centroids(it)%time  = time0 + it*timestepsize
         out%centroids(it)%m(1)   = psm%params(6)  * tfactor
         out%centroids(it)%m(2)   = psm%params(7)  * tfactor
         out%centroids(it)%m(3)   = psm%params(8)  * tfactor
         out%centroids(it)%m(4)   = psm%params(9)  * tfactor
         out%centroids(it)%m(5)   = psm%params(10) * tfactor
         out%centroids(it)%m(6)   = psm%params(11) * tfactor
      end do
      
    end subroutine psm_to_tdsm_table_point_lp
    
    !---------------------------------------
    
    subroutine psm_write_info_file_point_lp( psm, fn )
      
      type(t_psm), intent(in)  :: psm
      type(varying_string), intent(in) :: fn
      
      integer :: unit
      real    :: north, east, depth

      north   = psm%params(2)
      east    = psm%params(3)
      depth   = psm%params(4)
      
      call claim_unit( unit )
      open( unit=unit, file=char(fn), status='unknown' )
      
      write (unit,"(a)") "origin"
      write (unit,*) psm%origin%lat, psm%origin%lon
      write (unit,*)
      
      write (unit,"(a)") "center"
      write (unit,*) north, east, depth
      write (unit,*)
      
!       write (unit,"(a)") "outline"
!       do i=1,36
!          write (unit,*) matmul(psm%rotmat_rup,(/radius*cos(-i/36.*2.*pi),radius*sin(-i/36.*2.*pi),0./)) &
!               +(/north,east,depth/)
!       end do
      
!       write (unit,*)
!       write (unit,"(a)") "rupture"
!       write (unit,*) matmul(psm%rotmat_rup,(/0.,0.,0./))+(/north,east,depth/)
!       write (unit,*) matmul(psm%rotmat_rup,(/radius,0.,0./))
!       write (unit,*)
      
!       write (unit,"(a)") "slip"
!       ! at position
!       write (unit,*) matmul(psm%rotmat_slip,(/0.,0.,-200./))+(/north, east, depth/)
!       ! slip vector
!       write (unit,*) matmul(psm%rotmat_slip,(/1000.,0.,0./))
!       ! at position
!       write (unit,*) matmul(psm%rotmat_slip,(/0.,0.,200./))+(/north, east, depth/)
!       ! slip vector
!       write (unit,*) matmul(psm%rotmat_slip,(/-1000.,0.,0./))
!       write (unit,*) 
      
!       write (unit,"(a)") "p-axis"
!       write (unit,*) psm%pax(:)
!       write (unit,*)
      
!       write (unit,"(a)") "t-axis"
!       write (unit,*) psm%tax(:)
!       write (unit,*)
      
        close( unit ) 
        call release_unit( unit )
      
    end subroutine psm_write_info_file_point_lp

    !---------------------------------------

    ! source time function is given explicitly via function "stf()"
    
    !stf(t) is suited for lp-event with main period = "prd" seconds
    
    !build in stf(t) = e^(-(t-4)^2/4)*1/(1+e^(-2t+4))*sin(2*PI/prd*t) for finite lp-excitation
    
    elemental function stf( reltime, prd, dur_exc )
      
      real, intent(in)  :: reltime, prd, dur_exc
      real :: stf, t1,t2,t3
      t1=2.
      t2=t1 + dur_exc-5
      t3=t2/4.
      
      stf = exp(-(reltime-t3)**2./(2*pi*dur_exc))*1./ &
           (1.+exp(-2.*(reltime-t1)))*1./(1.+exp(0.5*(reltime-t2)))*sin(2.*pi/prd*reltime)
      
    end function stf
    !---------------------------------------
    



end module
