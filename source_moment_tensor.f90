! $Id$ 
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

module source_moment_tensor

    use constants
    use orthodrome
    use parameterized_source
    use discrete_source
    use piecewise_linear_function
    use euler
    use better_varying_string
    use unit
    use util
    
    implicit none 

    private
    
    integer, public, parameter :: n_source_params_moment_tensor = 11
    integer, public, parameter :: n_grid_dims_moment_tensor = 1
    
    public psm_to_tdsm_moment_tensor, psm_write_info_file_moment_tensor, psm_set_moment_tensor
    public psm_get_param_name_moment_tensor, psm_get_param_unit_moment_tensor, psm_get_param_id_moment_tensor
    public psm_cleanup_moment_tensor
        
    type(varying_string), private, dimension(:), allocatable :: psm_param_names_moment_tensor, psm_param_units_moment_tensor
    
    real, dimension(n_source_params_moment_tensor), public, parameter :: psm_params_norm_moment_tensor = &
    (/ 1., 10000., 10000., 10000., 7e18, 7e18, 7e18, 7e18, 7e18, 7e18, 1. /)
    
    real, private, parameter :: big = huge(psm_params_norm_moment_tensor(1))
    
  ! logical, physical or computational accuracy limits for params
    real, dimension(n_source_params_moment_tensor), public, parameter :: psm_params_min_hard_moment_tensor = &
    (/ -big, -100000., -100000., 0., 0., 0., 0., 0., 0., 0., 0./)
    real, dimension(n_source_params_moment_tensor), public, parameter :: psm_params_max_hard_moment_tensor = &
    (/ big, 100000., 100000., 1000000., 7e25, 7e25, 7e25, 7e25, 7e25, 7e25, 10./)
    
  ! these limits give range of realistic, non-redundant or practical params
    real, dimension(n_source_params_moment_tensor), public, parameter :: psm_params_min_soft_moment_tensor = &
    (/ -20., -10000., -10000., 0., 0., 0., 0., 0., 0., 0., 0./)
    real, dimension(n_source_params_moment_tensor), public, parameter :: psm_params_max_soft_moment_tensor = &
    (/ 20., 10000., 10000., 150000.,  7e25, 7e25, 7e25, 7e25, 7e25, 7e25, 5./)
    
  ! defaults for external use
    real, dimension(n_source_params_moment_tensor), public, parameter :: psm_params_default_moment_tensor = &
    (/ 0.,0., 0., 10000., 0., 0., 0., 7e18, 0., 0., 1./)
    
  ! M0 = 7e18 Nm <=> Mw = 6.5
    
  contains
  
    subroutine psm_cleanup_moment_tensor()
    
        integer :: i
        
        if (allocated( psm_param_names_moment_tensor )) then
            do i=1,n_source_params_moment_tensor
                call delete( psm_param_names_moment_tensor(i) )
            end do
            deallocate( psm_param_names_moment_tensor )
        end if
        if (allocated( psm_param_units_moment_tensor )) then
            do i=1,n_source_params_moment_tensor
                call delete( psm_param_units_moment_tensor(i) )
            end do
            deallocate( psm_param_units_moment_tensor )
        end if
        
    end subroutine
    
    subroutine psm_init_param_names_moment_tensor()
    
        if (.not. allocated( psm_param_names_moment_tensor )) then
            allocate( psm_param_names_moment_tensor( n_source_params_moment_tensor ) )
            psm_param_names_moment_tensor(1) = "time"
            psm_param_names_moment_tensor(2) = "north-shift"
            psm_param_names_moment_tensor(3) = "east-shift"
            psm_param_names_moment_tensor(4) = "depth"
            psm_param_names_moment_tensor(5) = "mxx"
            psm_param_names_moment_tensor(6) = "myy"
            psm_param_names_moment_tensor(7) = "mzz"
            psm_param_names_moment_tensor(8) = "mxy"
            psm_param_names_moment_tensor(9) = "mxz"
            psm_param_names_moment_tensor(10) = "myz"
            psm_param_names_moment_tensor(11) = "rise-time"
        end if
        
        if (.not. allocated( psm_param_units_moment_tensor )) then
            allocate( psm_param_units_moment_tensor( n_source_params_moment_tensor ) )
            psm_param_units_moment_tensor(1) = "s"
            psm_param_units_moment_tensor(2) = "m"
            psm_param_units_moment_tensor(3) = "m"
            psm_param_units_moment_tensor(4) = "m"
            psm_param_units_moment_tensor(5) = "Nm"
            psm_param_units_moment_tensor(6) = "Nm"
            psm_param_units_moment_tensor(7) = "Nm"
            psm_param_units_moment_tensor(8) = "Nm"
            psm_param_units_moment_tensor(9) = "Nm"
            psm_param_units_moment_tensor(10) = "Nm"
            psm_param_units_moment_tensor(11) = "s"
        end if
        
    end subroutine
  
    subroutine psm_get_param_name_moment_tensor( iparam, name )
    
        integer, intent(in)                 :: iparam
        type(varying_string), intent(inout) :: name
        
        call psm_init_param_names_moment_tensor()
        
        name = ""
        if (iparam < 1 .or. n_source_params_moment_tensor < iparam) return
        name = psm_param_names_moment_tensor(iparam)
        
    end subroutine
    
    subroutine psm_get_param_unit_moment_tensor( iparam, unit )
    
        integer, intent(in)                 :: iparam
        type(varying_string), intent(inout) :: unit
        
        call psm_init_param_names_moment_tensor()
        
        unit = ""
        if (iparam < 1 .or. n_source_params_moment_tensor < iparam) return
        unit = psm_param_units_moment_tensor(iparam)
        
    end subroutine
   
    subroutine psm_get_param_id_moment_tensor( name, iparam )
    
        type(varying_string), intent(in)     :: name
        integer, intent(out)                 :: iparam

        integer :: i
        
        call psm_init_param_names_moment_tensor()
        iparam = 0
        do i=1, n_source_params_moment_tensor
            if (name .eq. psm_param_names_moment_tensor(i)) then
                iparam = i
                return
            end if
        end do
    
    end subroutine
  
    subroutine psm_set_moment_tensor( psm,  params, normalized_ )
    
        type(t_psm), intent(inout)     :: psm
        real, dimension(:), intent(in) :: params
        logical, intent(in), optional  :: normalized_
        
        logical :: must_reset_grid
        logical :: normalized
        
        normalized = .false.
        if (present(normalized_)) normalized = normalized_
       
        must_reset_grid = .false.
        if (.not. allocated(psm%grid_size) .or. (psm%sourcetype .ne. psm_moment_tensor)) then
            must_reset_grid = .true.
        end if
        
        call resize( psm%params, 1, n_source_params_moment_tensor )
        call resize( psm%grid_size, 1, n_grid_dims_moment_tensor )
        call resize( psm%params_norm, 1, n_source_params_moment_tensor )
        
        psm%params_norm(:) = psm_params_norm_moment_tensor(:)
        
        if (must_reset_grid) then
            psm%grid_size = 1
        end if
        
        if (size(params,1) .ne. size(psm%params)) call die("wrong number of source parameters in psm_set_moment_tensor()")
                
        if (normalized) then
            psm%params = params * psm%params_norm
        else
            psm%params = params 
        end if
        
    end subroutine

    subroutine psm_to_tdsm_moment_tensor( psm, tdsm, shortest_doi, ok )
    
      ! translate a specific psm to tdsm
      ! shortest_duration is the shortest duration of interest
      ! automatically determines shape of output grid
      
        type(t_psm), intent(inout) :: psm
        type(t_tdsm), intent(inout) :: tdsm
        real, intent(in) :: shortest_doi
        logical, intent(out) :: ok
        
        real :: maxdt, risetime, tbeg, dt, ta, tb
        integer :: it, nt
        real, dimension(3) :: point 
        real, dimension(6) :: m
        type(t_plf) :: stf
        real, dimension(:), allocatable  :: wt, toff

        ok = .true.
        point(:) = psm%params(2:4)
        m(:) = psm%params(5:10)
        risetime = psm%params(11)
        
        maxdt = shortest_doi
        
      ! set temporal extension of source grid according to risetime
        nt = floor(risetime/maxdt) + 1
      ! need at least 2 points to get partial deriv. depend on risetime
        if (nt .le. 1) nt = 2
  
        psm%grid_size(1) = nt

      ! stf is a boxcar of length risetime
        call plf_make( stf, (-risetime)/2., 0., &
                            (-risetime)/2., 1./risetime, &
                            (risetime)/2.,  1./risetime, &
                            (risetime)/2.,  0. )
                
        tbeg = stf%f(1,1)
        dt = risetime/nt
        
        allocate(wt(nt),toff(nt))
      ! weight on each time interval
      ! and 
        do it=1,nt
            ta = tbeg+dt*(it-1)
            tb = tbeg+dt*it
            call plf_integrate_and_centroid( stf, ta, tb, wt(it), toff(it) )
        end do

        if (allocated( tdsm%centroids )) deallocate(tdsm%centroids)
        allocate( tdsm%centroids(nt) )
        do it=1,nt
            tdsm%centroids(it)%north = point(1)
            tdsm%centroids(it)%east = point(2)
            tdsm%centroids(it)%depth = point(3)
            tdsm%centroids(it)%time = toff(it)
            tdsm%centroids(it)%m(:) = m(:)*wt(it)
        end do
        deallocate(wt,toff)
        
    end subroutine
    
    subroutine psm_write_info_file_moment_tensor( psm, fn )
        
        type(t_psm), intent(in)  :: psm
        type(varying_string), intent(in) :: fn

        integer :: unit
        real    :: north, east, depth

        north = psm%params(2)
        east = psm%params(3)
        depth = psm%params(4)
       
        call claim_unit( unit )
        open( unit=unit, file=char(fn), status='unknown' )
        
        write (unit,"(a)") "origin"
        write (unit,*) psm%origin%lat, psm%origin%lon
        write (unit,*)

        write (unit,"(a)") "center"
        write (unit,*) north, east, depth
        write (unit,*)
        
        
        close( unit ) 
        call release_unit( unit )
    
    end subroutine

end module
