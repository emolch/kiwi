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

module parameterized_source

! see also source.f90

    use util
    use orthodrome
    use geometry
    use crust2x2
    
    implicit none
    
    private
    public eikonal_grid_destroy
    
    public psm_destroy
    public psm_reset_dependents
    public psm_set_default_constraints
    public psm_set_constraints    
    public psm_set_crustal_thickness_limit
    public psm_get_crustal_thickness
    public psm_point_in_constraints
    public psm_set_origin_and_time
    public psm_get_params
    public psm_set_mask
    public psm_get_subparams
    public psm_get_subparams_norm
  
  ! source types
    integer, public, parameter :: psm_bilat = 1
    integer, public, parameter :: psm_circular = 2
    integer, public, parameter :: psm_point_lp = 3
    integer, public, parameter :: psm_eikonal = 4
    integer, public, parameter :: psm_mt_eikonal = 5
    integer, public, parameter :: psm_moment_tensor = 6
    
    type, public :: t_eikonal_grid
        real, dimension(2) :: first
        real, dimension(2) :: last
        real, dimension(2) :: delta
        real, dimension(2) :: initialpoint
        integer, dimension(2) :: ndims
        real :: minspeed
        real, dimension(:,:), allocatable :: speed, times, durations, weights
        real, dimension(:,:,:), allocatable :: points
    end type
    
    type, public :: t_psm
        
        integer                           :: sourcetype = 0
    
        double precision                  :: ref_time           ! reference time
        type(t_geo_coords)                :: origin             ! reference origin at surface
        
        real                              :: moment = 1.0       ! may be used to apply moment after seismogram generation
        real                              :: risetime = 0.0     ! may be used to apply (constant) risetime after seismogram generation
        real, dimension(:),allocatable    :: params             ! the parameters of the source
        
        
        logical, dimension(:),allocatable :: params_mask        ! mask, which can be used
                                                                ! to set which params may
                                                                ! vary
                                                                
        real, dimension(:), allocatable   :: params_norm
                                                                
      ! rotation matrices to get slip and rupture vectors
        real, dimension(3,3)              :: rotmat_rup, rotmat_slip
      ! p-axis and t-axis
        real, dimension(2)                :: pax, tax
        
      ! caches last used grid dimensions; for example: (nx,ny,nt) for bilat source
        integer, dimension(:),allocatable :: grid_size
        
        real                              :: crustal_thickness_limit
        type(t_halfspace), dimension(:), allocatable :: constraints  ! not all source models respect these
        
        type(t_eikonal_grid)                         :: egrid         ! not all source models set these
        type(t_eikonal_grid)                         :: cgrid 
    end type
    
  contains
  
    subroutine eikonal_grid_destroy( egrid )
        type(t_eikonal_grid), intent(inout) :: egrid
        if (allocated(egrid%speed)) deallocate(egrid%speed)
        if (allocated(egrid%times)) deallocate(egrid%times)
        if (allocated(egrid%durations)) deallocate(egrid%durations)
        if (allocated(egrid%weights)) deallocate(egrid%weights)
        if (allocated(egrid%points)) deallocate(egrid%points)
    end subroutine
  
    subroutine psm_destroy( psm )
        type(t_psm), intent(inout)      :: psm
        if (allocated(psm%params)) deallocate(psm%params)
        if (allocated(psm%params_mask)) deallocate(psm%params_mask)
        if (allocated(psm%grid_size)) deallocate(psm%grid_size)
        if (allocated(psm%params_norm)) deallocate(psm%params_norm)
        if (allocated(psm%constraints)) deallocate(psm%constraints)
        call eikonal_grid_destroy( psm%egrid )
        call eikonal_grid_destroy( psm%cgrid )
        psm%sourcetype = 0
        psm%moment = 1.0
        psm%risetime = 0.0
    end subroutine
  
    subroutine psm_reset_dependents( psm )
        type(t_psm), intent(inout)      :: psm
        psm%moment = 1.0
        psm%risetime = 0.0
    end subroutine
  
    subroutine psm_set_default_constraints( self )
        
        type(t_psm), intent(inout)  :: self
        real :: thickness
        
        if (allocated(self%constraints)) deallocate(self%constraints)
        allocate( self%constraints(2) )

        call psm_get_crustal_thickness( self, thickness )        

      ! surface constraint:
        self%constraints(1)%point = (/0.,0.,1500./)
        self%constraints(1)%normal = (/0.,0.,-1./)
        
      ! crust bottom constraint:
        self%constraints(2)%point = (/0.,0.,thickness/)
        self%constraints(2)%normal = (/0.,0.,1./)
        
    end subroutine

    subroutine psm_set_constraints( self, points, normals )

        type(t_psm), intent(inout)       :: self
        real, dimension(:,:), intent(in) :: points
        real, dimension(:,:), intent(in) :: normals
        
        integer :: i,n

        n = size(points,2)

        if (allocated(self%constraints)) deallocate(self%constraints)
        allocate( self%constraints(n) )

        do i=1,n
            self%constraints(i)%point(:) = points(:,i)
            self%constraints(i)%normal(:) = normals(:,i)
        end do
        
        

    end subroutine


    pure function psm_point_in_constraints( self, point )
        type(t_psm), intent(in)  :: self
        real, dimension(3), intent(in) :: point
        logical :: psm_point_in_constraints
        integer :: icon 
        
        psm_point_in_constraints = .true.
        do icon=1,size(self%constraints)
            if (.not. point_in_halfspace( point, self%constraints(icon) )) then
                psm_point_in_constraints = .false.
                return
            end if
        end do
    end function
        
    subroutine psm_set_origin_and_time(psm, origin, ref_time)
    
        type(t_psm), intent(inout)      :: psm
        type(t_geo_coords), intent(in)  :: origin
        double precision, intent(in)    :: ref_time
        
        psm%origin = origin
        psm%ref_time = ref_time
        if ( crust2x2_loaded ) then
            call psm_set_default_constraints( psm )
        end if
    end subroutine
    
    subroutine psm_set_crustal_thickness_limit( self, thickness_limit )
        type(t_psm), intent(inout)      :: self
        real, intent(in)                :: thickness_limit
        self%crustal_thickness_limit = thickness_limit
        if ( crust2x2_loaded ) then
            call psm_set_default_constraints( self )
        end if
    end subroutine

    subroutine psm_get_crustal_thickness( self, thickness )

        type(t_psm), intent(in)      :: self
        real, intent(out)            :: thickness

        real                         :: vp, vs, vrho
        type(t_crust2x2_1d_profile)  :: profile
       
        call crust2x2_get_profile(r2d(self%origin), profile)
        call crust2x2_get_profile_averages(profile, vp, vs, vrho, thickness)
        if (self%crustal_thickness_limit > 0) then
            thickness = min(self%crustal_thickness_limit, thickness)
        end if
        
    end subroutine
    
    subroutine psm_get_params( psm, params, normalized_  )
    
        type(t_psm), intent(inout)     :: psm
        real, intent(inout), dimension(:), allocatable :: params
        logical, intent(in), optional     :: normalized_
        
        logical :: normalized

        normalized = .false.
        if (present(normalized_)) normalized = normalized_
        
        call resize( params, 1, size(psm%params) )
        
        if (normalized) then
            params(:) = psm%params(:) / psm%params_norm
        else
            params(:) = psm%params(:)
        end if
    
    end subroutine
    
    subroutine psm_set_mask( psm, params_mask )
    
        type(t_psm), intent(inout)          :: psm
        logical, dimension(:), intent(in)   :: params_mask

        if (size(params_mask) /= size(psm%params_mask)) then
            call die( "wrong length of source params mask" )
        end if
        
        psm%params_mask(:) = params_mask(:)
        
    end subroutine
    
    subroutine psm_get_subparams( psm, subparams, normalized_ )
    
      ! get a continuous array of all parameters, where params_mask is true
     
        type(t_psm), intent(in)                        :: psm
        real, intent(inout), dimension(:), allocatable :: subparams
        logical, intent(in), optional     :: normalized_
        
        integer :: iparam, isub
        logical :: normalized
        
        normalized = .false.
        if (present(normalized_)) normalized = normalized_
        
        if (size(psm%params) /= size(psm%params_mask)) then
            call die( "wrong sized arrays in psm_get_subparams()")
        end if
        
        call resize( subparams, 1, count(psm%params_mask) )
        
        isub = 1
        do iparam=1, size(psm%params)
            if (psm%params_mask(iparam)) then
                if (normalized) then
                    subparams(isub) = psm%params(iparam) / psm%params_norm(iparam)
                else
                    subparams(isub) = psm%params(iparam)
                end if
                isub = isub+1
            end if
        end do

    end subroutine
    
    subroutine psm_get_subparams_norm( psm, subparams_norm )

        type(t_psm), intent(in)                        :: psm
        real, intent(inout), dimension(:), allocatable :: subparams_norm
        
        integer :: iparam, isub
        
        if (size(psm%params) /= size(psm%params_mask)) then
            call die( "wrong sized arrays in psm_get_subparams()")
        end if
        
        call resize( subparams_norm, 1, count(psm%params_mask) )

        isub = 1
        do iparam=1, size(psm%params)
            if (psm%params_mask(iparam)) then
                subparams_norm(isub) = psm%params_norm(iparam)
                isub = isub+1
            end if
        end do

    end subroutine
    
    
    
end module
