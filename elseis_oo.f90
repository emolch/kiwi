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


module elseis_oo
    
    ! simple interface to elseis to wrap elementary seismogram in an object oriented way
    use util
    use elseis
    use differentiation
    use integration

    implicit none

    private 
    public set_coords, set_material, set_npq, set_stf, add_elseis_mt, add_elseis_sf, es_destroy
    
    type, public :: elseis_t
        integer                         :: n = 1
        integer                         :: p = 1
        integer                         :: q = 1
        real                            :: rho, alpha, beta
        real                            :: dt
        real, pointer, dimension(:)     :: stf
        real, allocatable, dimension(:) :: dstf, istf, istftau
        logical                         :: nf_flag, ff_flag
        
        real, dimension(3)              :: gamma
        real                            :: r = 1.
        
        real, dimension(5)              :: radiation_factor_mt, material_factor_mt, factor_mt
        real, dimension(5)              :: max_factor_mt
        
        real, dimension(3)              :: radiation_factor_sf, material_factor_sf, factor_sf
        real, dimension(3)              :: max_factor_sf
        
        
    end type elseis_t
    
  contains
    
    pure subroutine set_coords( es, coord, nfflag, ffflag )
    
        type(elseis_t), intent(inout)  :: es
        real, dimension(3), intent(in) :: coord
        logical, intent(in)            :: nfflag, ffflag
        
        real :: r
        real, dimension(3) :: gamma
        es % nf_flag = nfflag
        es % ff_flag = ffflag
        
        call make_direction_cosine( coord, gamma, r)
        
        if (r /= es % r) then
            es % r = r
            call update_factors( es )
        end if 
        
        if (any(gamma /= es % gamma)) then
            es % gamma(:) = gamma(:)
            call update_radiation_factors( es )
        end if 
        
    end subroutine
    
     pure subroutine set_material( es, rho, alpha, beta )
        
        type(elseis_t), intent(inout) :: es
        real, intent(in)              :: rho, alpha, beta
        
        es % rho = rho
        es % alpha = alpha
        es % beta = beta
        call update_material_factors( es )
        
    end subroutine set_material
    
    
    pure subroutine set_npq( es, n, p, q )
       
        type(elseis_t), intent(inout) :: es
        integer, intent(in)           :: n, p, q
        
        if (n /= es%n .or. p /= es%p .or. q /= es%q) then
            es%n = n
            es%p = p
            es%q = q
            call update_radiation_factors( es )
        end if
        
    end subroutine set_npq
    
    pure subroutine add_elseis_mt( es, toffset, weight, seismogram )
    
        type(elseis_t), intent(in)        :: es
        real, intent(in)                  :: toffset, weight
        real, dimension(:), intent(inout) :: seismogram
        
        call elseis_mt( es%factor_mt, es%r, es%alpha, es%beta, toffset, es%dt, &
                        es%stf, es%dstf, es%istf, es%istftau, &
                        es%nf_flag, es%ff_flag, seismogram, weight )
        
    end subroutine add_elseis_mt
    
    pure subroutine add_elseis_sf( es, toffset, weight, seismogram )
    
        type(elseis_t), intent(in)        :: es
        real, intent(in)                  :: toffset, weight
        real, dimension(:), intent(inout) :: seismogram
        
        call elseis_sf( es%factor_mt, es%r, es%alpha, es%beta, toffset, es%dt, &
                        es%stf, es%istf, es%istftau, &
                        es%nf_flag, es%ff_flag, seismogram, weight )
        
    end subroutine add_elseis_sf
    
    
    subroutine set_stf( es, stf, dt )
    
        type(elseis_t), intent(inout)           :: es
        real, dimension(:), intent(in), target  :: stf
        real, intent(in)                        :: dt
        
        integer :: lstf
        
        lstf = size(stf)
        es%dt = dt
        es%stf => stf
        
        if (allocated(es%dstf)) then
            if ( lstf /= size(es%dstf) ) then
                deallocate( es%istf )
                deallocate( es%istftau )
                deallocate( es%dstf )
            end if
        end if
        
        if ( .not. allocated( es%dstf ) ) then
            allocate( es%dstf(lstf) )
            allocate( es%istf(lstf) )
            allocate( es%istftau(lstf) )
        end if
        
        call make_istfs( dt, es%stf, es%istf, es%istftau )
        call differentiate( dt, es%stf, es%dstf )

    end subroutine set_stf

    subroutine es_destroy( es )
        type(elseis_t), intent(inout)           :: es
    
        if (allocated(es%dstf)) then
                deallocate( es%istf )
                deallocate( es%istftau )
                deallocate( es%dstf )
        end if
    
    end subroutine
            
    pure subroutine update_material_factors( es )
    
        type(elseis_t), intent(inout) :: es
     
        call material_factors_mt( es%rho, es%alpha, es%beta, es%material_factor_mt )
        call material_factors_sf( es%rho, es%alpha, es%beta, es%material_factor_sf )
        
        call update_factors( es )
        
    end subroutine update_material_factors
    
    pure subroutine update_radiation_factors( es )
        
        type(elseis_t), intent(inout) :: es
 
        call radpat_mt( es%gamma, es%n, es%p, es%q, es%radiation_factor_mt )
        call radpat_sf( es%gamma, es%n, es%p,       es%radiation_factor_sf )
        call update_factors( es )
        
    end subroutine update_radiation_factors
    
    
    pure subroutine update_factors( es )
    
        type(elseis_t), intent(inout) :: es
        
        call factors_mt( es%material_factor_mt, es%radiation_factor_mt, es%r, es%factor_mt )
        call factors_sf( es%material_factor_sf, es%radiation_factor_sf, es%r, es%factor_sf )

    end subroutine update_factors
    
    
end module elseis_oo
