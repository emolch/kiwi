! $Id: source_all.f90 684 2007-12-03 11:31:59Z sebastian $ 
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

module source_all

! abstraction layer: interface to different source types
! fortran cannot do polymorphy so this ugly thing is needed

! constants, type definition and generic psm subroutines are in parameterized_source.f90
! implementations for different sources are in source_bilat, source_circular, ...
! the mimic of dynamic dispatching is done here
    
    use orthodrome
    use better_varying_string
    use unit
    use util
    use parameterized_source
    use discrete_source
    
! implementations: 

    use source_circular
    use source_bilat
    use source_point_lp
    use source_eikonal
    use source_moment_tensor
    
    implicit none
    
    private

    public psm_set
    public psm_to_tdsm
    public psm_write_info_file
    public psm_set_subparams
    public psm_get_n_source_params
    public psm_get_source_name
    public psm_get_source_id
    public psm_get_param_name
    public psm_get_param_unit
    public psm_get_param_id
    public psm_cleanup
    public psm_get_param_limits
    public psm_get_param_defaults
  
    integer, public , parameter                            :: nsourcetypes   = 5
    integer, public , parameter, dimension(nsourcetypes)   :: sourcetypes = &
     (/ psm_bilat, psm_circular, psm_point_lp, psm_eikonal, psm_moment_tensor /)
    
  contains
  
    subroutine psm_cleanup()
        call psm_cleanup_bilat()
        call psm_cleanup_circular()
        call psm_cleanup_point_lp()
        call psm_cleanup_eikonal()
        call psm_cleanup_moment_tensor()
    end subroutine
  
    
    subroutine psm_get_source_name( sourcetype, name )
        integer, intent(in)                  :: sourcetype
        type(varying_string), intent(inout)  :: name
        name = ""
        if (sourcetype .eq. psm_bilat)          name = "bilateral"
        if (sourcetype .eq. psm_circular)       name = "circular"
        if (sourcetype .eq. psm_point_lp)       name = "point_lp"
        if (sourcetype .eq. psm_eikonal)        name = "eikonal"
        if (sourcetype .eq. psm_moment_tensor)  name = "moment_tensor"
    end subroutine
    
    subroutine psm_get_source_id( name, sourcetype )
        type(varying_string), intent(in)      :: name
        integer, intent(out)                  :: sourcetype
        sourcetype = 0
        if (name .eq. "bilateral")      sourcetype = psm_bilat
        if (name .eq. "circular")       sourcetype = psm_circular
        if (name .eq. "point_lp")       sourcetype = psm_point_lp
        if (name .eq. "eikonal")        sourcetype = psm_eikonal
        if (name .eq. "moment_tensor")  sourcetype = psm_moment_tensor
    end subroutine
        
    integer function psm_get_n_source_params(sourcetype)
        
        integer, intent(in) :: sourcetype
    
        select case(sourcetype)
            
            case (psm_bilat) 
                psm_get_n_source_params = n_source_params_bilat
            
            case (psm_circular) 
                psm_get_n_source_params = n_source_params_circular
                
            case (psm_point_lp) 
                psm_get_n_source_params = n_source_params_point_lp
            
            case (psm_eikonal) 
                psm_get_n_source_params = n_source_params_eikonal

            case (psm_moment_tensor) 
                psm_get_n_source_params = n_source_params_moment_tensor
                
        end select
        
    end function
    
    subroutine psm_get_param_name( sourcetype, iparam, name )
        
        integer, intent(in)                     :: sourcetype, iparam
        type(varying_string), intent(inout)     :: name
    
        select case(sourcetype)
            
            case (psm_bilat) 
                call psm_get_param_name_bilat( iparam, name )
                
            case (psm_circular) 
                call psm_get_param_name_circular( iparam, name )
        
            case (psm_point_lp) 
                call psm_get_param_name_point_lp( iparam, name )
        
            case (psm_eikonal) 
                call psm_get_param_name_eikonal( iparam, name )

            case (psm_moment_tensor) 
                call psm_get_param_name_moment_tensor( iparam, name )
        
        end select

    end subroutine
    
    subroutine psm_get_param_unit( sourcetype, iparam, unit )
        
        integer, intent(in)                     :: sourcetype, iparam
        type(varying_string), intent(inout)     :: unit
    
        select case(sourcetype)
            
            case (psm_bilat) 
                call psm_get_param_unit_bilat( iparam, unit )
                
            case (psm_circular) 
                call psm_get_param_unit_circular( iparam, unit )
        
            case (psm_point_lp) 
                call psm_get_param_unit_point_lp( iparam, unit )
            
            case (psm_eikonal) 
                call psm_get_param_unit_eikonal( iparam, unit )

            case (psm_moment_tensor) 
                call psm_get_param_unit_moment_tensor( iparam, unit )
        
        end select

    end subroutine
    
    subroutine psm_get_param_id( sourcetype, name, iparam )
        
        integer, intent(in)                 :: sourcetype
        type(varying_string), intent(in)    :: name
        integer, intent(out)                :: iparam
    
        select case(sourcetype)
        
            case (psm_bilat) 
                call psm_get_param_id_bilat( name, iparam )
            
            case (psm_circular) 
                call psm_get_param_id_circular( name, iparam )
                
            case (psm_point_lp) 
                call psm_get_param_id_point_lp( name, iparam )
            
            case (psm_eikonal) 
                call psm_get_param_id_eikonal( name, iparam )

            case (psm_moment_tensor) 
                call psm_get_param_id_moment_tensor( name, iparam )
                
        end select

    end subroutine
        
    subroutine psm_set( psm, sourcetype, params, normalized_ )
    
        type(t_psm), intent(inout)     :: psm
        integer, intent(in)            :: sourcetype
        real, dimension(:), intent(in) :: params
        logical, intent(in), optional  :: normalized_
        
        logical :: normalized
        
        normalized = .false.
        if (present(normalized_)) normalized = normalized_
        
        select case(sourcetype)
            
            case (psm_bilat)
                call psm_set_bilat( psm, params, normalized )
            case (psm_circular)
                call psm_set_circular( psm, params, normalized )
            case (psm_point_lp)
                call psm_set_point_lp( psm, params, normalized )
            case (psm_eikonal)
                call psm_set_eikonal( psm, params, normalized )
            case (psm_moment_tensor)
                call psm_set_moment_tensor( psm, params, normalized )
            
        end select
       
        call resize( psm%params_mask, 1, size(psm%params) )
        if (psm%sourcetype /= sourcetype) then
            psm%params_mask(:) = .true.
        end if
        
        psm%sourcetype = sourcetype
        
    end subroutine
    
    
    subroutine psm_get_param_limits( sourcetype, limits_min, limits_max, hard_ )
    
        integer, intent(in)                  :: sourcetype
        real, intent(inout), dimension(:), allocatable :: limits_min
        real, intent(inout), dimension(:), allocatable :: limits_max
        logical, intent(in), optional :: hard_
        
        integer :: nparams
        logical :: hard
        
        hard = .false.
        if (present(hard_)) hard = hard_
        
        nparams = psm_get_n_source_params(sourcetype)
        call resize( limits_min, 1, nparams )
        call resize( limits_max, 1, nparams )
        
        select case(sourcetype)
            
            case (psm_bilat)
                if (hard) then
                    limits_min(:) = psm_params_min_hard_bilat(:)
                    limits_max(:) = psm_params_max_hard_bilat(:)
                else
                    limits_min(:) = psm_params_min_soft_bilat(:)
                    limits_max(:) = psm_params_max_soft_bilat(:)
                end if
                
            case (psm_circular)
                if (hard) then
                    limits_min(:) = psm_params_min_hard_circular(:)
                    limits_max(:) = psm_params_max_hard_circular(:)
                else
                    limits_min(:) = psm_params_min_soft_circular(:)
                    limits_max(:) = psm_params_max_soft_circular(:)
                end if
                
            case (psm_point_lp)
                if (hard) then
                    limits_min(:) = psm_params_min_hard_point_lp(:)
                    limits_max(:) = psm_params_max_hard_point_lp(:)
                else
                    limits_min(:) = psm_params_min_soft_point_lp(:)
                    limits_max(:) = psm_params_max_soft_point_lp(:)
                end if
            
            case (psm_eikonal)
                if (hard) then
                    limits_min(:) = psm_params_min_hard_eikonal(:)
                    limits_max(:) = psm_params_max_hard_eikonal(:)
                else
                    limits_min(:) = psm_params_min_soft_eikonal(:)
                    limits_max(:) = psm_params_max_soft_eikonal(:)
                end if

            case (psm_moment_tensor)
                if (hard) then
                    limits_min(:) = psm_params_min_hard_moment_tensor(:)
                    limits_max(:) = psm_params_max_hard_moment_tensor(:)
                else
                    limits_min(:) = psm_params_min_soft_moment_tensor(:)
                    limits_max(:) = psm_params_max_soft_moment_tensor(:)
                end if
                
        end select
        
    end subroutine
    
    subroutine psm_get_param_defaults( sourcetype, defaults )
    
        integer, intent(in)                  :: sourcetype
        real, intent(inout), dimension(:), allocatable :: defaults
        
        integer :: nparams
        
        nparams = psm_get_n_source_params(sourcetype)
        call resize( defaults, 1, nparams )
        
        select case(sourcetype)
            
            case (psm_bilat)
                defaults(:) = psm_params_default_bilat(:)
               
            case (psm_circular)
                defaults(:) = psm_params_default_circular(:)
            
            case (psm_point_lp)
                defaults(:) = psm_params_default_point_lp(:)
            
            case (psm_eikonal)
                defaults(:) = psm_params_default_eikonal(:)

            case (psm_moment_tensor)
                defaults(:) = psm_params_default_moment_tensor(:)
        
        end select
        
    end subroutine
    
    
    
    subroutine psm_set_subparams( psm, subparams, normalized_ )
    
      ! set params where params_mask is true with values from subparams array
      
        type(t_psm), intent(inout)        :: psm
        real, intent(in), dimension(:)    :: subparams
        logical, intent(in), optional     :: normalized_
    
        real, dimension(:), allocatable :: paramscopy
        integer :: iparam, isub
        logical :: normalized
        
        normalized = .false.
        if (present(normalized_)) normalized = normalized_
        
        if (size(psm%params) /= size(psm%params_mask) .or. &
            count(psm%params_mask) /= size(subparams)) then
            call die( "wrong sized arrays in psm_set_subparams()")
        end if
        
        call psm_get_params( psm, paramscopy, normalized )
        
        isub = 1
        do iparam=1, size(paramscopy)
            if (psm%params_mask(iparam)) then
                paramscopy(iparam) = subparams(isub)
                isub = isub+1
            end if
        end do
        
        select case(psm%sourcetype)
            
            case (psm_bilat)
                call psm_set_bilat( psm, paramscopy, normalized )
            case (psm_circular)
                call psm_set_circular( psm, paramscopy, normalized )
            case (psm_point_lp)
                call psm_set_point_lp( psm, paramscopy, normalized )
            case (psm_eikonal)
                call psm_set_eikonal( psm, paramscopy, normalized )
            case (psm_moment_tensor)
                call psm_set_moment_tensor( psm, paramscopy, normalized )
            
        end select

        call resize(paramscopy,1,0)
        
    end subroutine
    
    subroutine psm_to_tdsm( psm, tdsm, shortest_doi )
    
      ! translate a specific psm to tdsm
      ! shortest_duration is the shortest duration of interest
      ! automatically determines shape of output grid
      
        type(t_psm), intent(inout)    :: psm
        type(t_tdsm), intent(out)     :: tdsm
        real, intent(in)              :: shortest_doi        
        
        call tdsm_destroy( tdsm )

      ! here we simply hand over to the appropriate implementation
       
        select case(psm%sourcetype)
            
            case (psm_bilat)
                call psm_to_tdsm_bilat( psm, tdsm, shortest_doi )
            case (psm_circular)
                call psm_to_tdsm_circular( psm, tdsm, shortest_doi )
            case (psm_point_lp)
                call psm_to_tdsm_point_lp( psm, tdsm, shortest_doi )
            case (psm_eikonal)
                call psm_to_tdsm_eikonal( psm, tdsm, shortest_doi )
            case (psm_moment_tensor)
                call psm_to_tdsm_moment_tensor( psm, tdsm, shortest_doi )
                
        end select
        
        tdsm%origin = psm%origin
      
   end subroutine
   
   subroutine psm_write_info_file( psm, fn )
      
      ! write some information about the model psm into file fn
    
        type(t_psm), intent(in)  :: psm
        type(varying_string), intent(in) :: fn
        
      ! here we simply hand over to the appropriate implementation

        select case(psm%sourcetype)
            
            case (psm_bilat)
                call psm_write_info_file_bilat( psm, fn )
            case (psm_circular)
                call psm_write_info_file_circular( psm, fn )
            case (psm_point_lp)
                call psm_write_info_file_point_lp( psm, fn )
            case (psm_eikonal)
                call psm_write_info_file_eikonal( psm, fn )
            case (psm_moment_tensor)
                call psm_write_info_file_moment_tensor( psm, fn )
                
        end select
        
    end subroutine
    
end module
