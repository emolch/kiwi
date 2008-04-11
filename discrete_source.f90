! $Id: discrete_source.f90 687 2008-02-13 14:39:42Z sebastian $ 
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


module discrete_source

    use unit
    use orthodrome
    use better_varying_string
    
    implicit none
    
  ! this module defines types for the discretized source
  
    type, public :: t_centroid
        real :: north, east, depth, time
        real, dimension(6) :: m                      ! (6) mxx myy mzz mxy mxz myz
    end type
    
    type, public :: t_tdsm 
      
      ! time-domain discrete source model
    
        double precision   :: ref_time           ! reference time
        type(t_geo_coords) :: origin             ! reference origin at surface
        
      ! table of moment tensor approximating the rupture
        type(t_centroid), allocatable, dimension(:) :: centroids
      
      ! mean centroid
        type(t_centroid) :: centroid

      ! a constant stf may be used to speedup seismogram calculation
        real, dimension(:), allocatable :: const_stf_shifts, const_stf_amplitudes
        
    end type

    private
    public tdsm_write_info_file, tdsm_destroy
    
  contains
    
    subroutine tdsm_write_info_file( dsm, fn )
        
        type(t_tdsm), intent(in)  :: dsm
        type(varying_string), intent(in) :: fn

        integer :: unit, ncentroids
       
        ncentroids = 0
        if (allocated(dsm%centroids)) then
            ncentroids = size(dsm%centroids)
        end if
        
        call claim_unit( unit )
        open( unit=unit, file=char(fn), status='unknown' )
        
        write (unit,"(a)") "ncentroids"
        write (unit,*) ncentroids
        write (unit,*)
        
        close( unit ) 
        call release_unit( unit )
    
    end subroutine

    subroutine tdsm_destroy( dsm )
    
        type(t_tdsm), intent(inout)  :: dsm
        
        if (allocated(dsm%centroids)) deallocate( dsm%centroids )
        if (allocated(dsm%const_stf_shifts)) deallocate(dsm%const_stf_shifts)
        if (allocated(dsm%const_stf_amplitudes)) deallocate(dsm%const_stf_amplitudes)
    end subroutine
    
end module
