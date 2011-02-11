! $Id: seismogram.f90 687 2008-02-13 14:39:42Z sebastian $ 
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


module seismogram

    use constants
    use gfdb
    use discrete_source
    use sparse_trace
    use orthodrome
    use receiver
    use util
    use comparator
    
    implicit none
    
    private
    public :: make_seismogram
    
  contains
  
    subroutine make_seismogram( source, receiver, greensf, interpolate_, xundersample_, zundersample_ )
    
        type(t_tdsm),intent(in)                  :: source
        type(t_receiver),intent(inout)           :: receiver
        type(t_gfdb),intent(inout)               :: greensf
        logical, intent(in), optional            :: interpolate_
        integer, intent(in), optional            :: xundersample_, zundersample_
        
      ! Calculate seismogram 'displacement' at 'receiver' by superposing
      ! greens funtions 'greensf' for 'source'.
      ! The time range of the seismogram is always extended, to what is needed 
      ! to see the full seismogram.
        
        real*8 :: azi, bazi, dist
        real*8 :: azi_orig, bazi_orig, dist_orig, lambda
        real :: dnorth, deast
        real, dimension(6) :: m
        real, dimension(6) :: f
        integer :: icentroid, i
        integer, dimension(2) :: ix,iz
        real :: depth, time
        real :: cl, sl
        real :: rshift
        real, dimension(:), allocatable :: rshifts
        type(t_trace), pointer :: tracep
        type(t_strip), dimension(2) :: displacement_temp
        type(t_strip), dimension(2) :: displacement_ar
        logical :: need_horizontal
        logical :: interpolate
        integer :: xundersample, zundersample
        real :: dix, diz
        
      ! index of away, right, down, north, east component
        integer :: ja, jr, jd, jn, je
      ! sign of each component
        real    :: sa, sr, sd, sn, se
        tracep => null()
        
        interpolate = .false.
        if (present(interpolate_)) interpolate = interpolate_
        
        xundersample = 1
        if (present(xundersample_)) xundersample = xundersample_

        zundersample = 1
        if (present(zundersample_)) zundersample = zundersample_

      ! get indices of components or zero, if they are not set...
        ja = receiver_component_index( receiver, C_AWAY )
        jr = receiver_component_index( receiver, C_RIGHT )
        jd = receiver_component_index( receiver, C_DOWN )
        jn = receiver_component_index( receiver, C_NORTH )
        je = receiver_component_index( receiver, C_EAST )
        
        sa = receiver_component_sign( receiver, C_AWAY )
        sr = receiver_component_sign( receiver, C_RIGHT )
        sd = receiver_component_sign( receiver, C_DOWN )
        sn = receiver_component_sign( receiver, C_NORTH )
        se = receiver_component_sign( receiver, C_EAST )
        
        need_horizontal = ja .ne. 0 .or. jr .ne. 0 .or. jn .ne. 0 .or. je .ne. 0 
      
      ! azimuths and distances of individual subfaults will be
      ! taken relative to these:
        call azibazi(source%origin,receiver%origin, azi_orig, bazi_orig)
        dist_orig = distance_accurate50m(source%origin, receiver%origin)
        
        do i=1,receiver%ncomponents
            if (allocated(receiver%displacement(i)%data)) then
                receiver%displacement(i)%data(:)  = 0.
            end if
        end do
        
        
        if (need_horizontal) then
          ! if horizontal components are wanted, preextend temporary arrays to 
          ! reduce probability of resizes  (usually only the first will actually do sth...)
            if (ja .ne. 0) &
                call strip_extend_to_same_span_5(receiver%displacement(ja), displacement_ar(1), &
                    displacement_temp(1), displacement_ar(2), displacement_temp(2) )
            if (jr .ne. 0) &
                call strip_extend_to_same_span_5(receiver%displacement(jr), displacement_ar(1), &
                    displacement_temp(1), displacement_ar(2), displacement_temp(2) )
            if (jn .ne. 0) &
                call strip_extend_to_same_span_5(receiver%displacement(jn), displacement_ar(1), &
                    displacement_temp(1), displacement_ar(2), displacement_temp(2) )
            if (je .ne. 0) &
                call strip_extend_to_same_span_5(receiver%displacement(je), displacement_ar(1), &
                    displacement_temp(1), displacement_ar(2), displacement_temp(2) )
                    
            do i=1,2
                if (allocated(displacement_ar(i)%data)) then
                    displacement_ar(i)%data(:) = 0.
                end if
            end do
        end if
        do icentroid=1,size(source%centroids)
        
            dnorth = source%centroids(icentroid)%north
            deast = source%centroids(icentroid)%east
            depth = source%centroids(icentroid)%depth
            time = source%centroids(icentroid)%time
            m(:) = source%centroids(icentroid)%m
            
            rshift = time/greensf%dt
            call approx_differential_azidist( dnorth, deast, azi_orig, bazi_orig, dist_orig, &
                                              azi, bazi, dist )
            
            
            call make_weights( real(azi), m, f )
                        
            if (interpolate) then
                call gfdb_get_indices_bilin( greensf, real(dist), depth-receiver%depth, &
                                             xundersample, zundersample, ix,iz, dix, diz )
            else
                call gfdb_get_indices( greensf, real(dist), depth-receiver%depth, ix(1),iz(1) )
                ix(2) = ix(1)+1
                iz(2) = iz(1)+1
                dix = 0.
                diz = 0.
            end if
            
          ! horizontal components
            if ( need_horizontal) then
                lambda = bazi-bazi_orig
                if (lambda /= 0.) then  
                
                  ! take into account differing backazimuths for different centroids
                
                    cl = real(cos(lambda))
                    sl = real(sin(lambda))
                    
                  ! make seismogram of horizontal components in coordinates 
                  ! centered at at the current centroid
                    if (allocated(displacement_temp(1)%data) ) &
                        displacement_temp(1)%data(:) = 0.
                    call gfdb_get_trace_bilin( greensf,ix,iz,1, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_temp(1), f(1), rtraceshift_=rshift )
                    call gfdb_get_trace_bilin( greensf,ix,iz,2, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_temp(1), f(2), rtraceshift_=rshift )
                    call gfdb_get_trace_bilin( greensf,ix,iz,3, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_temp(1), f(3), rtraceshift_=rshift )
                    if (greensf%ng == 10) then
                        call gfdb_get_trace_bilin( greensf,ix,iz,9, dix, diz, tracep )
                        if (.not. associated(tracep)) cycle
                        call trace_multiply_add( tracep, displacement_temp(1), f(6), rtraceshift_=rshift )
                    end if

                    if (allocated(displacement_temp(2)%data) ) &
                        displacement_temp(2)%data(:) = 0.
                    call gfdb_get_trace_bilin( greensf,ix,iz,4, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_temp(2), f(4), rtraceshift_=rshift )
                    call gfdb_get_trace_bilin( greensf,ix,iz,5, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_temp(2), f(5), rtraceshift_=rshift )
                    
                  ! rotate horizontal compontents into coordinates centered at the master centroid
                    call strip_extend_to_same_span_4( displacement_temp(1), displacement_temp(2), &
                                                    displacement_ar(1), displacement_ar(2) )
        
                  ! add rotated seismogram of horizontal compontents to the master seismogram
                    displacement_ar(1)%data(:) = displacement_ar(1)%data(:) + &
                        cl*displacement_temp(1)%data(:) - sl*displacement_temp(2)%data(:)
                    displacement_ar(2)%data(:) = displacement_ar(2)%data(:) + &
                        cl*displacement_temp(2)%data(:) + sl*displacement_temp(1)%data(:)
                
                else 
                
                  ! not taking into account differing backazimuths for different centroids;
                  ! result can be accumulated directly in displacement_ar
                
                    call gfdb_get_trace_bilin( greensf,ix,iz,1, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_ar(1), f(1), rtraceshift_=rshift )
                    call gfdb_get_trace_bilin( greensf,ix,iz,2, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_ar(1), f(2), rtraceshift_=rshift )
                    call gfdb_get_trace_bilin( greensf,ix,iz,3, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_ar(1), f(3), rtraceshift_=rshift )
                    if (greensf%ng == 10) then
                        call gfdb_get_trace_bilin( greensf,ix,iz,9, dix, diz, tracep )
                        if (.not. associated(tracep)) cycle
                        call trace_multiply_add( tracep, displacement_ar(1), f(6), rtraceshift_=rshift )
                    end if
                    
                    call gfdb_get_trace_bilin( greensf,ix,iz,4, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_ar(2), f(4), rtraceshift_=rshift )
                    call gfdb_get_trace_bilin( greensf,ix,iz,5, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, displacement_ar(2), f(5), rtraceshift_=rshift )
                
                end if
            end if
            
          ! vertical component
            if (jd .ne. 0) then 
                call gfdb_get_trace_bilin( greensf,ix,iz,6, dix, diz, tracep )
                if (.not. associated(tracep)) cycle
                call trace_multiply_add( tracep, receiver%displacement(jd), f(1)*sd, rtraceshift_=rshift )
                call gfdb_get_trace_bilin( greensf,ix,iz,7, dix, diz, tracep )
                if (.not. associated(tracep)) cycle
                call trace_multiply_add( tracep, receiver%displacement(jd), f(2)*sd, rtraceshift_=rshift )
                call gfdb_get_trace_bilin( greensf,ix,iz,8, dix, diz, tracep )
                if (.not. associated(tracep)) cycle
                call trace_multiply_add( tracep, receiver%displacement(jd), f(3)*sd, rtraceshift_=rshift )
                
                if (greensf%ng == 10) then
                    call gfdb_get_trace_bilin( greensf,ix,iz,10, dix, diz, tracep )
                    if (.not. associated(tracep)) cycle
                    call trace_multiply_add( tracep, receiver%displacement(jd), f(6)*sd, rtraceshift_=rshift )
                end if

            end if
        end do

        if (need_horizontal) then
        
          ! put to final destination
            if (ja .ne. 0) then
                call strip_extend_to_same_span_2(receiver%displacement(ja), displacement_ar(1))
                receiver%displacement(ja)%data(:) = displacement_ar(1)%data(:) * sa
            end if
            if (jr .ne. 0) then
                call strip_extend_to_same_span_2(receiver%displacement(jr), displacement_ar(2))
                receiver%displacement(jr)%data(:) = displacement_ar(2)%data(:) * sr
            end if
    
          ! rotate (away,right) to (north,east) if needed
            if (jn .ne. 0 .or. je .ne. 0) then
                cl = real(cos(bazi_orig+pi))
                sl = real(sin(bazi_orig+pi))
                call strip_extend_to_same_span_2(displacement_ar(1),displacement_ar(2) )
                call rotate( displacement_ar(1)%data, displacement_ar(2)%data, cl, sl )
                if (jn .ne. 0) then
                    call strip_extend_to_same_span_2(receiver%displacement(jn), displacement_ar(1))
                    receiver%displacement(jn)%data(:) = displacement_ar(1)%data(:) * sn
                end if
                if (je .ne. 0) then
                    call strip_extend_to_same_span_2(receiver%displacement(je), displacement_ar(2))
                    receiver%displacement(je)%data(:) = displacement_ar(2)%data(:) * se
                end if
                
            end if
            
            call strip_destroy( displacement_temp(1) )
            call strip_destroy( displacement_temp(2) )
            call strip_destroy( displacement_ar(1) )
            call strip_destroy( displacement_ar(2) )
        end if
        do i=1,receiver%ncomponents
            if (any(isnan(receiver%displacement(i)%data(:)))) &
                call warn('NaN value in syntetic seismogram')
            if (any(abs(receiver%displacement(i)%data(:)) >= huge(depth))) &
                call warn('huge value in syntetic seismogram')
        end do
        
        
                      
        call gfdb_housekeeping( greensf )
        
    end subroutine
    
    elemental subroutine rotate( a, b, cos_angle, sin_angle)

        real, intent(inout) :: a, b
        real, intent(in) :: cos_angle, sin_angle

        real :: aa

        aa = cos_angle*a - sin_angle*b
        b = cos_angle*b + sin_angle*a
        a = aa

    end subroutine
   
    pure subroutine make_weights( azimuth, m, f )
    
        real, intent(in) :: azimuth
        real, intent(in), dimension(:) :: m    ! (6)
        real, intent(out), dimension(:) :: f   ! (6)
        
        real :: sa, ca, s2a, c2a
        
        sa = sin(azimuth)
        ca = cos(azimuth)
        s2a = sin(2.*azimuth)
        c2a = cos(2.*azimuth)
        
        f(1) = m(1)*ca**2 + m(2)*sa**2 + m(4)*s2a
        f(2) = m(5)*ca + m(6)*sa
        f(3) = m(3)
        f(4) = 0.5*(m(2)-m(1))*s2a + m(4)*c2a
        f(5) = m(6)*ca - m(5)*sa
        f(6) = m(1)*sa**2 + m(2)*ca**2 - m(4)*s2a   ! needed by near field terms
    
    end subroutine 
    
    pure subroutine intersection( a, b, c )
    
        integer, dimension(2), intent(in) :: a, b
        integer, dimension(2), intent(out) :: c
        
        c(1) = max(a(1),b(1))
        c(2) = min(a(2),b(2))
        
    end subroutine

end module
