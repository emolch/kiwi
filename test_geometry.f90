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

program test_geometry

    use util
    use constants
    use euler
    use orthodrome
    use geometry

    implicit none
    
    type(t_halfspace) :: halfspace
    real, dimension(3,4) :: points = reshape( (/0.,2.,-1.,  0.,-2.,-1.,  0.,0.,-1.,  0.,0.,0./), (/3,4/) )
    logical, dimension(4) :: expected = (/.true., .false., .true., .true./)
    real, dimension(3) :: piercingpoint
    logical :: between_ab, parallel
    integer :: i
    type(t_circle) :: circle
    real :: dip, strike
    type(t_polygon) :: polygon, polygon_trimmed
    real, dimension(3,7) :: expected_circ = reshape( (/0.14987442, 2.4953687, &
    2.6585152, -1.9344299, 0.9903534, 3.0681345, -2.5620692, -1.2604182, &
    1.9204066, -1.2604178, -2.5620692, 0.07959348, -1.1043297, -2.5185432, 0., &
    2.3468528, 0.9326396, 0., 2.12132, 2.1213207, 1.0000004 /), (/3,7/) )
    type(t_polygon) :: square
    
    
    
    call test_begin("test_geometry")
    
    halfspace%point = (/0.,0.,-1./)
    halfspace%normal = (/0.,-1.,-1./)
    
    do i=1,4
        if (.not. (expected(i) .eqv. point_in_halfspace( points(:,i), halfspace )) ) then
            call test_fail("point_in_halfspace")
        end if
    end do
    
    call get_piercingpoint( points(:,1), points(:,2), halfspace, piercingpoint, between_ab, parallel )
    if (.not. (all(piercingpoint == (/0.,0.,-1./) ) .and. &
              (between_ab .eqv. .true.) .and. (parallel .eqv. .false.))) then
        call test_fail("piercingpoint 1")
    end if
    
    
    call get_piercingpoint( (/0.,2.,-1./), (/0.,1.,-2./), halfspace, piercingpoint, between_ab, parallel )
    if (.not. (all(piercingpoint == (/0.,1.,-2./) ) .and. &
              (between_ab .eqv. .false.) .and. (parallel .eqv. .false.))) then
        call test_fail("piercingpoint 2")
    end if
    
    call get_piercingpoint( (/0.,2.,5./), (/0.,1.,1./), halfspace, piercingpoint, between_ab, parallel )
    if (.not. (all(near(piercingpoint,(/0.,0.4,-1.4/),0.0001)) .and. &
              (between_ab .eqv. .false.) .and. (parallel .eqv. .false.))) then
        call test_fail("piercingpoint 3")
    end if
    
    call get_piercingpoint( (/0.,1.,0./), (/0.,2.,-1.0001/), halfspace, piercingpoint, between_ab, parallel )
    if (.not. (all(piercingpoint == (/0.,0.,0./) ) .and. &
              (between_ab .eqv. .false.) .and. (parallel .eqv. .true.))) then
        call test_fail("piercingpoint 4")
    end if
    
    dip = d2r(45.)
    strike = d2r(45.)
    call init_euler(dip,strike,0.,circle%transform)
    circle%transform = circle%transform * 3.
    circle%center = (/0.,0.,1./)
    call circle_to_polygon( circle, 7, polygon )
    halfspace%point = (/0.,0.,0./)
    halfspace%normal = (/0.,0.,-1./)
    call trim_polygon( polygon, halfspace, polygon_trimmed )
    if (.not. all(near( polygon_trimmed%points,expected_circ, 0.00001))) then
        call test_fail("trimmed circle")
    end if
   
    
    dip = d2r(13.)
    strike = d2r(100.)
    call init_euler(dip,strike,0.,circle%transform)
    circle%transform = circle%transform * 1.
    circle%center = (/-5.,-2.,1./)
    call circle_to_polygon( circle, 3600, polygon )
    if (.not. near( polygon_area(polygon), pi, 0.0001)) then
        call test_fail("pi estimation")
    end if
    
    allocate( square%points(3,4) )
    square%points(:,:) = reshape( (/0.,0.,1.,  0.,2.,1.,  2.,2.,1.,  2.,0.,1./), (/3,4/) )
    if (.not. near( polygon_area(square), 4., 0.00001)) call test_fail("square area 1")
    
    
    square%points(:,:) = reshape( (/1.,0.,0.,  1.,0.,2.,  1.,2.,2.,  1.,2.,0./), (/3,4/) )
    if (.not. near( polygon_area(square), 4., 0.00001)) call test_fail("square area 2")
    
    square%points(:,:) = reshape( (/0.,1.,0.,  2.,1.,0.,  2.,1.,2.,  0.,1.,2./), (/3,4/) )
    if (.not. near( polygon_area(square), 4., 0.00001)) call test_fail("square area 3")
    
    call polygon_destroy( square )
    call polygon_destroy( polygon )
    call polygon_destroy( polygon_trimmed )
    call test_end()
    call cleanup()
    
    
contains
    
    elemental logical function near( a,b, eps )
        real, intent(in) :: a,b,eps
        near = abs(a-b) < eps
    end function
    
end program
