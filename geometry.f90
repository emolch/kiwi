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

module geometry

    use constants
    
    implicit none
    
    private
    
    type, public :: t_halfspace
        real, dimension(3) :: point
        real, dimension(3) :: normal
    end type
    
    type, public :: t_circle
      ! Untransformed circle has unit radius and lies in x-y plane.
        real, dimension(3)   :: center
        real, dimension(3,3) :: transform
    end type
    
    type, public :: t_polygon
        real, dimension(:,:), allocatable :: points
    end type
    
    public point_in_halfspace
    public get_piercingpoint
    public nearest_point_on_polygon
    public circle_to_polygon
    public trim_polygon
    public copy_polygon
    public polygon_box
    public polygon_area
    public polygon_destroy
    
    interface trim_polygon
        module procedure trim_polygon_one
        module procedure trim_polygon_more
    end interface
    
  contains

    pure function point_in_halfspace( point, halfspace )
    
      ! Decide if point is in the halfspace.
      !
      ! Returns true if the point is in the halfspace, from which the normal 
      ! vector is pointing away, or if it lies on the borderplane.
      
        real, dimension(3), intent(in) :: point
        type(t_halfspace), intent(in)  :: halfspace
        logical :: point_in_halfspace
        
        point_in_halfspace = &
           (dot_product( halfspace%normal, halfspace%point(:)-point(:) ) >= 0.0)
            
    end function
    
    pure subroutine get_piercingpoint( point_a, point_b, halfspace, piercingpoint, &
                                  between_ab, parallel, a_inside_, b_inside_ )
    
      ! Get the piercing point of line between a and b with halfspace border.
        
        real, dimension(3), intent(in) :: point_a
        real, dimension(3), intent(in) :: point_b
        type(t_halfspace), intent(in)  :: halfspace
        real, dimension(3), intent(out) :: piercingpoint
        logical, intent(out)            :: between_ab, parallel
        logical, intent(out), optional  ::  a_inside_, b_inside_
        logical ::  a_inside, b_inside
        
        real :: lambda_a, lambda_b, lambda_ab
        real, dimension(3) :: ab
        
        ab(:) = point_b(:) - point_a(:)
        
        lambda_a = dot_product( halfspace%normal, halfspace%point(:) - point_a(:) )
        lambda_b = dot_product( halfspace%normal, halfspace%point(:) - point_b(:) )
        lambda_ab = dot_product( halfspace%normal, ab )
        
        a_inside = (lambda_a >= 0.)
        b_inside = (lambda_b >= 0.)
        if (present(a_inside_)) a_inside_ = a_inside
        if (present(b_inside_)) b_inside_ = b_inside
        between_ab = ((a_inside .and. .not. b_inside) .or. &
                      (b_inside .and. .not. a_inside))
       
        parallel = ( lambda_ab*lambda_ab < dot_product(ab,ab) / 2**(digits(lambda_ab)) )
        
        if (parallel .and. between_ab) then
            ! AB lies (almost) in plane
            if (abs(lambda_a) <= abs(lambda_b)) then
                piercingpoint(:) = point_a(:)
            else 
                piercingpoint(:) = point_b(:)
            end if
            return
        end if
        
        if (parallel .and. .not. between_ab) then
            piercingpoint = (/0.,0.,0./)
            return
        end if
        
        piercingpoint(:) = point_a(:) + ab(:)*lambda_a/lambda_ab
            
    end subroutine
    
    pure function nearest_point_on_polygon( polygon, point ) result(nearest_point)
         
        type(t_polygon), intent(in)    :: polygon
        real, dimension(3), intent(in) :: point
        real, dimension(3)             :: nearest_point
        
        type(t_halfspace)   :: halfspace
        real, dimension(3)  :: piercingpoint
        logical             :: does_pierce, parallel 
        integer             :: npoints, ipoint, jpoint
        real                :: distsquare, comp_distsquare
        
        nearest_point = point
        if (.not. allocated(polygon%points)) return
        if (size(polygon%points,2) == 0) return
        
        npoints = size(polygon%points,2)
        halfspace%point = point    
        distsquare = huge(distsquare)
        nearest_point = polygon%points(:,1)
        if (size(polygon%points,2) == 1) return
        
        do ipoint=1,npoints
            jpoint = mod(ipoint,npoints)+1
            halfspace%normal = polygon%points(:,jpoint) - polygon%points(:,ipoint)
            call get_piercingpoint( polygon%points(:,ipoint), polygon%points(:,jpoint), &
                                    halfspace, piercingpoint, &
                                    does_pierce, parallel )
            
            comp_distsquare = vecsqr(piercingpoint-point)
            if (does_pierce .and. comp_distsquare < distsquare) then
                distsquare = comp_distsquare
                nearest_point = piercingpoint
            end if
            
            comp_distsquare = vecsqr(polygon%points(:,ipoint)-point)
            if (comp_distsquare < distsquare) then
                distsquare = comp_distsquare
                nearest_point = polygon%points(:,ipoint)
            end if
        end do
    
    end function
    
    pure function vecsqr( a )
        real, dimension(3), intent(in)  :: a
        real :: vecsqr
        vecsqr = dot_product(a,a)
    end function
    
    pure subroutine circle_to_polygon( circle, npoints, polygon )
    
      ! Create a polygon approximating the circle.
    
        type(t_circle), intent(in)          :: circle
        integer, intent(in)                 :: npoints
        type(t_polygon), intent(inout)      :: polygon
        integer :: i
        
        if (allocated(polygon%points)) deallocate(polygon%points)
        allocate(polygon%points(3,npoints))
        
        do i=1,npoints
            polygon%points(:,i) = matmul(circle%transform,(/cos(i*2.*pi/npoints), &
                                                            sin(i*2.*pi/npoints), &
                                                            0./)) + circle%center(:)
        end do
        
    end subroutine
    
    pure subroutine trim_polygon_one( polygon, halfspace, polygon_trimmed )
    
      ! Cut off parts of polygon, which lie outside of halfspace.
      
        type(t_polygon), intent(in)              :: polygon
        type(t_halfspace), intent(in)            :: halfspace
        type(t_polygon), intent(inout)           :: polygon_trimmed
        
        real, dimension(3,size(polygon%points,2))  :: piercingpoints
        logical, dimension(size(polygon%points,2)) :: does_pierce, point_inside
        integer :: npoints, npoints_trimmed, ipoint, jpoint
        logical :: parallel, b_inside
        
        npoints = size(polygon%points,2)
        
        if (allocated(polygon_trimmed%points)) deallocate(polygon_trimmed%points)
        
        npoints_trimmed = 0
        do ipoint=1,npoints
            jpoint = mod(ipoint,npoints)+1
            call get_piercingpoint( polygon%points(:,ipoint), polygon%points(:,jpoint), &
                                    halfspace, piercingpoints(:,ipoint), &
                                    does_pierce(ipoint), parallel, & 
                                    point_inside(ipoint), b_inside )
        
            if (point_inside(ipoint)) npoints_trimmed = npoints_trimmed + 1
            if (does_pierce(ipoint)) npoints_trimmed = npoints_trimmed + 1
        end do
        
        allocate(polygon_trimmed%points(3,npoints_trimmed))
        
        jpoint = 1
        do ipoint=1,npoints
            if (point_inside(ipoint)) then
                polygon_trimmed%points(:,jpoint) = polygon%points(:,ipoint)
                jpoint = jpoint + 1
            end if
            if (does_pierce(ipoint)) then
                polygon_trimmed%points(:,jpoint) = piercingpoints(:,ipoint)
                jpoint = jpoint + 1
            end if
        end do 
        
    end subroutine
    
    pure subroutine trim_polygon_more( polygon, halfspaces, polygon_trimmed )
    
      ! Cut off parts of polygon, which lie outside of halfspaces.
      
        type(t_polygon), intent(in)                 :: polygon
        type(t_halfspace), dimension(:), intent(in) :: halfspaces
        type(t_polygon), intent(inout)              :: polygon_trimmed
        
        type(t_polygon) :: temp
        integer :: icon
        
        call copy_polygon( polygon, temp )
        do icon=1,size(halfspaces)
            if (icon /= 1) call copy_polygon( polygon_trimmed, temp )
            call trim_polygon( temp, halfspaces(icon), polygon_trimmed )
        end do

    end subroutine
    
    pure subroutine polygon_box( polygon, mi, ma )
    
         type(t_polygon), intent(in) :: polygon
         real, intent(out), dimension(3) :: mi, ma
         
         mi = minval(polygon%points,2)
         ma = maxval(polygon%points,2)
    
    end subroutine
    
    pure subroutine copy_polygon( a, b )
    
        type(t_polygon), intent(in)    :: a
        type(t_polygon), intent(inout) :: b
        
        if (allocated(b%points)) deallocate(b%points)
        allocate(b%points(3,size(a%points,2)))
        b%points = a%points
        
    end subroutine
    
    pure function polygon_area( polygon ) result(area)
    
        type(t_polygon), intent(in) :: polygon
        real :: area
        
        real :: area_xy, area_yz, area_zx
        integer :: ip, jp, np
        
        
        np = size(polygon%points,2)
        
        area_xy = 0.
        area_yz = 0.
        area_zx = 0.
        
        if (.not. allocated(polygon%points)) return
        
        np = size(polygon%points,2)
        if (np <= 2) return
        
        do ip=1,np
            jp = mod(ip,np)+1
            area_xy = area_xy + (polygon%points(1,ip)-polygon%points(1,jp))* &
                                (polygon%points(2,ip)+polygon%points(2,jp))*0.5
            area_yz = area_yz + (polygon%points(2,ip)-polygon%points(2,jp))* &
                                (polygon%points(3,ip)+polygon%points(3,jp))*0.5
            area_zx = area_zx + (polygon%points(3,ip)-polygon%points(3,jp))* &
                                (polygon%points(1,ip)+polygon%points(1,jp))*0.5
        end do
        
        area = sqrt(area_xy**2 + area_yz**2 + area_zx**2)
        
    end function
    
    
    
    pure subroutine polygon_destroy( poly )
    
        type(t_polygon), intent(inout) :: poly
        if (allocated(poly%points)) deallocate( poly%points )
    
    end subroutine
     
end module

