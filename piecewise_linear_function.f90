! $Id: piecewise_linear_function.f90 703 2008-04-03 15:51:31Z sebastian $ 
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


module piecewise_linear_function

! this module may be used to deal with piecewise linear functions.
! this may be useful to make boxcars, ramps, triangles, trapezoids ...
! the function is assumed to jump to zero at the endpoints.

    use constants    

    implicit none

    type, public :: t_plf
    
      ! these coordinates define the function
        real, allocatable, dimension(:,:)  :: f
      
      ! number of points
        integer               :: n = 0
        
    end type
    
    public plf_make, plf_integrate, plf_destroy, plf_span
    private trapezoid_area
    public ip_cos, ip_linear, ip_zero_one
    
    interface plf_make
        module procedure plf_make_from_args
        module procedure plf_make_from_arrays
    end interface
    
    interface plf_taper_array
        module procedure plf_taper_array_r
        module procedure plf_taper_array_c
    end interface

  contains
    
    pure subroutine plf_make_from_args( s, x1,y1, x2,y2, x3,y3, x4,y4 )         
      ! make 2-4 point plf
      
        real,intent(in) :: x1,y1, x2,y2
        real,intent(in),optional :: x3,y3, x4,y4
        type(t_plf),intent(inout) :: s
        
        call plf_destroy( s )
        allocate( s%f(2,4) )
        
        s%f(1,1)=x1; s%f(1,2)=x2;
        s%f(2,1)=y1; s%f(2,2)=y2; 
        s%n = 2
        
        if (present(x3) .and. present(y3)) then
            s%f(1,3)=x3;
            s%f(2,3)=y3;
            s%n = 3
        end if
        
        if (present(x4) .and. present(y4)) then
            s%f(1,4)=x4
            s%f(2,4)=y4
            s%n = 4
        end if
        
    end subroutine
    
    pure subroutine plf_make_from_arrays( s, x, y )
    
        type(t_plf),intent(inout) :: s
        real, dimension(:), intent(in) :: x, y
        integer :: n
        
        call plf_destroy( s )
        
        n = min(size(x),size(y))
        allocate( s%f(2,n) )
        s%f(1,:) = x(1:n)
        s%f(2,:) = y(1:n)
        s%n = n
        
    end subroutine
    
    pure subroutine plf_destroy( s )
         type(t_plf),intent(inout) :: s
   
        if (allocated(s%f)) deallocate( s%f )
    
    end subroutine
    
    pure function plf_defined( s )
        type(t_plf),intent(in) :: s
        logical :: plf_defined
        plf_defined = allocated(s%f)
        
    end function
    
    pure subroutine plf_copy( s, d )
        type(t_plf),intent(in) :: s
        type(t_plf),intent(inout) :: d
   
        call plf_destroy( d )
        allocate( d%f(2,s%n) )
        d%f(:,:) = s%f(:,:)
        d%n = s%n
 
    end subroutine
    
    pure function plf_span( s )
        type(t_plf), intent(in) :: s
        real, dimension(2) :: plf_span
        
        if (allocated(s%f)) then
            plf_span(1) = s%f(1,1)
            plf_span(2) = s%f(1,s%n)
        else
            plf_span = (/0.,-1./)
        end if

    end function

    pure function plf_integrate( s, a, b ) result(area)
    
      ! get area between x=a and x=b
      
        type(t_plf), intent(in) :: s
        real, intent(in)        :: a,b
        real                    :: area
        real                    :: x0,x1,y0,y1
        integer                 :: i
        
        area = 0.
        if (.not. allocated( s%f )) return
        if (b .le. s%f(1,1)) return
        if (a .ge. s%f(1,s%n)) return
        do i=1,s%n-1
            if (a .ge. s%f(1,i+1)) cycle
            if (b .le. s%f(1,i)) return
            x0 = max( a, s%f(1,i) )
            x1 = min( b, s%f(1,i+1) )
            y0 = s%f(2,i)
            if (x0 .ne. s%f(1,i)) y0 = ip_linear( s%f(1,i), s%f(2,i), s%f(1,i+1), s%f(2,i+1), a)
            y1 = s%f(2,i+1)
            if (x1 .ne. s%f(1,i+1)) y1 = ip_linear( s%f(1,i), s%f(2,i), s%f(1,i+1), s%f(2,i+1), b)
            area = area+trapezoid_area( x0,y0, x1,y1 )
        end do
        
    end function
    
    pure subroutine plf_integrate_and_centroid( s, a, b, area, centroid )
    
      ! get area between x=a and x=b
      
        type(t_plf), intent(in) :: s
        real, intent(in)        :: a,b
        real,intent(out)        :: area, centroid
        real                    :: x0,x1,y0,y1, c, areathis
        integer                 :: i
        
        area = 0.
        centroid = (a+b)/2.
        c = 0.
        if (.not. allocated( s%f )) return
        if (b .le. s%f(1,1)) return
        if (a .ge. s%f(1,s%n)) return
        do i=1,s%n-1
            if (a .ge. s%f(1,i+1)) cycle
            if (b .le. s%f(1,i)) exit
            x0 = max( a, s%f(1,i) )
            x1 = min( b, s%f(1,i+1) )
            y0 = s%f(2,i)
            if (x0 .ne. s%f(1,i)) y0 = ip_linear( s%f(1,i), s%f(2,i), s%f(1,i+1), s%f(2,i+1), a)
            y1 = s%f(2,i+1)
            if (x1 .ne. s%f(1,i+1)) y1 = ip_linear( s%f(1,i), s%f(2,i), s%f(1,i+1), s%f(2,i+1), b)
            areathis = trapezoid_area( x0,y0, x1,y1 )
            c = c + areathis * trapezoid_centroid( x0,y0, x1,y1 )
            area = area + areathis
        end do
        centroid = c/area
    end subroutine
 
    subroutine plf_taper_array_r( s, array, span, dx, ip_method )
    
        type(t_plf), intent(in) :: s
        integer, dimension(2), intent(in) :: span
        real, dimension(span(1):), intent(inout) :: array
        real, intent(in) :: dx
        interface
            pure function ip_method( x0,y0, x1,y1, xi ) result(yi)
                real, intent(in) :: x0,y0,x1,y1,xi
                real :: yi
            end function
        end interface
        
        integer :: ibeg, iend, i, j, ibegatleast
        
        ibeg = floor(s%f(1,1)/dx)
              
      ! fill zeroes before start of plf
        if (span(1) <= ibeg) then
            array(span(1):min(ibeg,span(2))) = 0.
        end if
      
      ! interpolate linear ramps  
        ibegatleast = span(1)
        do i=1,s%n-1
            ibeg = max(floor(s%f(1,i)/dx)+1,span(1),ibegatleast)
            iend = min(floor(s%f(1,i+1)/dx),span(2))
            if (ibeg <= iend) then
                forall (j=ibeg:iend)
                    array(j) = array(j) * ip_method( s%f(1,i), s%f(2,i), &
                                               s%f(1,i+1), s%f(2,i+1), j*dx )
                end forall
            end if
            ibegatleast = iend+1
        end do
        
      ! fill zeroes until end
        iend = floor(s%f(1,s%n)/dx)+1 
        if (span(2) >= iend) then
            array(max(iend,span(1)):span(2)) = 0.
        end if
        
    end subroutine

    subroutine plf_taper_array_c( s, array, span, dx, ip_method )
    
        type(t_plf), intent(in) :: s
        integer, dimension(2), intent(in) :: span
        complex, dimension(span(1):), intent(inout) :: array
        real, intent(in) :: dx
        interface
            pure function ip_method( x0,y0, x1,y1, xi ) result(yi)
                real, intent(in) :: x0,y0,x1,y1,xi
                real :: yi
            end function
        end interface
        
        
        integer :: ibeg, iend, i, j, ibegatleast
        
        ibeg = floor(s%f(1,1)/dx)
              
      ! fill zeroes before start of plf
        if (span(1) <= ibeg) then
            array(span(1):min(ibeg,span(2))) = 0.
        end if
      
      ! interpolate linear ramps  
        ibegatleast = span(1)
        do i=1,s%n-1
            ibeg = max(floor(s%f(1,i)/dx)+1,span(1),ibegatleast)
            iend = min(floor(s%f(1,i+1)/dx),span(2))
            if (ibeg <= iend) then
                forall (j=ibeg:iend)
                    array(j) = array(j) * ip_method( s%f(1,i), s%f(2,i), &
                                               s%f(1,i+1), s%f(2,i+1), j*dx )
                end forall
            end if
            ibegatleast = iend+1
        end do
        
      ! fill zeroes until end
        iend = floor(s%f(1,s%n)/dx)+1 
        if (span(2) >= iend) then
            array(max(iend,span(1)):span(2)) = 0.
        end if
        
    end subroutine
   
    
    pure function trapezoid_centroid( x0,y0, x1,y1 ) result(x)
        real, intent(in) :: x0,y0,x1,y1
        real :: x
      ! returns center of mass of trapezoid
        if (y0+y1 == 0.) then
            x = (x0+x1)/2.
        else
            x = (x0*(2.*y0+y1)+x1*(y0+2.*y1))/(3.*(y0+y1))
        end if
    end function
    
    pure function trapezoid_area( x0,y0, x1,y1 ) result(a)
        real, intent(in) :: x0,y0,x1,y1
        real :: a
        a = (y0+y1)*(x1-x0)/2.
    end function
    
    pure function ip_linear( x0,y0, x1,y1, xi ) result(yi)
        real, intent(in) :: x0,y0,x1,y1,xi
        real :: yi
        yi = y0+(y1-y0)/(x1-x0)*(xi-x0)
    end function
    
    pure function ip_cos( x0,y0, x1,y1, xi ) result(yi)
        real, intent(in) :: x0,y0,x1,y1,xi
        real :: yi
        if (y1/=y0) then
            yi = y0+(y1-y0)*(0.5-0.5*cos( (xi-x0)/(x1-x0)*pi) )
        else
            yi = y0
        end if
    end function

    pure function ip_zero_one( x0,y0, x1,y1, xi ) result(yi)
        real, intent(in) :: x0,y0,x1,y1,xi
        real :: yi
        if ((y0 == 0.) .and. (y1 == 0.)) then
            yi = 0.     + 0.*(x0+x1+xi)  ! zero ! (to get rid of compiler warnings)
        else
            yi = 1.
        end if

    end function
    
end module
