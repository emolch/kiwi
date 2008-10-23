! $Id: orthodrome.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


module orthodrome

  ! This module provides methods to calculate azimuths and distances between
  ! points on a sphere or ellipsoid of rotation.
  
  ! It also provides a method to calculate relative azimuths and distances for
  ! small perturbations of the point coordinates.
  
  ! All angles and coordinates in the subroutines are required in radians.
  
    use constants
    use util
    
    implicit none
    
   type, public :: t_geo_coords
        real*8 :: lat, lon
    end type
    
    private
    public approx_differential_azidist
    public arcdistance
    public distance
    public distance_accurate50m
    public azimuth
    public azibazi
    public azidist
    
    public r2d
    public d2r
    
  ! degrees to radians
    interface d2r
        module procedure d2r_tgc
        module procedure d2r_r
        module procedure d2r_d
    end interface
  
  ! radians to degrees
    interface r2d
        module procedure r2d_tgc
        module procedure r2d_r
        module procedure r2d_d
    end interface
    
  ! these parameters control what approximation is used when calculating relative
  ! azimuths and distances:
  !
  ! flat geometry approximation is used for distances up to this:
  !  real*8, parameter, public :: max_distance_flat_approx = 250000. ! [m]
    real*8, parameter, public :: max_distance_flat_approx = -1. ! [m] ! turned off.
   
   
  ! constant azimuth approximation is used when distance/dr ratio is greater than this:
  !  real*8, parameter, public :: min_ratio_const_azimuth_approx = 50.
    real*8, parameter, public :: min_ratio_const_azimuth_approx = huge(min_ratio_const_azimuth_approx) ! turned off.
    
    
  contains
    
    subroutine approx_differential_azidist( delta_x, delta_y, azimuth, backazimuth, dist,&
                                      new_azimuth, new_backazimuth, new_dist )
    
        real, intent(in) :: delta_x, delta_y ! north, east
        real*8, intent(in) :: azimuth, backazimuth, dist
        real*8, intent(out) :: new_azimuth, new_backazimuth, new_dist
        
      ! given two points a and b connected with distance 'dist' in [m]
      ! and azimuth 'azimuth', calculate azimuth and arc-distance to b
      ! from new point with relative cartesian coordinates r=(delta_x,delta_y)
      ! the new azimuth and distance are put to 'new_dist'
      ! and 'new_azimuth'
      
      ! the following two approximations are tried in this order:
      
      ! 1) if dist < max_distance_flat_approx
      !    the calculation is performed assuming a flat geometry
      
      ! 2) if r/dist > min_ratio_const_azimuth_approx
      !    it is assumed that the azimuth does not change
               
      ! 3) exact formula for spherical earth are used

      
        real*8 :: new_dist_x, new_dist_y
        real*8 :: r, a,b,c
        real*8 :: alpha, beta, gamma, lambda
      
        if (dist < max_distance_flat_approx) then
          ! flat approximation for local distances
            new_dist_x = dist*cos(azimuth)-delta_x
            new_dist_y = dist*sin(azimuth)-delta_y
            new_azimuth = atan2(new_dist_y, new_dist_x)
            new_backazimuth = backazimuth + (new_azimuth-azimuth)
            new_dist = sqrt(new_dist_x**2 + new_dist_y**2)
        else
        
            r = sqrt(delta_x**2 + delta_y**2)
                        
            if (dist/r > min_ratio_const_azimuth_approx) then
              ! 'const-azimuth'-approximation for larger distances
                new_azimuth = azimuth
                new_backazimuth = backazimuth
                new_dist = dist - (delta_x*cos(azimuth) + delta_y*sin(azimuth))
            else 
              ! "exact" formulas assuming spherical earth
                a = r/earthradius
                b = dist/earthradius
                lambda = atan2(delta_y,delta_x)
                gamma = azimuth-lambda
                
                c = acos( clip(cos(a)*cos(b)+sin(a)*sin(b)*cos(gamma),-1D00,1D00) )
                alpha = asin( clip(sin(a)*sin(gamma)/sin(c),-1D00,1D00) )
                beta = asin( clip(sin(b)*sin(gamma)/sin(c),-1D00,1D00) )
                
                ! put alpha and beta to the correct quandrant...
                if (cos(a)-cos(b)*cos(c) < 0) then
                    if (alpha > 0) then
                        alpha = pi_-alpha
                    else 
                        alpha = -pi_-alpha
                    end if 
                end if
                if (cos(b)-cos(a)*cos(c) < 0) then
                    if (beta > 0) then
                        beta = pi_-beta
                    else 
                        beta = -pi_-beta
                    end if
                end if
                
                new_dist = c*earthradius
                new_backazimuth = wrap(backazimuth + alpha,-pi_,pi_)
                new_azimuth = wrap(lambda - pi_ - beta,-pi_,pi_)
                
            end if
            
        end if
        
    end subroutine
    
    elemental function clip( x, mi, ma )
        real*8, intent(in) :: x, mi, ma
        real*8 :: clip
        
        clip = min(max(mi,x),ma)
        
    end function
    
    elemental function wrap( x, mi, ma )
        real*8, intent(in) :: x, mi, ma
        real*8 :: wrap
        wrap = x - floor((x-mi)/(ma-mi)) * (ma-mi)
    end function
    
    elemental function arcdistance( a,b )
        
        type(t_geo_coords), intent(in) :: a,b
        real*8 :: arcdistance
      
      ! returns arc-distance in radians between points a and b.
      ! coordinates must be given in radians.
        
        arcdistance = acos(cosdelta(a,b))
        
    end function
    
    elemental function distance( a, b )
       
       type(t_geo_coords), intent(in) :: a,b
       real*8 :: distance
    
        distance = arcdistance(a, b)*earthradius
    
    end function
    
    elemental function distance_accurate50m( a, b )
    
      ! more accurate distance calculation based on a spheroid of rotation
      
      ! returns distance in [m] between points a and b
      ! coordinates must be given in radians.
      
      ! should be accurate to 50 m using WGS84
      !     earth_oblateness = 1./298.257223563
      !     earthradius_equator = 6378.14 * 1000.

      ! from wikipedia :  http://de.wikipedia.org/wiki/Orthodrome
      ! based on: Meeus, J.: Astronomical Algorithms, S 85, Willmann-Bell, 
      !           Richmond 2000 (2nd ed., 2nd printing), ISBN 0-943396-61-1
       
        type(t_geo_coords), intent(in) :: a,b
        real*8 :: distance_accurate50m
    
        real*8 :: f,g,l,s,c,w,r,d,h1,h2
        
        f = (a%lat + b%lat) / 2.
        g = (a%lat - b%lat) / 2.
        l = (a%lon - b%lon) / 2.
        
        s = sin(g)**2 * cos(l)**2 + cos(f)**2 * sin(l)**2
        c = cos(g)**2 * cos(l)**2 + sin(f)**2 * sin(l)**2

        w = atan( sqrt( s/c ) )
        r = sqrt(s*c)/w
        d = 2.*w*earthradius_equator
        h1 = (3.*r-1.)/(2.*c)
        h2 = (3.*r+1.)/(2.*s)
        
        distance_accurate50m = d * (1.+ earth_oblateness * h1 * sin(f)**2 * cos(g)**2 - &
                                        earth_oblateness * h2 * cos(f)**2 * sin(g)**2)
        
    end function
    
    elemental function azimuth( a, b )
        
        type(t_geo_coords), intent(in) :: a,b
        real*8 :: azimuth
        
      ! get azimuth of point b as seen at a.
      ! coordinates must be given in radians.
      ! result is in the range ]-pi,pi]
        
        azimuth = atan2( cos(a%lat) * cos(b%lat) * sin( b%lon - a%lon ), &
                         sin(b%lat) - sin(a%lat) * cosdelta(a,b) )
                         
    end function
    
    elemental subroutine azibazi( a, b, azimuth, backazimuth )
        
        type(t_geo_coords), intent(in) :: a,b
        real*8, intent(out) :: azimuth, backazimuth
        
        real*8 :: t, sb, sa, cd
        
      ! get azimuth of point b as seen at a.
      ! and azimuth of point a as seen at b (backazimuth).
      ! coordinates must be given in radians.
      ! result is in the range ]-pi,pi]
      
        t = cos(a%lat) * cos(b%lat) * sin( b%lon - a%lon )
        sb = sin(b%lat)
        sa = sin(a%lat)
        cd = cosdelta(a,b)
        
        azimuth     = atan2(  t, sb - sa * cd )
        backazimuth = atan2( -t, sa - sb * cd )
        
    end subroutine
    
    elemental subroutine azidist( a,b, azimuth, arcdist )
    
        type(t_geo_coords), intent(in) :: a,b
        real*8,intent(out) :: azimuth, arcdist
        
      ! get arc-distance and azimuth in one shot
      ! coordinates must be given in radians.
      
        real*8 :: cd
        
        cd = cosdelta(a,b)
        
        azimuth = atan2( cos(a%lat) * cos(b%lat) * sin( b%lon - a%lon ), &
                         sin(b%lat) - sin(a%lat) * cd )
        arcdist = acos(cd)
        
    end subroutine
    
    elemental function cosdelta( a, b )

        type(t_geo_coords), intent(in) :: a,b
        real*8 :: cosdelta
        
      ! needed by azimuth and arcdistance
    
        cosdelta = sin(a%lat) * sin(b%lat) + cos(a%lat) * cos(b%lat) * cos(b%lon-a%lon)
        
    end function
     
    elemental function d2r_tgc( deg ) result(rad)
    
        type(t_geo_coords), intent(in) :: deg
        type(t_geo_coords) :: rad
        
        rad%lat = 2./360.*pi*deg%lat
        rad%lon = 2./360.*pi*deg%lon
        
    end function
    
    elemental function r2d_tgc( rad ) result(deg)
    
        type(t_geo_coords), intent(in) :: rad
        type(t_geo_coords) :: deg
        
        deg%lat = 360./2./pi*rad%lat
        deg%lon = 360./2./pi*rad%lon
        
    end function
    
    elemental function d2r_r( deg ) result(rad)
    
        real, intent(in) :: deg
        real :: rad
        
        rad = 2./360.*pi*deg
        
    end function
    
    elemental function r2d_r( rad ) result(deg)
    
        real, intent(in) :: rad
        real :: deg
        
        deg = 360./2./pi*rad
        
    end function
    
    elemental function d2r_d( deg ) result(rad)
    
        real*8, intent(in) :: deg
        real*8 :: rad
        
        rad = 2./360.*pi*deg
        
    end function
    
    elemental function r2d_d( rad ) result(deg)
    
        real*8, intent(in) :: rad
        real*8 :: deg
        
        deg = 360./2./pi*rad
        
    end function
   
end module
