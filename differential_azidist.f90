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


program differential_azidist

    use util
    use constants
    use orthodrome
    use better_varying_string
    use varying_string_getarg

    implicit none

    
    type(t_geo_coords) :: a,b,r
    real :: northab, eastab, northba, eastba
    real*8 :: distab, distar, distbr
    real*8 :: aziab, aziba, aziar, azira, azibr, azirb
    real*8 :: approx_distbr, approx_azibr, approx_azirb
    real*8 :: approx_distar, approx_aziar, approx_azira
    integer :: ilat, ilon, nlat, nlon
    real :: dlat, dlon
    
    a%lat = 53.556867
    a%lon = 9.994622
    b%lat = 53.139743
    b%lon = 9.560050
    
    nlat = 180
    nlon = 360
    dlat = 180./nlat
    dlon = 360./nlon
    do ilat = 1,nlat
        do ilon = 1,nlon
            r%lat = -90+(ilat-1)*dlat
            r%lon = -180+(ilon-1)*dlon
    !call vs_getarg(1, str)
    !r%lat = str
    !call vs_getarg(2, str)
    !r%lon = str
    
            distab = distance_accurate50m( d2r(a), d2r(b) )
            distar = distance_accurate50m( d2r(a), d2r(r) )
            distbr = distance_accurate50m( d2r(b), d2r(r) )
            call azibazi( d2r(a), d2r(b), aziab, aziba )
            call azibazi( d2r(a), d2r(r), aziar, azira )
            call azibazi( d2r(b), d2r(r), azibr, azirb )
            
            northab = real(cos(aziab)*distab)
            eastab = real(sin(aziab)*distab)
            call approx_differential_azidist( northab, eastab, &
                                            aziar, azira, distar, &
                                            approx_azibr, approx_azirb, approx_distbr )
            !print *, "distance between b and r: ", distbr, approx_distbr, (distbr-approx_distbr)
            !print *, "azimuth to r at b: ", r2d(azibr), r2d(approx_azibr), r2d(azibr-approx_azibr) 
            !print *, "azimuth to b at r: ", r2d(azirb), r2d(approx_azirb), r2d(azirb-approx_azirb)
            
            northba = real(cos(aziba)*distab)
            eastba = real(sin(aziba)*distab)
            call approx_differential_azidist( northba, eastba, &
                                            azibr, azirb, distbr, &
                                            approx_aziar, approx_azira, approx_distar )
            !print *, "distance between a and r: ", distar, approx_distar, (distar-approx_distar)
            !print *, "azimuth to r at a: ", r2d(aziar), r2d(approx_aziar), r2d(aziar-approx_aziar) 
            !print *, "azimuth to a at r: ", r2d(azira), r2d(approx_azira), r2d(azira-approx_azira)
            
            print *, r%lat, r%lon, (distar-approx_distar), r2d(piwrap(aziar-approx_aziar)), &
            r2d(piwrap(azira-approx_azira)), (distbr-approx_distbr), r2d(piwrap(azibr-approx_azibr)) &
            , r2d(piwrap(azirb-approx_azirb))
!            print *, r%lat, r%lon,     distance_accurate50m( d2r(a), d2r(r) ) - distance(d2r(a), d2r(r) )
    
        end do
    end do
    
  contains
    
    elemental function wrap( x, mi, ma )
        real*8, intent(in) :: x, mi, ma
        real*8 :: wrap
        wrap = x - floor((x-mi)/(ma-mi)) * (ma-mi)
    end function
    
    
    elemental function piwrap( x )
        real*8, intent(in) :: x
        real*8 :: piwrap
        piwrap = x - floor((x+pi_)/(2.*pi_)) * 2.*pi_
    end function
    
    
    
    
    
end program