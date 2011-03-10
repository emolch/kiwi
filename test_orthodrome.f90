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

program test_orthodrome

    use util
    use constants
    use orthodrome

    implicit none

    type(t_geo_coords) :: a,b
    real*8 :: azi,dist,bazi
    real :: km
        
    call test_begin("test_orthodrome")
    
    a%lat = 0.
    a%lon = 0.
    
    b%lat = 90.
    b%lon = 0.
    if (.not. near(0D00, r2d(azimuth( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("azimuth 1")
    end if
    if (.not. near(90D00, r2d(arcdistance( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("arcdistance 1")
    end if
    
    b%lat = 0.
    b%lon = 100.
    if (.not. near(90D00, r2d(azimuth( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("azimuth 2")
    end if
    if (.not. near(100D00, r2d(arcdistance( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("arcdistance 2")
    end if
    
    a%lon = 10.
    b%lat = 0.
    b%lon = -10.
    if (.not. near(-90D00, r2d(azimuth( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("azimuth 3")
    end if
    if (.not. near(20D00, r2d(arcdistance( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("arcdistance 3")
    end if
    
    a%lat = 10.
    a%lon = 0.
    b%lat = -10.
    b%lon = -0.001
    if (.not. near(-180D00, r2d(azimuth( d2r(a), d2r(b) )), 0.1D00)) then
        call test_fail("azimuth 4")
    end if 
    if (.not. near(20D00, r2d(arcdistance( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("arcdistance 4")
    end if
    
    a%lat = 10.
    a%lon = 0.
    b%lat = -10.
    b%lon = 0.
    if (.not. near(180D00, r2d(azimuth( d2r(a), d2r(b) )), 0.1D00)) then
        call test_fail("azimuth 5")
    end if
    if (.not. near(20D00, r2d(arcdistance( d2r(a), d2r(b) )), 0.001D00)) then
        call test_fail("arcdistance 5")
        print *, r2d(arcdistance( d2r(a), d2r(b) ))
    end if

    km = 1000.
    
  ! a real world test: distance hamburg-munich
    a%lat = 53.556867
    a%lon = 9.994622
    b%lat = 48.139743
    b%lon = 11.560050
    dist = distance_accurate50m( d2r(a), d2r(b) )
    if (.not. near(dist/km, 612.59D00, 0.05D00)) then
        print *, dist
        call test_fail( "hamburg-munich: " )
    end if
    
  ! geomatikum-home
    a%lat = 53.568391
    a%lon = 9.973672
    b%lat = 53.580049
    b%lon = 9.956507
    dist = distance_accurate50m( d2r(a), d2r(b) )
    if (.not. near(dist/km, 1.73D00, 0.05D00)) then
        call test_fail("geomatikum-home")
    end if
    
  ! a real world test: distance berlin-tokio
    a%lat = 52.5167
    a%lon = 13.4000
    b%lat = 35.7000
    b%lon = 139.7667
    dist = distance_accurate50m( d2r(a), d2r(b) )
    if (.not. near(dist/km, 8941.20671D00, 0.05D00)) then
        call test_fail( "berlin-tokio")
    end if
    
    
    call approx_differential_azidist( 10.*km,10.*km, &
                                      d2r(90D00), d2r(-90D00), 20D00*km, &
                                      azi, bazi, dist )
                                      
    if (.not. (near(r2d(azi),135D00,1D00) .and. &
        near(dist,10D00*km*sqrt(2.),1D00))) then
        call test_fail("differential azidist 1")
    end if
    
    call approx_differential_azidist( 10.*km,10.*km, &
                                      d2r(90D00), d2r(-90D00), 1111.95D00*km, &
                                      azi, bazi, dist )
    if (.not. (near(r2d(azi),90D00,1D00) .and. &
        near(dist,1111.95D00*km-10D00*km,50D00))) then
        call test_fail("differential azidist 2")
    end if
    
    call test_end()
    call cleanup()
    
contains
    
    elemental logical function near( a,b, eps )
        real*8, intent(in) :: a,b,eps
        near = abs(a-b) < eps
    end function
    
end program
