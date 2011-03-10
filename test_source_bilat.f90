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

program test_source_bilat

    use util
    use discrete_source
    use constants
    use orthodrome
    
    use source
    
    implicit none
  
  ! make a simple source with rupturing from west to east, with a dip of 45 deg.
  ! hanging wall to south.
  ! length 2 km; width 1 km; depth 1km; rupture velocity: 3000 m/s;  risetime 1 s
  ! slip perpendicular to rupture direction, thrust fault
  
    integer :: i
    type(t_psm)        :: psm
    type(t_tdsm)       :: tdsm
    real :: mxx,myy,mzz,mxy,mxz,myz, epsm, msum
    logical :: thrusttestfail, ok
    double precision :: ref_time
    
    
    call test_begin("test_source_bilat")
    
    ref_time = 0.0
    
    call psm_set_origin_and_time( psm, t_geo_coords(0.,0.), ref_time  )
    
    call psm_set( psm, psm_bilat, &
                              (/    0.,   & ! time
                                    0.,    & ! north of reference
                                    0.,    & ! east of reference
                                    1000., & ! depth
                                    1.,    & ! moment
                                    90., & ! strike
                                    45., & ! dip
                                    90.,    & ! rake
                                    0., & ! rupture direction
                                    2000., & ! length_a
                                    0.,    & ! length_b
                                    1000., & ! width
                                    2000., & ! rupvel
                                    1. /) )   ! risetime
    
    call psm_to_tdsm( psm, tdsm, 0.5, ok )
    
    epsm = 1./ size(tdsm%centroids) / 100
    msum = 0.
    thrusttestfail = .false.
    do i=1,size(tdsm%centroids)
   
        mxx = tdsm%centroids(i)%m(1)
        myy = tdsm%centroids(i)%m(2)
        mzz = tdsm%centroids(i)%m(3)
        mxy = tdsm%centroids(i)%m(4)
        mxz = tdsm%centroids(i)%m(5)
        myz = tdsm%centroids(i)%m(6)
                
        msum = msum + mxx
      ! 45deg dip has has -mxx == mzz all others zero
        if ( near(mxx,0.,epsm) .or. &
             near(mzz,0.,epsm) .or. &
             .not. near(mxx,-mzz,epsm) .or. &
             mxx > 0. .or. &
             .not. near(myy,0.,epsm) .or. &
             .not. near(mxy,0.,epsm) .or. &
             .not. near(mxz,0.,epsm) .or. &
             .not. near(myz,0.,epsm) ) then
            thrusttestfail = .true.
        end if
    end do
    if (thrusttestfail) call test_fail( "thrust1")
    if (.not. near(-1.,msum,0.01)) &
        call test_fail("thrustsum1")
    
    
    
    
    call psm_set( psm, psm_bilat, &
                              (/    0.,   & ! time
                                    0.,    & ! north of reference
                                    0.,    & ! east of reference
                                    1000., & ! depth
                                    1.,    & ! magnitude
                                    45., & ! strike
                                    90., & ! dip
                                    0.,    & ! rake
                                    0., & ! rupture direction
                                    2000., & ! length_a
                                    0.,    & ! length_b
                                    1000., & ! width
                                    2000., & ! rupvel
                                    1. /) )   ! risetime
    
    call psm_to_tdsm( psm, tdsm, 0.5, ok )
    
    
    
    epsm = 1./ size(tdsm%centroids) / 100
    msum = 0.
    thrusttestfail = .false.
    do i=1,size(tdsm%centroids)
   
        mxx = tdsm%centroids(i)%m(1)
        myy = tdsm%centroids(i)%m(2)
        mzz = tdsm%centroids(i)%m(3)
        mxy = tdsm%centroids(i)%m(4)
        mxz = tdsm%centroids(i)%m(5)
        myz = tdsm%centroids(i)%m(6)
               
        msum = msum + mxx
      ! 90deg dip with strike 45 has -mxx == myy all others zero
        if ( near(mxx,0.,epsm) .or. &
             near(myy,0.,epsm) .or. &
             .not. near(mxx,-myy,epsm) .or. &
             mxx > 0. .or. &
             .not. near(mzz,0.,epsm) .or. &
             .not. near(mxy,0.,epsm) .or. &
             .not. near(mxz,0.,epsm) .or. &
             .not. near(myz,0.,epsm) ) then
            thrusttestfail = .true.
        end if
    end do
    if (thrusttestfail) call test_fail( "thrust2")
    if (.not. near(-1.,msum,0.01)) &
        call test_fail("thrustsum2")
    
    call test_end()
    call cleanup()
    
  contains
  
    elemental logical function near( a,b, eps )
        real, intent(in) :: a,b,eps
        near = abs(a-b) < eps
    end function

end program
