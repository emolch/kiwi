! $Id: test_comparator.f90 669 2007-09-14 12:26:15Z sebastian $ 
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


program test_comparator


    use util
    use sparse_trace
    use comparator
    use piecewise_linear_function
    
    implicit none

    type(t_strip) :: s1, s2, stap
    type(t_probe) :: a, b
    real :: n, eps, dt
    integer, dimension(2) :: shiftrange
    real, dimension(-5:5) :: cross_corr
    type(t_plf) :: taper

    dt = 1.
    call probe_init( a,dt )
    call probe_init( b,dt )
    
    call test_begin("test_comparator")

    call strip_init( (/-1,2/), (/0.,0.,5.,1./), s1 )
    call strip_init( (/1,4/), (/5.,1.,1.,1./), s2 )
    
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    n = probes_norm( a, b )
    if (n /= 0.) call test_fail("1")
    n = probes_norm( a, b, L1NORM )
    if (n /= 0.) call test_fail("1b")
    
    call strip_init( (/-1,2/), (/0.,0.,0.,1./), s1 )
    call strip_init( (/1,4/), (/0.,1.,1.,1./), s2 )
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    n = probes_norm( a, b )
    if (n /= 0.) call test_fail("2")
    n = probes_norm( a, b, L1NORM )
    if (n /= 0.) call test_fail("2b")
    
    eps = 0.000001;
    
    call strip_init( (/-4,-1/), (/1.,0.,0.,0./), s1 )
    call strip_init( (/1,4/), (/1.,0.,0.,0./), s2 )
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    n = probes_norm( a,b )
    if (.not. near(n,sqrt(2.),eps)) call test_fail("3")
    n = probes_norm( a,b, L1NORM )
    if (.not. near(n,2.,eps)) call test_fail("3b")
    
    call strip_init( (/0,3/), (/1.,2.,1.,0./), s1 )
    call strip_init( (/1,4/), (/1.,1.,0.,1./), s2 )
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    n = probes_norm( a,b )
    if (.not. near(n ,sqrt(3.),eps)) call test_fail("4")
    n = probes_norm( a,b, L1NORM )
    if (.not. near(n ,3.,eps)) call test_fail("4b")
    
    call strip_init( (/0,3/), (/0.,1.,2.,1./), s1 )
    call strip_init( (/1,2/), (/1.,2./), s2 )
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    n = probes_norm( a,b )
    if (.not. near(n,1.,eps)) call test_fail("5")
    n = probes_norm( a,b, L1NORM )
    if (.not. near(n,1.,eps)) call test_fail("5b")

    
    call strip_init( (/0,4/), (/0.,1.,2.,1.,0./), s1 )
    call strip_init( (/10,14/), (/0.,1.,2.,1.,0./), s2 )
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    n = probes_norm( a,b, AMPSPEC_L2NORM )
    if (.not. near(n,0.,eps)) call test_fail("6")

    call strip_init( (/1,5/), (/0.,1.,2.,1.,0./), s1 )
    call strip_init( (/2,6/), (/0.,1.,2.,1.,0./), s2 )
    call probe_set_array( a, s1 )
    call probe_set_array( b, s2 )
    shiftrange = (/-5,5/)
    call probes_windowed_cross_corr( a, b, shiftrange, cross_corr )
    if (any(cross_corr .ne. (/0., 0., 1., 4., 6., 4., 1., 0., 0., 0., 0./))) then
        call test_fail("cross correlation")
    end if

    call plf_make( taper, 2.5,1., 3.5,1. )
    call probe_set_taper( a, taper )
    call probe_get_tapered(a, stap)
    call probe_set_taper( b, taper )
    call probes_windowed_cross_corr( a, b, shiftrange, cross_corr )
    if (any(cross_corr .ne. (/0., 0., 0., 2., 4., 2., 0., 0., 0., 0., 0./))) then
        call test_fail("tapered cross correlation")
    end if

    call test_end()
    
    call cleanup()
    
  contains
      
    elemental logical function near( a,b, eps )
        real, intent(in) :: a,b,eps
        near = abs(a-b) < eps
    end function

    
    
end program
