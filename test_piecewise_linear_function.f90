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

program test_piecewise_linear_function
    
    use util
    use piecewise_linear_function

    implicit none

    type(t_plf) :: func
    real :: a,c
        
    call test_begin("test_piecewise_linear_function")
    
    call plf_make( func, 0.,0., 1.,1., 2.,1., 3.,0. )
    if (2. .ne. plf_integrate( func, -1.,3. )) then
        call test_fail("1")
    end if
    
    if (1. .ne. plf_integrate( func, -1.,1.5 )) then
        call test_fail("2")
    end if
    
    if (6./8.+1. .ne. plf_integrate( func, 0.5,2.5 )) then
        call test_fail("3")
    end if
    
    if (3./8. .ne. plf_integrate( func, 2.,2.5 )) then
        call test_fail("4")
    end if
    
    if (0. .ne. plf_integrate( func, 2.,2. )) then
        call test_fail("5")
    end if
    
    if (0. .ne. plf_integrate( func, 2.5,2.5 )) then
        call test_fail("6")
    end if
    
    if (1./8.-1./32. .ne. plf_integrate( func, 2.5,2.75 )) then
        call test_fail("7")
    end if
    
    if (1. .ne. plf_integrate( func, 1.,2. )) then
        call test_fail("8")
    end if
    
    call plf_integrate_and_centroid( func, -1., 6., a, c )
    if (a .ne. 2. .and. c .ne. 1.5) then
        call test_fail("9")
    end if
    
    call plf_integrate_and_centroid( func, 0., 0.5, a, c )
    if (a .ne. 1./8. .and. c .ne. 1./3.) then
        call test_fail("10")
    end if
    
    call plf_integrate_and_centroid( func, 0., 2., a, c )
    if (a .ne. 3./2. .and. c .ne. 1. + 2./9.) then
        call test_fail("11")
    end if
    
    
    
    call test_end()
    call cleanup()
    
end program
