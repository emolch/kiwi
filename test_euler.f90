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


program test_euler

    use util
    use euler
    use constants
    
    implicit none
    
    real, dimension(3,3) :: rot
    real, dimension(3,3) :: sys = reshape( (/ 1,0,0, 0,1,0, 0,0,1 /), (/3,3/) )
    real, dimension(3,3) :: rotsys, rotalpha, rotbeta, rotgamma
    real :: eps = 0.001
    integer :: i

    call test_begin("test_euler")
    
    ! expected results for pi/2 rotations
    rotbeta = reshape( (/ 0,1,0, -1,0,0, 0,0,1 /), (/3,3/) )
    rotalpha = reshape( (/ 1,0,0, 0,0,1, 0,-1,0 /), (/3,3/) )
    rotgamma = reshape( (/ 0,1,0, -1,0,0, 0,0,1 /), (/3,3/) )

    
    call init_euler( pi/2., 0., 0., rot )
    do i=1,3
        rotsys(:,i) = matmul(rot,sys(:,i))
    end do
    if (any( .not. near(rotsys,rotalpha,eps) )) call test_fail("alpha")
    
    
    call init_euler( 0., pi/2., 0., rot )
    do i=1,3
        rotsys(:,i) = matmul(rot,sys(:,i))
    end do
    if (any( .not. near(rotsys,rotbeta,eps) )) call test_fail("beta")
    
    
    call init_euler( 0., 0., pi/2., rot )
    do i=1,3
        rotsys(:,i) = matmul(rot,sys(:,i))
    end do
    if (any( .not. near(rotsys,rotgamma,eps) )) call test_fail("gamma")

    
    call test_end()

    call cleanup()
    
contains
    
    elemental logical function near( a,b, eps )
        real, intent(in) :: a,b,eps
        near = abs(a-b) < eps
    end function
            
end program
