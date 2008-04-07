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

program test_eikonal

    use util
    use heap
    use eikonal
    
    implicit none

    integer, parameter :: nx = 500
    integer, parameter :: ny = 1000
    
    real, dimension(nx,ny) :: speed, times
    
    real, dimension(2) :: delta, initialpoint
    real :: eps
    
    call test_begin("test_eikonal")
    
    
    delta(:) = (/50./nx, 50./ny/)
    speed(:,:) = 2.
    initialpoint(1) = 0.
    initialpoint(2) = 25.
    
    eps = max(delta(1),delta(2))/minval(speed)
    
    call eikonal_solver_fmm(speed, (/0.,0./), delta, initialpoint, times )
    
    if (.not. near(times(1,1), 12.5, eps) ) then
        call test_fail("1,1")
    end if
    if (.not. near(times(1,ny), 12.5, eps) ) then
        call test_fail("1,ny")
    end if
    if (.not. near(times(nx,1), 27.95, eps) ) then
        call test_fail("nx,1")
    end if
    if (.not. near(times(nx,ny), 27.95, eps) ) then
        call test_fail("nx,ny")
    end if
        
    !call dumpfield( "chocko", times )
    
    call test_end()

    call cleanup()
  
  contains
  
    elemental logical function near( a,b, eps )
        real, intent(in) :: a,b,eps
        near = abs(a-b) < eps
    end function
    
    subroutine dumpfield( fn, a )

        character(len=*) :: fn
        real, dimension(:,:), intent(in) :: a

        integer :: i,j

        open(10,file=fn,status='unknown')

        do j=1,size(a,2)
            do i=1,size(a,1)
                write(10,*) i,j, a(i,j)
            end do
        end do

        close(10)

   end subroutine
  
end program
