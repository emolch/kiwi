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

program eikonal_benchmark

    use better_varying_string
    use util
    use heap
    use eikonal
    
    implicit none

    integer :: nx
    integer :: ny
    
    real, dimension(:,:), allocatable :: speed, times
    
    real, dimension(2) :: delta, initialpoint
    real :: time_begin, time_end, duration
    
    g_pn = 'eikonal_benchmark'  

    delta(:) = (/1.,1./)
    initialpoint(1) = 0.
    initialpoint(2) = 0.

    do nx=100,2000,100
        do ny=nx,2000,100

            if (allocated(speed)) deallocate(speed)
            if (allocated(times)) deallocate(times)
            allocate(speed(nx,ny))
            allocate(times(nx,ny))
            speed = 2.0
            call cpu_time( time_begin )
            call eikonal_solver_fmm(speed, (/0.,0./), delta, initialpoint, times )
            call cpu_time( time_end )
            duration = time_end - time_begin
            print *, nx, ny, duration
        end do
    end do

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
