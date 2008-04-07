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

program test_heap

    use better_varying_string
    use util
    use heap
    
    implicit none
    
    integer, parameter :: nn = 1000000
    real, dimension(nn) :: keys
    integer, dimension(nn) :: back
    
    type(t_index_heap) :: h
    integer :: i,j
    
    call test_begin("test_heap")
    
    call initheap( h, nn )
    
    back(:) = 0
    do i=1,nn
        keys(i) = nn-(i-1)
    end do
    
    do i=1,nn
        call pushheap(h,i,keys,back)
    end do
    
    do i=1,nn/1000
        if (back(h%iheap(i)) .ne. i) then
            call test_fail("backpointer "//i)
            exit
        end if
    end do
    
    do i=1,nn/1000
        call popheap(h,j,keys,back)
        if (abs( nn-nn+i - keys(j)) .gt. 0.00001) then
            call test_fail("popheap "//i)
            exit
        end if
    end do
    
    call destroyheap( h )
    
    call test_end()
    call cleanup()
    
end program
