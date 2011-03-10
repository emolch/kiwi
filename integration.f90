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

module integration

    implicit none
    
    private
    
    public :: antiderivate, integrate
    
    contains
    
    pure subroutine antiderivate( dt, f, ff, method )
    
        real, intent(in)                     :: dt
        real, dimension(:), intent(in)       :: f
        real, dimension(:), intent(out)      :: ff
        integer, optional, intent(in)        :: method

    ! this subroutine integrates a sampled function f given at 
    ! equitistantly spaced points with the spacing dt and
    ! puts the result to ff.
    
    ! the optional argument method is currently ignored.
        
    ! the integration is simply done by
    ! ff(i+1) = ff(i) + (f(i+1) + f(i))/2 * dt
        
        integer :: i,n
        
        if (present(method)) return

        n = min(size(f),size(ff))
        
        if (size(f) < 2) then
            ff(:) = 0.0
            return
        end if
        
        ff(1) = 0.0
        do i=1,n-1
            ff(i+1) = ff(i) + (f(i+1) + f(i))/2 * dt
        end do
        
    end subroutine antiderivate
    
    pure function integrate( dt, f, method ) result( ff )
        
        real, intent(in)                 :: dt
        real, dimension(:), intent(in)   :: f
        integer, optional, intent(in)    :: method
        real :: ff
        
        integer :: i,il,iu
        
        if (present(method)) return
     
        if (size(f) < 2) then
            ff = 0.0
            return
        end if
        
        il = lbound(f,1)
        iu = ubound(f,1)
        
        ff = (f(il)+f(iu))/2.0
        do i=il+1,iu-1
            ff = ff + f(i)
        end do
        ff = ff * dt
        
        
    end function integrate

end module integration
