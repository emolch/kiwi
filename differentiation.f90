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

module differentiation

    implicit none
    
    private
    
    public :: differentiate
    
    contains
    
    subroutine differentiate( dt, f, df, method )
    
        real, intent(in)                     :: dt
        real, dimension(:), intent(in)       :: f
        real, dimension(:), intent(out)      :: df
        integer, optional, intent(in)        :: method

    ! this subroutine differentiates a sampled function f given at 
    ! equitistantly spaced points with the spacing dt and
    ! puts the result to df.
    
    ! the optional argument method is currently ignored.
    
    ! the differentiation is currently done by doing central differences:
    !        f'(t) = (f(t+dt)-f(t-dt)) / (2 dt)
    ! boundary points are done using forward or backward differences
        
        integer :: i,j,il,iu,jl,ju
        
        if (present(method)) return
        
        if (size(f) /= size(df)) then
            stop "sizes of input and output array for differentiation do not match."
        end if
        
        if (size(f) < 2) then
            stop "sizes of arrays for differentiation are too short."
        end if
        
        jl = lbound(df,1)
        ju = ubound(df,1)
        il = lbound(f,1)
        iu = ubound(f,1)
        
        j=jl+1
        do i=il+1,iu-1
            df(j) = ( f(i+1) - f(i-1) ) / (dt*2)
            j = j+1
        end do
        
        df(jl) = ( f(il+1)-f(il) ) / dt
        df(ju) = ( f(iu)-f(iu-1) ) / dt
    
    end subroutine differentiate
    

end module differentiation
