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

module eikonal

    use util
    use heap    

    implicit none
    
    private 
    public eikonal_solver_fmm
    
  contains

    subroutine eikonal_solver_fmm( speed, origin, delta, initialpoint, times )
      
      ! Fast marching level set method solver (sethian 1996)
  
      ! On a 2D-grid with spacing delta(2) and speeds 'speed' and a signal 
      ! starting at initialpoint(2), calculate arrival times 'times'.
      
        real, dimension(:,:), intent(in)  :: speed
        real, dimension(2),   intent(in)  :: origin
        real, dimension(2),   intent(in)  :: delta
        real, dimension(2),   intent(in)  :: initialpoint
        real, dimension(:,:), intent(out) :: times
        
        integer :: nalive
        
        integer, dimension(size(speed,1),size(speed,2)) :: backpointers
        
      ! backpointers are set to two different invalid values:
        integer, parameter :: FARAWAY = -1 
        integer, parameter :: ALIVE = 0
        
        type(t_index_heap) :: heap
        
        real :: dx,dy
        integer :: ix, iy, nx, ny
        integer :: imin
        real :: infinity
        
        infinity = huge(infinity) * 0.1
        
        nx = size(speed,1)
        ny = size(speed,2)
        dx = delta(1)
        dy = delta(2)
        
        if (size(times,1) /= nx .or. size(times,2) /= ny) then
            call die("eikonal_solver_fmm(): sizes of times and speed arrays do not match.")
        end if
        
      ! initialize 
        
        backpointers(:,:) = FARAWAY
        
        ix = int((initialpoint(1)-origin(1))/dx) + 1
        iy = int((initialpoint(2)-origin(2))/dy) + 1
        
        if (ix < 1) ix = 1
        if (nx < ix) ix = nx
        if (iy < 1) iy = 1
        if (ny < iy) iy = ny
        
        times(:,:) = infinity
      
      ! alive points
      
        times(ix,iy) = 0.0
        if (nx==1 .and. ny==1) return
        
        backpointers(ix,iy) = ALIVE
        nalive = 1
      
      ! narrowband points
      
        call initheap( heap, nx*ny )
        
        if (1 < ix)  times(ix-1,iy) = dx/speed(ix-1,iy)
        if (ix < nx) times(ix+1,iy) = dx/speed(ix+1,iy)
        if (1 < iy)  times(ix,iy-1) = dy/speed(ix,iy-1)
        if (iy < ny) times(ix,iy+1) = dy/speed(ix,iy+1)
                
        if (1 < ix)  call pushheap( heap, ind(ix-1,iy), times, backpointers )
        if (ix < nx) call pushheap( heap, ind(ix+1,iy), times, backpointers )
        if (1 < iy)  call pushheap( heap, ind(ix,iy-1), times, backpointers )
        if (iy < ny) call pushheap( heap, ind(ix,iy+1), times, backpointers )        
    
      ! propagate solution
      
        do while (nalive <= nx*ny)
        
          ! get the narrowband element with the smallest value for T
	      ! and make it alive
   
            call popheap( heap, imin, times, backpointers )
            if (imin == 0) exit
            ix = mod(imin-1,nx) + 1
            iy = (imin-1)/nx + 1
            backpointers(ix,iy) = ALIVE
            nalive = nalive + 1
            
          ! update neighbors
            if (1 < ix)  call update_neighbor(ix-1,iy)
            if (ix < nx) call update_neighbor(ix+1,iy)
            if (1 < iy)  call update_neighbor(ix,iy-1)
            if (iy < ny) call update_neighbor(ix,iy+1)
            
        end do
        
      contains
        
        subroutine update_neighbor(ix,iy)
        
            integer, intent(in) :: ix, iy
            real :: a,b,c,d,aa,cc
            real :: t, s, told
            integer :: i
            
            i = (iy-1)*nx+ix
            
            if (backpointers(ix,iy) .eq. ALIVE) return
            if (backpointers(ix,iy) .eq. FARAWAY) then
                call pushheap( heap, i, times, backpointers )
            end if
            
            a = infinity
            b = infinity
            c = infinity
            d = infinity
            
            told = times(ix,iy)
            
          ! values of t at neighbors
            if (1 < ix)  a = times(ix-1,iy)
            if (ix < nx) b = times(ix+1,iy)
            if (1 < iy)  c = times(ix,iy-1)
            if (iy < ny) d = times(ix,iy+1)
            
          ! pick maximum solution
            t = 0.
            aa = min(a,b)
            cc = min(c,d)
            if (max(aa,cc) /= infinity) then
                s = dx**2*dy**2*(dx**2+dy**2-((aa-cc)*speed(ix,iy))**2)
                if (s>=0.) then
                    t = max(t, ( (aa*dy**2+cc*dx**2)*speed(ix,iy)+sqrt(s) ) / &
                               ( speed(ix,iy)*(dx**2+dy**2) ) &
                           )
                end if
            end if
            if (min(c,d) == infinity) then
                if (a < infinity) t = max(t, a+dx/speed(ix,iy))
                if (b < infinity) t = max(t, b+dx/speed(ix,iy))
            end if
            if (min(a,b) == infinity) then
                if (c < infinity) t = max(t, c+dy/speed(ix,iy))
                if (d < infinity) t = max(t, d+dy/speed(ix,iy))
            end if
            
            if (t == 0.) then  ! fallback condition, should rarely happen.
                               ! may happen at sharp edges in speed function.
                t = infinity
                if (a < infinity) t = min(t, a+dx/speed(ix,iy))
                if (b < infinity) t = min(t, b+dx/speed(ix,iy))
                if (c < infinity) t = min(t, c+dy/speed(ix,iy))
                if (d < infinity) t = min(t, d+dy/speed(ix,iy))
            end if
            
            if (t /= 0. .and. told /= t) then
                call updateheap( heap, ind(ix,iy), t, times, backpointers )
            end if
        
        end subroutine
        
        function ind(ixl,iyl) 
            integer, intent(in) :: ixl, iyl
            integer :: ind
            
            ind = (iyl-1)*nx+ ixl 
        
        end function
        
    end subroutine

end module

