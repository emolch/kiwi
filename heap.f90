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

module heap

  ! This module implements an index heap with backpointer tracking.
  !
  ! The concepts used here are IT folklore and can be found in any text about 
  ! heapsort.
  !
  ! The only addition is, that sort keys may be changed, while simultaneously 
  ! maintaining the heap. (see updateheap())
  
    implicit none

    type, public :: t_index_heap
        integer, dimension(:), allocatable :: iheap
        integer :: n
    end type

    private
    
    public initheap
    public destroyheap
    
    public pushheap
    public popheap
    public updateheap
    public buildheap
    public downheap
    public upheap
    
  contains
  
    pure subroutine initheap( heap, maxsize )
    
        type(t_index_heap), intent(inout) :: heap
        integer, intent(in) :: maxsize
    
        call destroyheap( heap )
        
        allocate( heap%iheap(maxsize) )
        heap%n = 0
    
    end subroutine
    
    pure subroutine destroyheap( heap )
    
        type(t_index_heap), intent(inout) :: heap
        if (allocated( heap%iheap )) deallocate( heap%iheap )
        heap%n = 0
    
    end subroutine
    
  
    pure subroutine pushheap( heap, keyindex, keys, backpointers )
    
        type(t_index_heap), intent(inout) :: heap
        integer, intent(in)               :: keyindex
        real, dimension(*), intent(in)    :: keys ! lookup sort keys in this array
      ! maintain backward lookup of heapindices, if given
        integer, dimension(*), intent(inout), optional :: backpointers
        
        if (heap%n+1 > size(heap%iheap)) then
            ! capacity exceeded
            return
        end if
        
        heap%n = heap%n + 1
        heap%iheap(heap%n) = keyindex 
        
        if (present(backpointers)) then
            backpointers(keyindex) = heap%n
            call upheap( heap, heap%n, keys, backpointers )
        else 
            call upheap( heap, heap%n, keys )
        end if
    
    end subroutine
    
    pure subroutine popheap( heap, keyindex, keys, backpointers )
    
        type(t_index_heap), intent(inout) :: heap
        integer, intent(out)              :: keyindex
        real, dimension(*), intent(in)    :: keys ! lookup sort keys in this array
      ! maintain backward lookup of heapindices, if given
        integer, dimension(*), intent(inout), optional :: backpointers
        
        if (heap%n .eq. 0) then
            keyindex = 0
            return
        end if
        
        
        call swap( heap%iheap(1), heap%iheap(heap%n) )
        if (present(backpointers)) then
            call swap( backpointers(heap%iheap(1)), backpointers(heap%iheap(heap%n)) )
            backpointers(heap%iheap(heap%n)) = 0
        end if
        
        keyindex = heap%iheap(heap%n)
        heap%n = heap%n - 1
        if (present(backpointers)) then
            call downheap( heap, 1, keys, backpointers )
        else
            call downheap( heap, 1, keys )
        end if
        
    end subroutine
    
    pure subroutine updateheap( heap, keyindex, newkey, keys, backpointers )
    
      ! change a key, while keeping heap up to date
      ! (backpointers are required here to locate the element in the heap)
    
        type(t_index_heap), intent(inout)    :: heap
        integer, intent(in)                  :: keyindex
        real, intent(in)                     :: newkey
        real, dimension(*), intent(inout)    :: keys ! lookup sort keys in this array
        integer, dimension(*), intent(inout) :: backpointers
        
        real :: oldkey
        
        oldkey = keys(keyindex)
        keys(keyindex) = newkey
        
        if (newkey .lt. oldkey) then ! move it upward
            call upheap( heap, backpointers(keyindex), keys, backpointers )
        end if
        
        if (newkey .gt. oldkey) then ! move it downward
            call downheap( heap, backpointers(keyindex), keys, backpointers )
        end if 
    
    end subroutine
    
    pure subroutine buildheap( heap, keys, backpointers )
    
        type(t_index_heap), intent(inout) :: heap
        real, dimension(*), intent(in)    :: keys ! lookup sort keys in this array
      ! maintain backward lookup of heapindices, if given
        integer, dimension(*), intent(inout), optional :: backpointers
        
        integer :: v
        
        do v=heap%n/2,1,-1
            if (present(backpointers)) then
                call downheap(heap, v, keys, backpointers )
            else 
                call downheap(heap, v, keys )
            end if
        end do
        
    end subroutine
  
    pure subroutine downheap( heap, element, keys, backpointers )
        
        type(t_index_heap), intent(inout) :: heap
        integer, intent(in)               :: element ! heapindex of element to be inserted
        real, dimension(*), intent(in)    :: keys ! lookup sort keys in this array
      ! maintain backward lookup of heapindices, if given
        integer, dimension(*), intent(inout), optional :: backpointers
        
        integer :: v,w
        
        v = element
        w = 2*(v-1) + 2 ! first descendant of v
        do while (w .le. heap%n) 
        
          ! if w has second descendant and it is greater, take that one
            if ( w+1 .le. heap%n .and. &
                 keys(heap%iheap(w+1)) .lt. keys(heap%iheap(w)) ) w = w+1
          
          ! w points to largest descendant
            if ( keys(heap%iheap(v)) .le. keys(heap%iheap(w)) ) return ! v has heap property
            
            call swap( heap%iheap(v), heap%iheap(w) )
            if (present(backpointers)) then
                call swap( backpointers(heap%iheap(v)), backpointers(heap%iheap(w)) )
            end if
            
            v = w
            w = 2*(v-1) + 2
        end do
        
    end subroutine

    pure subroutine upheap( heap, element, keys, backpointers )
    
        type(t_index_heap), intent(inout) :: heap
        integer, intent(in)               :: element ! heapindex of element to be inserted
        real, dimension(*), intent(in)    :: keys ! lookup sort keys in this array
      ! maintain backward lookup of heapindices, if given
        integer, dimension(*), intent(inout), optional :: backpointers
        
        integer :: v,u
    
        v = element
        do while (v .gt. 1)
            u = (v-2)/2 + 1 ! parent
            if ( keys(heap%iheap(u)) .le. keys(heap%iheap(v)) ) return ! u has heap prop
            
            call swap( heap%iheap(u), heap%iheap(v) )
            if (present(backpointers)) then
                call swap( backpointers(heap%iheap(u)), backpointers(heap%iheap(v)) )
            end if
            
            v = u
        end do
        
    end subroutine
    
    pure subroutine swap( a,b )
        integer, intent(inout) :: a,b
        integer :: t
        t = a
        a = b
        b = t
    end subroutine
    
end module