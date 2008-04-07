! $Id: test_gfdb.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


program test_gfdb

    use util
    use gfdb
    use sparse_trace
    use better_varying_string
    
    implicit none
    
    type(t_gfdb) :: db
    type(t_trace) :: tr1, tr2, tr
    type(t_strip) :: st1, st2
    type(t_trace), pointer ::trp => null()
    type(t_trace), pointer ::trp2 => null()
    
    call test_begin("test_gfdb")
    
    call gfdb_init( db, var_str("/tmp/test-db"), 2, 3, 2, 8, 1., 1., 1. )
    
  ! reopen and look if contents are set corrently
    call gfdb_init( db, var_str("/tmp/test-db") )
    if (any ( (/db%nchunks,db%nx,db%nz,db%ng/) /= (/2,3,2,8/) ))&
        call test_fail("reopen index")
    
  ! save a trace
    call strip_init( (/11,17/), (/1.,0.,0.,0.,1.,1.,1./), st1 )
    call strip_init( (/23,24/), (/1.,1./), st2 )
    call trace_pack( st1, tr1 )
    call trace_pack( st2, tr2 )
    call trace_join( tr1, tr2, tr )
    
    call gfdb_save_trace( db, 1, 2, 2, tr )
    
  ! recover it
    call gfdb_get_trace( db, 1, 2, 2, trp )
    
    if (.not. associated(trp)) then
        call test_fail("recover 0")
    else
        if (trace_is_empty(trp)) then
            call test_fail("recover 1")
        else 
    
            if (trp%nstrips /= 2) then
                call test_fail("recover 1")
            else 
                if (any(strip_span(trp%strips(1)) /= (/11,17/))) &
                    call test_fail("recover 2")
                if (any(strip_span(trp%strips(2)) /= (/23,24/))) &
                    call test_fail("recover 3")
                if (any(trp%strips(1)%data /= st1%data)) &
                        call test_fail("recover 4")
                if (any(trp%strips(2)%data /= st2%data)) &
                        call test_fail("recover 5")
            end if 
        end if
    end if
    
  ! recover it again
    call gfdb_get_trace( db, 1, 2, 2, trp2 )
    if (.not. associated(trp,trp2)) &
       call test_fail("cached get")
   
    
  ! try to recover non-existent trace
  !  call gfdb_get_trace( db, 2, 2, 2, trp )
    
  !  if (.not. associated(trp)) then
  !      call test_fail("recover nonexistent")
  !  else 
  !      if (.not. trace_is_empty(trp)) then
  !          call test_fail("recover nonexistent")
  !      end if
  !  end if
    
  ! out of bounds must return null
  !  call gfdb_get_trace( db, 9,9,9, trp )
  !  if ( associated(trp) ) &
  !      call test_fail("recover bad")
  ! ... worked, but commented out, because it emmits a warning
    
    call gfdb_destroy( db )
  
    call test_end()
    call cleanup()
        
end program
