! $Id: test_io.f90 668 2007-09-11 15:09:11Z sebastian $ 
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


program test_io


    use util
    use sparse_trace
    use comparator
    use seismogram_io
    
    implicit none

    type(t_strip) :: strip1
    integer :: nerr
    real(kind=8) :: toffset
    real :: deltat
    real, allocatable, dimension(:) :: seismogram
    
    call test_begin("test_io")
    call strip_init( (/1,5/), (/0.,0.,1.,0.,0./), strip1 )
    deltat = 0.5
    toffset = 0.
    call writeseismogram("/tmp/test_io.table", "*", strip1%data, toffset, deltat, '','','','', nerr )
    if (nerr /= 0) then
        call test_fail("write 1")
    end if
    
    call readseismogram("/tmp/test_io.table", "*", seismogram, toffset, deltat, nerr )
    if (nerr /= 0) then
        call test_fail("read 1")
    end if
     if (toffset /= 0. .or. deltat /= 0.5) then
        call test_fail("read 1a")
    end if
    
    if (any(seismogram .ne. strip1%data)) then
        call test_fail("read 1b")
    end if 
    
    
    if (toffset /= 0. .or. deltat /= 0.5) then
        call test_fail("read 1")
    end if
    
    if (any(seismogram .ne. strip1%data)) then
        call test_fail("read 1")
    end if 
    
    
    call strip_init( (/1,5/), (/0.,0.,1.,2.,0./), strip1 )
    
!     call writeseismogram("test.sac", "*", strip1%data, toffset, deltat, nerr )
!     if (nerr /= 0) then
!         call test_fail("write 2")
!     end if
!     
!     call readseismogram("test.sac", "*", seismogram, toffset, deltat, nerr )
!     if (nerr /= 0) then
!         call test_fail("read  2")
!     end if
!     
!     if (toffset /= 0. .or. deltat /= 0.5) then
!         call test_fail("read 2a")
!     end if
!     
!     if (any(seismogram .ne. strip1%data)) then
!         call test_fail("read 2b")
!     end if 

    call writeseismogram("/tmp/test_io.mseed", "*", strip1%data, toffset, deltat, '','','','', nerr )
    if (nerr /= 0) then
        call test_fail("write 3")
    end if
    
    call readseismogram("/tmp/test_io.mseed", "*", seismogram, toffset, deltat, nerr )
    if (nerr /= 0) then
        call test_fail("read  3")
    end if
    
    if (toffset /= 0. .or. deltat /= 0.5) then
        call test_fail("read 3a")
    end if
    
    if (any(seismogram .ne. strip1%data)) then
        call test_fail("read 3b")
    end if 
    
    
    call test_end()
    call cleanup()
    
end program
