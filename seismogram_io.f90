! $Id: seismogram_io.f90 680 2007-11-20 09:08:13Z sebastian $ 
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

module seismogram_io

  ! this module provides methods for writing seismograms
  ! at the moment ascii table output and sac output are intended,
  ! but other formats may be added

    
    use better_varying_string
    use read_table
    use unit
    use util
    
    implicit none
    
    private
    
    public :: writeseismogram, readseismogram
    
    interface writeseismogram
        module procedure writeseismogram_vs
        module procedure writeseismogram_c
    end interface
    
    interface readseismogram
        module procedure readseismogram_c
    end interface
    
  contains
  
    subroutine writeseismogram_vs( filename, fileformat, seismogram, toffset, deltat, nerr )
        
        type(varying_string), intent(in) :: filename, fileformat
        real, dimension(:), intent(in)   :: seismogram
        real, intent(in)                 :: toffset, deltat
        integer, intent(out)             :: nerr
        
        call writeseismogram_c( char(filename), char(fileformat), seismogram, toffset, deltat, nerr )
        
    end subroutine
  
    subroutine writeseismogram_c( filename, fileformat, seismogram, toffset, deltat, nerr )
    
        character(len=*), intent(in)     :: filename, fileformat
        real, dimension(:), intent(in)   :: seismogram
        real, intent(in)                 :: toffset, deltat
        integer, intent(out)             :: nerr

        type(varying_string)             :: fileformat_

      ! dump 1-component seismogram to file
        
        integer :: iunit, i, nlen
        real(kind=8) :: dtoffset, ddeltat
        real, dimension(size(seismogram)) :: seismogram_copy
        character(len=len(filename)+1) :: filename_cstr
             
      ! look at the filename extension what kind of output is wanted
        nerr = 0
        fileformat_ = trim(fileformat)
        if (fileformat == '*') then
            fileformat_ = 'table'
            if ( len(filename) >= 4 ) then
                if (filename(len(filename)-3:len(filename)) == '.sac') &
                    fileformat_ = 'sac'
            end if
            if ( len(filename) >= 6 ) then
                if (filename(len(filename)-5:len(filename)) == '.mseed') &
                    fileformat_ = 'mseed'
            end if
        end if
      
        if (fileformat_ == 'sac') then
            seismogram_copy(:) = seismogram(:)
            
            filename_cstr = filename//char(0)
            nlen = size(seismogram_copy)
            call wsac1( filename_cstr, seismogram_copy, nlen, &
                        toffset, deltat, nerr )
	   !  call die("sac output currently disabled, because there is no 64bit libsacio.")
        end if
        
        if (fileformat_ == 'mseed') then
            seismogram_copy(:) = seismogram(:)            
            
            dtoffset = toffset
            ddeltat = deltat
            filename_cstr = filename//char(0)
            nlen = size(seismogram_copy)
            call writemseed( filename_cstr, seismogram_copy, nlen, &
                dtoffset, ddeltat, nerr )
            
        end if        
        
        if (fileformat_ == 'table') then
        
            call claim_unit( iunit )
            open( unit=iunit, file=filename, status='unknown', iostat=nerr )
            
            if ( nerr /= 0 ) then
                call release_unit( iunit )
                return
            end if
            
            do i=1,size(seismogram)
                write (iunit,*) toffset+(i-1)*deltat, seismogram(i)
            end do
            
            close( iunit )
            call release_unit( iunit )
            
        end if
        
    end subroutine
    
    subroutine readseismogram_c( filename, fileformat, seismogram, toffset, deltat, nerr )
        
        character(len=*), intent(in)                  :: filename, fileformat
        real, dimension(:), allocatable, intent(inout)   :: seismogram
        real, intent(out) :: toffset, deltat
        integer, intent(out) :: nerr
        
        type(varying_string) :: vsfn
        integer :: nlen, maxlen
        real, dimension(:),allocatable :: scratch
        real, dimension(:,:),allocatable :: tablebuf
        type(varying_string) :: fileformat_
        character(len=len(filename)+1) :: filename_cstr
        real(kind=8) :: dtoffset, ddeltat
        
      ! this is not a very brilliant way to determine the file type, but
      ! anyway, look at the filename extension if this is sac data...
        fileformat_ = fileformat
        if (fileformat == '*') then
            fileformat_ = 'table'
            if ( len(filename) >= 4 ) then
                if (filename(len(filename)-3:len(filename)) == '.sac') &
                    fileformat_ = 'sac'
            end if
            if ( len(filename) >= 6 ) then
                if (filename(len(filename)-5:len(filename)) == '.mseed') &
                    fileformat_ = 'mseed'
            end if
        end if
        
        if (fileformat_ == 'sac') then
        
          ! read data into buffer 'scratch'...
        
          ! couldn't find a way to get the correct data length, without reading the
          ! full thing... so this stupid iterative attempt is used.
          !   - getnhv('NPTS',...) always returned maxlen when there was more data in file
          !   - nerr is not set to -803 when there is more data, opposing to what
          !     the manual tells
          !      (at least with my copy of libsacio) 
          ! what a crap!
            
            nerr = -1
            maxlen = 1024
            do while (nerr<0)
                if (allocated(scratch)) deallocate(scratch)
                allocate( scratch(maxlen) )
                
                filename_cstr = filename//char(0)
                call rsac1(filename_cstr, scratch, nlen, toffset, deltat, maxlen, nerr)
          !        call die("sac input currently disabled, because there is no 64bit libsacio")
                if (nerr > 0) then
                    call error( "rsac1 returned an error" )
                    if (allocated(scratch)) deallocate(scratch)
                    return
                end if 
                if (maxlen == nlen) nerr = -1

                maxlen = maxlen*2
            end do
            
            if (allocated(seismogram)) deallocate(seismogram)
            allocate( seismogram(nlen) )
            seismogram(:) = scratch(1:nlen)
            deallocate(scratch)
        end if
    
        if (fileformat_ == 'mseed') then
          
            filename_cstr = filename//char(0)
            
            call readmseed1( filename_cstr, nlen, dtoffset, ddeltat, nerr )
            if (nerr .ne. 0) then
                call error( "readmseed1 returned an error" )
                return
            end if
            
            if (allocated(seismogram)) deallocate(seismogram)
            allocate( seismogram(nlen) )
            
            call readmseed2( seismogram )
            
            toffset = real(dtoffset)
            deltat = real(ddeltat)
            
        end if
        
        if (fileformat_ == 'table') then
            vsfn = filename
            call readtable_file( tablebuf, vsfn, min_cols=2, min_rows=2, nerr=nerr )
            if (nerr /= 0) then
                return
            end if 
            nlen = size(tablebuf,2)
            if (allocated(seismogram)) deallocate(seismogram)
            allocate( seismogram(nlen) )
            seismogram(:) = tablebuf(2,:)
            toffset = tablebuf(1,1)
            deltat = (tablebuf(1,nlen)-tablebuf(1,1))/(nlen-1)
            deallocate( tablebuf )
            
        end if
        
    end subroutine
    

end module
