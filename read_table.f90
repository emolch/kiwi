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

module read_table
    
    use util
    use unit
    use read_line
    use better_varying_string
    
    implicit none

    private    
  ! type defs
  
    integer, parameter, private :: CHUNK_ROWS = 256
    logical :: myok
    
    type, private :: chunk
        private
        real, dimension(:,:), allocatable :: content
        integer                           :: nrows
        type(chunk), pointer              :: next
    end type chunk

  ! globals
  
    real, allocatable, dimension(:), private :: a
    integer, private                         :: n_detected_cols
    logical, private                         :: need_col_detect
    
        
  ! access specs
  
    public  :: readtable, readtable_file 
    private :: line_found, count_words, whitespace, initchunk
        
   contains
     
     subroutine readtable_file( field, fn, min_cols, min_rows, nerr  )
        real, dimension(:,:), allocatable, intent(out) :: field
        type(varying_string), intent(in) :: fn
        integer, intent(in), optional :: min_cols, min_rows
        integer, intent(out), optional :: nerr
       
    ! this reads a table of real numbers from file filename to field, which is
    ! allocated to the proper size.
    ! it is checked, that the table has at least min_cols columns and min_rows rows
    ! if these args are present
    ! if nerr is present, nerr is set to 0 if all is ok or somthing /= 0 if not
    ! g_errstr is set 
        
        integer :: ifile, iostat
        
        if ( present(nerr) ) nerr = 0
        
        call claim_unit( ifile )
        open( unit=ifile, file=char(fn), status='old', iostat=iostat )
        if (iostat /= 0) then
            call error('can''t open file "' // fn // '"')
            if ( present(nerr) ) then
                nerr = 1
                return
            else 
                call die()
            end if
        end if
        
        call readtable( field, iunit=ifile )
        close( ifile )
        call release_unit( ifile )
        
        if ( present(min_cols) ) then
            if (size(field,1) <  min_cols) then
                call error( 'expected at least ' // min_cols // ' column(s) in file "' // fn // '"' )
                if ( present(nerr) ) then
                    nerr = 1
                    return
                else 
                    call die()
                end if
            end if
        end if
        
        if ( present(min_rows) ) then
            if (size(field,1) <  min_rows) then
                call error( 'expected at least ' // min_rows // ' row(s) in file "' // fn // '"' )
                if ( present(nerr) ) then
                    nerr = 1
                    return
                else 
                    call die()
                end if
            end if
        end if
        
    end subroutine

    subroutine readtable( field, n_wanted_cols, iunit )
        real, allocatable, dimension(:,:), intent(inout) :: field
        integer, intent(in), optional                    :: n_wanted_cols
        integer, intent(in), optional                    :: iunit
     
      ! read ascii table into 2D-Array field
      ! field is allocated/reallocated if needed
      ! the number of wanted columns may be specified in n_wanted_cols
      ! the data is read from stdin or iunit if present 
     
        type(chunk), pointer  :: first, current, new
        
        integer               :: iostat
        logical               :: ok
        
        integer               :: iline, ncols, nchunks, nrows, i, ioffset
        
        ncols = 0
        need_col_detect = .false.
        if ( present( n_wanted_cols ) ) then
            ncols = n_wanted_cols
            allocate ( a(ncols) )
        else
            need_col_detect = .true.
        end if
        
        first => null()

        iline = 0
        nchunks = 0
        line_loop : do
            myok = .false.
            if ( present(iunit) ) then
                call readline( line_found, iostat, ok, iunit)
            else
                call readline( line_found, iostat, ok )
            end if
            if (ok .and. myok) then
                iline = iline + 1
                if (iline == 1) then
                    if (ncols == 0) ncols = n_detected_cols
                    allocate( first )
                    call initchunk( first, ncols )
                    current => first
                    nchunks = 1
                end if
                if (current%nrows == CHUNK_ROWS) then
                    allocate( new )
                    call initchunk( new, ncols )
                    current%next => new
                    current => new
                    nchunks = nchunks + 1
                end if
                current%nrows = current%nrows + 1 
                current%content(1:ncols, current%nrows) = a(1:ncols)
            end if
            
            if (iostat == IOSTAT_EOF) exit line_loop
        end do line_loop
        if ( allocated(a) ) deallocate( a )
        
        if ( nchunks == 0 ) then
            if (allocated(field)) then
                deallocate( field )
            end if
            return
        end if
        nrows = (nchunks-1)*CHUNK_ROWS + current%nrows
        
        if ( allocated( field ) ) then
            if ( ncols /= size(field,1) .or. nrows /= size(field,2) ) then
                deallocate( field )
            end if
        end if
        if ( .not. allocated( field ) )  allocate( field( ncols, nrows ) )
        
        current => first
        do i=1,nchunks
            ioffset = (i-1)*CHUNK_ROWS
            if (current%nrows >= 1) then
                field( 1:ncols, ioffset+1:ioffset+current%nrows ) = current%content( 1:ncols, 1:current%nrows )
            end if
            new => current%next
            deallocate(current%content)
            deallocate(current)
            current => new
        end do
        
    end subroutine readtable
    
    subroutine initchunk( newchunk, ncols)
        type(chunk), pointer :: newchunk
        integer, intent(in)  :: ncols
        
        allocate( newchunk%content(ncols, CHUNK_ROWS) )
        newchunk%nrows = 0
        
    end subroutine initchunk

    
    subroutine line_found( buffer, ok )
    
        character(len=*), intent(in) :: buffer
        logical, intent(out)         :: ok
        
        ! called by readline for every non-comment-line it finds
        
        integer :: iostat
        
        if ( need_col_detect ) then
            n_detected_cols = count_words( buffer )
            if ( .not. allocated( a ) ) allocate( a(n_detected_cols) )
            need_col_detect = .false.
        end if
        ok = .false.
        read (unit=buffer,fmt=*,iostat=iostat) a
        if (iostat == 0) then
            ok = .true.
            myok = .true.
        end if
    end subroutine line_found
    
end module read_table
