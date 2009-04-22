! $Id: sparse_trace.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


module sparse_trace

    use util

    implicit none

    private
    
    integer, parameter :: maxgap = 5

    integer, parameter :: real_kind = 4

    type, public :: t_strip  
      ! continuous strip of data at a specific offset
        real(kind=real_kind), dimension(:), allocatable :: data
    
    end type
    
    type, public :: t_trace
    
      ! this contains a single elementary seismogram trace efficiently,
      ! so that it may be broken up into several time windows, of where the 
      ! function is non-zero.
      ! the end point is repeated.
      
      ! functions are written, such that reading of traces is fast
      ! while writing or modifing traces may me slow
      
        integer :: nstrips
        integer, dimension(2) :: span
        
        type(t_strip), dimension(:), allocatable :: strips
        
    end type
    
    interface strip_copy
        module procedure strip_copy_s
        module procedure strip_copy_a
    end interface 
    
    public strip_init, strip_span, strip_copy, strip_extend, strip_destroy
    public strip_intersection, strip_length, strip_dataspan
    public strip_extend_to_same_span_5
    public strip_extend_to_same_span_4
    public strip_extend_to_same_span_2
    public strip_fold
    public strip_nullify
    
    public trace_join, trace_pack, trace_unpack, trace_multiply_add, trace_multiply_add_nogrow
    public trace_copy, trace_destroy, trace_is_empty
    public trace_from_storable, trace_to_storable
    public trace_create_simple, trace_create_simple_nodata
    public trace_size_bytes

  contains
  
    subroutine strip_init( span, data, strip )
    
        integer, dimension(:), intent(in) :: span  ! (2)
        real, dimension(:), intent(in) :: data
        type(t_strip), intent(inout) :: strip
        
        integer :: length
        
        length = span(2)-span(1) + 1
        if (size(data) /= length) then
            stop 'strip_init(): length of data does not match span'
        end if
        
        call resize( strip%data, span(1), length )
        strip%data = data
        
    end subroutine
    
    pure function strip_length( strip ) result(length)
        type(t_strip), intent(in) :: strip
        integer :: length
        
        if (.not. allocated(strip%data)) then
            length=0
            return
        end if
        length=size(strip%data)
    end function
    
    subroutine strip_destroy( strip )
        type(t_strip), intent(inout) :: strip
        call resize( strip%data, 1, 0 )
    end subroutine

    pure subroutine strip_nullify( strip )
        type(t_strip), intent(inout) :: strip
        if (allocated(strip%data)) then
            strip%data(:) = 0.0
        end if
    end subroutine
    
    pure function strip_span( strip )
        type(t_strip), intent(in) :: strip
        integer, dimension(2) :: strip_span
        strip_span(1) = lbound(strip%data,1)
        strip_span(2) = ubound(strip%data,1)
    end function
    
    pure subroutine strip_intersection( s1, s2, span )
    
        type(t_strip), intent(in) :: s1, s2
        integer, dimension(2), intent(out) :: span
        integer, dimension(2)  :: span1,span2
        
        span1 = strip_span( s1 )
        span2 = strip_span( s2 )
        
        call intersection(span1, span2, span)
    
    end subroutine

    subroutine trace_join( a, b, c )
      
      ! join two sparse traces a and b inefficiently
      ! by copying everything to c
      
        type(t_trace), intent(in) :: a,b
        type(t_trace), intent(inout) :: c
        
        call trace_destroy( c )
        
        if (.not. allocated(a%strips)) then 
            call trace_copy( b, c )
            return
        end if
        
        if (.not. allocated(b%strips)) then
            call trace_copy( a, c )
            return
        end if
        
        if (a%span(2) >= b%span(1)) then
            stop 'trace_join(): span overlap detected'
        end if
        
        c%nstrips = a%nstrips + b%nstrips
        c%span(1) = a%span(1)
        c%span(2) = b%span(2)
        
        allocate( c%strips( c%nstrips ) )
        
        call strip_copy( a%strips(:), c%strips(1:a%nstrips) )        
        call strip_copy( b%strips(1:b%nstrips), &
                         c%strips(a%nstrips+1:a%nstrips+b%nstrips) )
        
    end subroutine
    
    subroutine trace_copy( source, dest )
    
        type(t_trace), intent(in)    :: source
        type(t_trace), intent(inout) :: dest    
    
        call trace_destroy( dest )
        if (allocated(source%strips)) then
            allocate(dest%strips(source%nstrips)) 
            dest%span = source%span
            dest%nstrips = source%nstrips
            call strip_copy(source%strips,dest%strips)
        end if
    
    end subroutine
    
    
    subroutine strip_copy_a( source, dest )

      ! copy strips source to dest
      
        type(t_strip), dimension(:), intent(in)  :: source
        type(t_strip), dimension(:), intent(inout) :: dest
        
        integer :: i

        if (size(source) /= size(dest)) then
            stop 'strip_copy(): sizes of source and dest do not match'
        end if
        
        do i=1,size( source )
            call resize( dest(i)%data, lbound(source(i)%data,1) , size(source(i)%data) )
            dest(i)%data = source(i)%data
        end do
        
    end subroutine
    
    subroutine strip_copy_s( source, dest )

      ! copy strips source to dest
      
        type(t_strip), intent(in)  :: source
        type(t_strip), intent(inout) :: dest
        
        call resize( dest%data, lbound(source%data,1) , size(source%data) )
        dest%data = source%data
        
    end subroutine
         
    subroutine strip_extend_to_same_span_5( a, b, c, d, e )
    
        type(t_strip), intent(inout)  :: a,b,c,d,e
        
        integer, dimension(2) :: span
        
        span(1) = huge(span(1))
        span(2) = -huge(span(2))
        if (allocated(a%data)) then
            span(1) = min(lbound(a%data,1), span(1))
            span(2) = max(ubound(a%data,1), span(2))
        end if
        if (allocated(b%data)) then
            span(1) = min(lbound(b%data,1), span(1))
            span(2) = max(ubound(b%data,1), span(2))
        end if
        if (allocated(c%data)) then
            span(1) = min(lbound(c%data,1), span(1))
            span(2) = max(ubound(c%data,1), span(2))
        end if
        if (allocated(d%data)) then
            span(1) = min(lbound(d%data,1), span(1))
            span(2) = max(ubound(d%data,1), span(2))
        end if
        if (allocated(e%data)) then
            span(1) = min(lbound(e%data,1), span(1))
            span(2) = max(ubound(e%data,1), span(2))
        end if
        
        if (span(1) < span(2)) then
            call strip_extend( a, span )
            call strip_extend( b, span )
            call strip_extend( c, span )
            call strip_extend( d, span )
            call strip_extend( e, span )
        end if
    
    end subroutine
    
     subroutine strip_extend_to_same_span_4( a, b, c, d )
    
        type(t_strip), intent(inout)  :: a,b,c,d
        
        integer, dimension(2) :: span
        
        span(1) = huge(span(1))
        span(2) = -huge(span(2))
        if (allocated(a%data)) then
            span(1) = min(lbound(a%data,1), span(1))
            span(2) = max(ubound(a%data,1), span(2))
        end if
        if (allocated(b%data)) then
            span(1) = min(lbound(b%data,1), span(1))
            span(2) = max(ubound(b%data,1), span(2))
        end if
        if (allocated(c%data)) then
            span(1) = min(lbound(c%data,1), span(1))
            span(2) = max(ubound(c%data,1), span(2))
        end if
        if (allocated(d%data)) then
            span(1) = min(lbound(d%data,1), span(1))
            span(2) = max(ubound(d%data,1), span(2))
        end if
        
        if (span(1) < span(2)) then
            call strip_extend( a, span )
            call strip_extend( b, span )
            call strip_extend( c, span )
            call strip_extend( d, span )
        end if
    
    end subroutine
    
    pure subroutine strip_extend_to_same_span_2( a, b )
    
        type(t_strip), intent(inout)  :: a,b
        
        integer, dimension(2) :: span
        
        span(1) = huge(span(1))
        span(2) = -huge(span(2))
        
        if (allocated(a%data)) then
            span(1) = min(lbound(a%data,1), span(1))
            span(2) = max(ubound(a%data,1), span(2))
        end if
        if (allocated(b%data)) then
            span(1) = min(lbound(b%data,1), span(1))
            span(2) = max(ubound(b%data,1), span(2))
        end if
        
        if (span(1) < span(2)) then
            call strip_extend( a, span )
            call strip_extend( b, span )
        end if
        
    end subroutine
    
    pure subroutine strip_extend( s, newspan )
    
      ! extend span of strip s to newspan, preserving its data  
      ! WARNING: this fails if newspan is smaller then the current span of s
        
        type(t_strip), intent(inout) :: s
        integer, dimension(2), intent(in) :: newspan
        
        real, dimension(:), allocatable :: temp
        integer, dimension(2) :: r
        
        if (allocated(s%data)) then
            r(1) = lbound(s%data,1)
            r(2) = ubound(s%data,1)
            allocate( temp(r(1):r(2)) )
            temp = s%data
        end if
        
        call resize( s%data, newspan(1),  newspan(2) - newspan(1) + 1 )
        
        if (allocated(temp)) then
            if (newspan(1) < r(1)) s%data(newspan(1):r(1)-1) = 0.
            if (newspan(2) > r(2)) s%data(r(2)+1:newspan(2)) = temp(r(2))
            s%data(r(1):r(2)) = temp(:)
            deallocate(temp)
        else
            s%data = 0.
        end if
        
    end subroutine
    
    function strip_dataspan( s )

        type(t_strip), intent(in) :: s
        integer, dimension(2) :: strip_dataspan

        integer, dimension(2) :: span
        real :: firstvalue, lastvalue
        integer :: i

        if (strip_length(s) == 0) then
            strip_dataspan = (/0,-1/)
            return
        end if

        span = strip_span(s)
        strip_dataspan = span

        firstvalue = 0.
        
        do i=span(1),span(2)
            strip_dataspan(1) = i
            if (s%data(i) /= firstvalue) exit 
        end do
        
        lastvalue = s%data(span(2))
        do i=span(2),span(1),-1
            if (s%data(i) /= lastvalue) exit 
            strip_dataspan(2) = i
        end do

    end function

    subroutine strip_fold( s, shifts, amplitudes )

        type(t_strip), intent(inout) :: s
        real, dimension(:), intent(in) :: shifts, amplitudes

        integer, dimension(2) :: dataspan
        type(t_trace) :: t
        integer :: i

        dataspan = strip_dataspan( s )
       
        if (dataspan(2) < dataspan(1)) then
            return
        end if

        call trace_create_simple( t, s%data(dataspan(1):dataspan(2)), dataspan )
        s%data(:) = 0.
        do i=1,size(shifts)
            call trace_multiply_add( t, s, factor_=amplitudes(i), rtraceshift_=shifts(i) )
        end do

        call trace_destroy(t)

    end subroutine

    pure subroutine trace_create_simple( trace, data, span )
    
        type(t_trace), intent(inout)     :: trace
        real, dimension(:), intent(in)   :: data
        integer, dimension(2),intent(in) :: span
  
        call trace_destroy( trace )
        
        allocate( trace%strips(1) )
        trace%nstrips = 1
        call resize( trace%strips(1)%data, span(1), span(2)-span(1)+1 )
        trace%strips(1)%data(:) = data(:)
        trace%span(:) = span(:)
        
    end subroutine
    
    pure subroutine trace_create_simple_nodata( trace, span )
    
        type(t_trace), intent(inout)     :: trace
        integer, dimension(2),intent(in) :: span
  
        call trace_destroy( trace )
        
        allocate( trace%strips(1) )
        trace%nstrips = 1
        call resize( trace%strips(1)%data, span(1), span(2)-span(1)+1 )
        trace%span(:) = span(:)
        
    end subroutine
    
    pure function span_contains( span, subspan )

        integer, dimension(2), intent(in) :: span, subspan
        logical :: span_contains

        span_contains =  (subspan(1) <= subspan(2) .and. span(1) <= subspan(1) .and. subspan(2) <= span(2))
        
    end function

    pure subroutine trace_pack( strip, trace, last_span_as_hint )
        
      ! this continuous trace into a sparse one
         
        type(t_strip), intent(in)     :: strip
        type(t_trace), intent(inout)    :: trace

      ! turns on hinting, of positions, where an empty trace is put:
        integer, dimension(2), intent(in), optional :: last_span_as_hint

        integer :: i,istrip,nstrips
        integer :: ibeg,iend
        logical :: interest
        integer :: gap
        integer, dimension(2) :: span
      
      ! this could be done more elegantly, utilizing
      ! a kind of growable data structure,
      ! but for now, as speed does not matter here,
      ! we do it as a two-pass
          
      ! first pass: count number of strips needed
        gap = 0
        interest = .false.
        istrip = 0
        do i=lbound(strip%data,1),ubound(strip%data,1)
            
            if (strip%data(i) /= 0.) then
                if (.not. interest) then
                    interest = .true.
                    istrip = istrip + 1
                end if
                gap = 0
            else
                if (interest) then
                    gap = gap + 1
                    if (gap > maxgap) then
                        interest = .false.
                    end if
                end if
            end if
            
        end do
        nstrips = istrip
        
        call trace_destroy( trace )

      ! special treatment of traces without data...
      ! save a single datapoint at a position which is not likely to
      ! stretch output seismograms too much.
        if ( nstrips == 0 ) then
            allocate( trace%strips(1) )
            span = strip_span( strip )
            ibeg = span(1)
            iend = ibeg
            if (present(last_span_as_hint)) then
                if (span_contains( span, last_span_as_hint )) then
                    ibeg = last_span_as_hint(1)
                    iend = ibeg
                end if
            end if
            call resize( trace%strips(1)%data, ibeg, iend-ibeg+1 ) ! add one of the zeros
            trace%strips(1)%data(:) = 0.0
            trace%nstrips = 1
            trace%span(1) = ibeg
            trace%span(2) = iend

            return
        end if

      ! adjust size of the sparse trace
        allocate( trace%strips(nstrips) )
        trace%nstrips = nstrips
      
      ! second pass: copy interesting data to the strips
        gap = 0
        interest = .false.
        istrip = 0
        do i=lbound(strip%data,1),ubound(strip%data,1)
            
            if (strip%data(i) /= 0.) then
                if (.not. interest) then
                    interest = .true.
                    ibeg = i
                    istrip = istrip + 1
                end if
                gap = 0
                iend = i
            else
                if (interest) then
                    gap = gap + 1
                    if (gap > maxgap) then
                        call resize( trace%strips(istrip)%data, ibeg, iend-ibeg+2 ) ! add one of the zeros
                        trace%strips(istrip)%data(:) = strip%data(ibeg:iend+1)
                        interest = .false.
                    end if
                end if
            end if
        end do
        if (interest) then
            if (gap > 0) then ! if there are zeros at the end, add one of them
                call resize( trace%strips(istrip)%data, ibeg, iend-ibeg+2 )
                trace%strips(istrip)%data(:) = strip%data(ibeg:iend+1)
            else
                call resize( trace%strips(istrip)%data, ibeg, iend-ibeg+1 )
                trace%strips(istrip)%data(:) = strip%data(ibeg:iend)
            end if    
        end if
        
        trace%span(1) = lbound(trace%strips(1)%data,1)
        trace%span(2) = ubound(trace%strips(nstrips)%data,1)
        
    end subroutine
    
    pure subroutine trace_unpack( trace, strip )
        
      ! this turns a sparse trace into a continuous one
      
        type(t_trace), intent(in)    :: trace
        type(t_strip), intent(inout)   :: strip
        
        integer :: istrip
        integer :: length
        integer, dimension(2) :: r
        
        length = trace%span(2)-trace%span(1)+1
        
        call resize( strip%data, trace%span(1), length )
        
        strip%data = 0.
        
        do istrip=1,trace%nstrips
            r(1) = lbound(trace%strips(istrip)%data,1)
            r(2) = ubound(trace%strips(istrip)%data,1)
            strip%data(r(1):r(2)) = trace%strips(istrip)%data(:)
        end do
        
    end subroutine
                                         
    pure subroutine trace_multiply_add( trace, strip, factor_, &
                                         itraceshift_, rtraceshift_ )
    
      ! This adds a sparse trace 'trace' multiplied by 'factor' on the continuous 'strip'.
      
      ! The output strip is extended to the appropriate range if neccessary.
      
      ! If 'itraceshift_' is present, the trace is shifted by that number of samples
      ! when it is added to strip.   strip(x) += trace(x-itraceshift_)
    
      ! If 'rtraceshift_' is present, the trace is shifted by a subsample distance,
      ! and linear interpolation is used to resample the data to the positions of
      ! strip. 
      
      ! rtraceshift and itraceshift should not be specified both.
      ! Just in case: rtraceshift takes precedence.
      
        type(t_trace), intent(in) :: trace
        type(t_strip), intent(inout) :: strip
        real, intent(in), optional :: factor_
        integer, intent(in), optional :: itraceshift_
        real, intent(in), optional :: rtraceshift_
        
        integer :: istrip
        integer, dimension(2) :: stripspan, span, r, contspan, itracespan_shift
        integer, dimension(2) :: outspan_needed
        integer :: itraceshift, outlength_needed
        real :: factor, weight_right, weight_left, lastval
        
        if (present(factor_)) then
            factor = factor_
        else
            factor = 1.
        end if
        
        
        itraceshift = 0
        
        if (present(itraceshift_)) then
            itraceshift = itraceshift_
        end if
        
        if (present(rtraceshift_)) then
            itraceshift = floor( rtraceshift_ )
          ! factors for the interpolation
            weight_right = rtraceshift_-itraceshift
            weight_left = 1.-weight_right
            weight_right = weight_right * factor
            weight_left = weight_left * factor
        end if
        
        itracespan_shift = trace%span + itraceshift
      
      ! this is the operation span
        span = itracespan_shift
        
        outspan_needed = span
      ! one sample more is needed if interpolation is wanted...
        if (present(rtraceshift_)) outspan_needed(2) = outspan_needed(2)+1
        
      ! extend output strip if it is too short
        if (allocated( strip%data)) then
            contspan(1) = min(outspan_needed(1),lbound(strip%data,1))
            contspan(2) = max(outspan_needed(2),ubound(strip%data,1))
            if (contspan(1) /= lbound(strip%data,1) .or. contspan(2) /= ubound(strip%data,1)) then
                call strip_extend( strip, contspan )
            end if
        else
            outlength_needed = outspan_needed(2)-outspan_needed(1)+1
            call resize( strip%data, outspan_needed(1), outlength_needed )
            strip%data = 0.
        end if
      
      ! multiply-add the intersecting strips
        do istrip=1,trace%nstrips
            stripspan(1) = lbound(trace%strips(istrip)%data,1)+itraceshift
            stripspan(2) = ubound(trace%strips(istrip)%data,1)+itraceshift
            
            if (stripspan(2) < span(1)) cycle
            if (stripspan(1) > span(2)) exit
            
            call intersection( stripspan, span, r )
            
            if (.not. present(rtraceshift_) ) then
               strip%data(r(1):r(2)) = strip%data(r(1):r(2)) &
                      + factor * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift)
            else  ! linear interpolation
                strip%data(r(1):r(2)) = strip%data(r(1):r(2)) &
                      + weight_left * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift)
                      
                if (istrip .eq. trace%nstrips) then !last point covered by repeat e.p. below
                    strip%data(r(1)+1:r(2)) = strip%data(r(1)+1:r(2)) &
                      + weight_right * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift-1)
                else 
                    strip%data(r(1)+1:r(2)+1) = strip%data(r(1)+1:r(2)+1) &
                      + weight_right * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift)
                end if                
            end if
            
          ! repeat end point if neccessary
           
            if ((istrip .eq. trace%nstrips) .and. (r(2)+1 .le. ubound(strip%data,1 ))) then
                lastval = trace%strips(istrip)%data(ubound(trace%strips(istrip)%data,1))
                if (lastval /= 0.) then
                    strip%data(r(2)+1:) = strip%data(r(2)+1:) + factor * lastval
                end if
            end if
            
        end do
        
    end subroutine

    
    pure subroutine trace_multiply_add_nogrow( trace, array, arrayspan, factor_, &
                                         itraceshift_, rtraceshift_ )
     
        type(t_trace), intent(in) :: trace
        integer, dimension(:), intent(in) :: arrayspan ! (2)
        real, dimension(arrayspan(1):), intent(inout) :: array
        real, intent(in), optional :: factor_
        integer, intent(in), optional :: itraceshift_
        real, intent(in), optional :: rtraceshift_
        
      ! same as trace_multiply_add but this version works on an array of fixed size.
      
        integer               :: istrip
        integer, dimension(2) :: stripspan, span, r, itracespan_shift
        integer               :: itraceshift
        real                  :: factor, weight_right, weight_left, lastval
        
        if (present(factor_)) then
            factor = factor_
        else
            factor = 1.
        end if
        
        itraceshift = 0
        
        if (present(itraceshift_)) then
            itraceshift = itraceshift_
        end if
        
        if (present(rtraceshift_)) then
            itraceshift = floor( rtraceshift_ )
          ! factors for the interpolation
            weight_right = rtraceshift_-itraceshift
            weight_left = 1.-weight_right
            weight_right = weight_right * factor
            weight_left = weight_left * factor
        end if
        
        itracespan_shift = trace%span + itraceshift
        
      ! determine operation span
        call intersection( arrayspan, itracespan_shift, span )
        if (span(2) < span(1)) return
        
      ! multiply-add the intersecting strips
        do istrip=1,trace%nstrips
            stripspan(1) = lbound(trace%strips(istrip)%data,1)+itraceshift
            stripspan(2) = ubound(trace%strips(istrip)%data,1)+itraceshift
            
            if (stripspan(2) < span(1)) cycle
            if (stripspan(1) > span(2)) exit
            
            call intersection( stripspan, span, r )
            
            if (.not. present(rtraceshift_) ) then
                array(r(1):r(2)) = array(r(1):r(2)) &
                  + factor * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift)
                  
            else  ! linear interpolation
                array(r(1):r(2)) = array(r(1):r(2)) &
                      + weight_left * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift)
                      
                if (istrip .eq. trace%nstrips .or. r(2)+1 > arrayspan(2) ) then   ! last point covered by repeat e.p. below
                    array(r(1)+1:r(2)) = array(r(1)+1:r(2)) &
                      + weight_right * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift-1)
                else
                    array(r(1)+1:r(2)+1) = array(r(1)+1:r(2)+1) &
                      + weight_right * trace%strips(istrip)%data(r(1)-itraceshift:r(2)-itraceshift)
                end if
            end if
            
         ! repeat end point if neccessary
           
            if ((istrip .eq. trace%nstrips) .and. (r(2)+1 .le. arrayspan(2))) then
                lastval = trace%strips(istrip)%data(ubound(trace%strips(istrip)%data,1))
                if (lastval /= 0.) then
                    array(r(2)+1:) = array(r(2)+1:) + factor * lastval
                end if
            end if
            
        end do
    
    end subroutine

    
    subroutine trace_to_storable( trace, packed, npacked, poffsets, offsets, noffsets )
    
      ! pack trace to a single, continuous array (packed)
      ! generate arrays with offsets in the packed array (poffsets)
      ! and offsets of the data strips (offsets)
      
      ! the output arrays are grown to fit the neccessary length,
      ! but never shortened.
      
      ! length of packed data is put to npacked
      ! number of offsets is put to noffsets
      
        type(t_trace), intent(in) :: trace
        real, dimension(:), allocatable, intent(inout) :: packed
        integer, dimension(:), allocatable, intent(inout) :: poffsets, offsets
        integer, intent(out) :: npacked, noffsets
        
        integer :: i, istrip
    
        noffsets = trace%nstrips
        
        if (.not. allocated(poffsets)) &
            call resize(poffsets,1,noffsets)
        if (noffsets > size(poffsets)) &
            call resize(poffsets,1,noffsets)
            
        if (.not. allocated(offsets)) &
            call resize(offsets,1,noffsets)
        if (noffsets > size(offsets)) &
            call resize(offsets,1,noffsets)
        
        i = 1
        do istrip=1,noffsets
            poffsets(istrip) = i
            offsets(istrip) = lbound(trace%strips(istrip)%data,1)
            i = i + size(trace%strips(istrip)%data)
        end do
        
        npacked = i-1
        
        if (.not. allocated(packed)) &
            call resize(packed,1,npacked)
            
        if (npacked > size(packed)) &
            call resize(packed,1,npacked)

        do istrip=1,noffsets
            packed(poffsets(istrip): &
                   poffsets(istrip)+size(trace%strips(istrip)%data)-1) = &
                trace%strips(istrip)%data(:)
        end do
        
    end subroutine
    
    subroutine trace_from_storable( trace, packed, poffsets, offsets )
    
      ! do inverse of trace_to_storable
      
        type(t_trace), intent(inout) :: trace
        real, dimension(:), intent(in) :: packed
        integer, dimension(:), intent(in) :: poffsets, offsets
        
        integer :: istrip, n
        
        call trace_destroy(trace)
    
        trace%nstrips = size(poffsets)
        allocate( trace%strips( trace%nstrips ) )
        
        do istrip=1,trace%nstrips
            if (istrip /= trace%nstrips) then
                n = poffsets(istrip+1)-poffsets(istrip)
            else
                n = size(packed) - poffsets(istrip) + 1
            end if
            call resize(trace%strips(istrip)%data,offsets(istrip),n)
            trace%strips(istrip)%data(:) = packed(poffsets(istrip):&
                                                  poffsets(istrip)+n-1)
        end do
        
        trace%span(1) = lbound(trace%strips(1)%data,1)
        trace%span(2) = ubound(trace%strips(trace%nstrips)%data,1)
 
    end subroutine
    
    pure function trace_size_bytes( trace )

        type(t_trace), intent(in) :: trace
        integer  :: trace_size_bytes

        integer :: istrip

        trace_size_bytes = 0
        if ( allocated(trace%strips) ) then
            do istrip=1,trace%nstrips
                trace_size_bytes = trace_size_bytes + &
                            size(trace%strips(istrip)%data) * real_kind
            end do
        end if

    end function

    elemental subroutine trace_destroy( trace )
    
      ! deallocate all strips in a trace
      
        type(t_trace), intent(inout) :: trace
        integer :: istrip
        
        if ( allocated(trace%strips) ) then
            
            do istrip=1,trace%nstrips
                call resize(trace%strips(istrip)%data,0,0 )
            end do
            deallocate( trace%strips )
        
        end if
        trace%nstrips = 0
        trace%span(:) = 0
        
    end subroutine
    
    pure logical function trace_is_empty( trace )
    
        type(t_trace), intent(in) :: trace    
    
        trace_is_empty = .not. allocated(trace%strips)
  
    end function
    
    pure subroutine intersection( a, b, c )
    
        integer, dimension(2), intent(in) :: a, b
        integer, dimension(2), intent(out) :: c
        
        c(1) = max(a(1),b(1))
        c(2) = min(a(2),b(2))
        
    end subroutine

end module
