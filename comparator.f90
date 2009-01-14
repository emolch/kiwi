! $Id: comparator.f90 703 2008-04-03 15:51:31Z sebastian $ 
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

module comparator

  ! this module provides methods to compare seismograms (or any real-valued array)
  ! to a reference seismogram

    use util
    use sparse_trace
    use piecewise_linear_function
    use better_varying_string
    
    implicit none
    
    private
    
    include 'fftw3.f'
    
    integer, public, parameter :: n_comparator_norms = 5
    
    integer, public, parameter :: L2NORM = 1
    integer, public, parameter :: L1NORM = 2
    integer, public, parameter :: AMPSPEC_L2NORM = 3
    integer, public, parameter :: AMPSPEC_L1NORM = 4
    integer, public, parameter :: SCALAR_PRODUCT = 5

  ! processing ids
    integer, parameter, public :: PLAIN = 1
    integer, parameter, public :: TAPERED = 2
    integer, parameter, public :: FILTERED = 3
    
    type(varying_string), private, dimension(:), allocatable :: comparator_norm_names
    
    integer(8), dimension(32)  :: fftwplans = 0
    integer(8), dimension(32)  :: fftwplans_inv = 0

    
    type, public :: t_probe
    
        real                                  :: dt, df
        integer, dimension(2)                 :: span
        integer, dimension(2)                 :: dataspan ! left outside of dataspan lie only zeros, right the last value is repeated.

        real(4), dimension(:), allocatable    :: array
        real(4), dimension(:), allocatable    :: array_tapered
        complex(4), dimension(:), allocatable :: spectrum
        complex(4), dimension(:), allocatable :: spectrum_filtered
        real(4), dimension(:), allocatable    :: amp_spectrum
        real(4), dimension(:), allocatable    :: amp_spectrum_filtered
        real(4), dimension(:), allocatable    :: array_filtered

        logical                               :: array_dirty = .true.
        logical                               :: array_tapered_dirty = .true.
        logical                               :: spectrum_dirty = .true.
        logical                               :: spectrum_filtered_dirty = .true.
        logical                               :: array_filtered_dirty = .true.

        real                                  :: paddingfactor
        type(t_plf)                           :: taper
        type(t_plf)                           :: filter
        real                                  :: factor = 1.

    end type
    
    public comparator_get_norm_name
    public comparator_get_norm_id
    
    public probe_init
    public probe_destroy
    public probe_set_array
    public probe_extend_span
    public probe_get_amp_spectrum
    public probe_get_plain
    public probe_get_tapered
    public probe_get_filtered
    public probe_get
    public probe_shift

    public probe_set_taper
    public probe_set_filter
    public probe_set_factor

  
    public probes_adjust_spans
    public probes_norm
    public probe_norm
    public probes_windowed_cross_corr
  
    public cleanup_comparator
  
  contains
  
    subroutine cleanup_comparator()
    
        integer :: i
        
        do i=1,32
            if (fftwplans(i) .ne. 0) then
                call sfftw_destroy_plan( fftwplans(i) )
                fftwplans(i) = 0
            end if
            if (fftwplans_inv(i) .ne. 0) then
                call sfftw_destroy_plan( fftwplans_inv(i) )
                fftwplans_inv(i) = 0
            end if
        end do
        
        if (allocated( comparator_norm_names )) then
            do i=1,n_comparator_norms
                call delete( comparator_norm_names(i) )
            end do
            deallocate( comparator_norm_names )
        end if
    
    end subroutine
    
    subroutine comparator_init_norm_names()
    
        if (.not. allocated( comparator_norm_names )) then
            allocate( comparator_norm_names( n_comparator_norms ) )
            comparator_norm_names(1) = "l2norm"
            comparator_norm_names(2) = "l1norm"
            comparator_norm_names(3) = "ampspec_l2norm"
            comparator_norm_names(4) = "ampspec_l1norm"
            comparator_norm_names(5) = "scalar_product"
        end if
    
    end subroutine
    
    subroutine comparator_get_norm_name( id, name )
    
        integer, intent(in)                 :: id
        type(varying_string), intent(inout) :: name
        
        call comparator_init_norm_names()
        
        name = ""
        if (id < 1 .or. n_comparator_norms < id) return
        name = comparator_norm_names(id)
        
    end subroutine
   
    subroutine comparator_get_norm_id( name, id )
    
        type(varying_string), intent(in)     :: name
        integer, intent(out)                 :: id

        integer :: i
        
        call comparator_init_norm_names()
        
        id = 0
        do i=1, n_comparator_norms
            if (name .eq. comparator_norm_names(i)) then
                id = i
                return
            end if
        end do
    
    end subroutine
  
    subroutine probe_init( self, dt, paddingfactor )
    
        type(t_probe), intent(inout) :: self
        real, intent(in)             :: dt
        real, intent(in), optional   :: paddingfactor
        
        call probe_destroy( self )
        
        self%dt = dt
        self%span = (/0,0/)
        self%dataspan = (/0,0/)
        call dirtyfy_array( self )
        self%paddingfactor = 2.
        if (present(paddingfactor)) self%paddingfactor = paddingfactor
        self%factor = 1.
        
    end subroutine
    
    subroutine probe_destroy( self )
    
        type(t_probe), intent(inout) :: self
        call dirtyfy_array(self)
        if (allocated( self%array )) deallocate(self%array)
        if (allocated( self%array_tapered )) deallocate(self%array_tapered)
        if (allocated( self%spectrum )) deallocate(self%spectrum)
        if (allocated( self%spectrum_filtered )) deallocate(self%spectrum_filtered)
        if (allocated( self%amp_spectrum )) deallocate(self%amp_spectrum)
        if (allocated( self%amp_spectrum_filtered )) deallocate(self%amp_spectrum_filtered)
        if (allocated( self%array_filtered )) deallocate(self%array_filtered)
        call plf_destroy( self%taper )
        call plf_destroy( self%filter )
        
    end subroutine
    
    subroutine probe_set_array( self, strip, allow_shrink_, span_hint, factor_ )
    
        type(t_probe), intent(inout) :: self
        type(t_strip), intent(in) :: strip
        logical, intent(in), optional :: allow_shrink_ ! defaults to .false.
        integer, dimension(2), intent(in), optional :: span_hint
        real, intent(in), optional :: factor_
        
        integer, dimension(2) :: newspan, tempspan
        logical :: allow_shrink
        integer :: datalength
        real :: factor
        
        allow_shrink = .false.
        if (present( allow_shrink_ )) allow_shrink = allow_shrink_

        factor = 1.
        if (present( factor_ )) factor = factor_        

        self%dataspan = strip_span( strip )

      ! resize strip such that length is a power of two
      ! add padding to both sides to reduce probabability, that a later resize is needed.
        if (allow_shrink .or. .not. allocated(self%array)) then
            newspan = strip_span(strip)
        else
            call union( strip_span(strip), self%span, newspan )
        end if
        if (present(span_hint)) then
            tempspan(:) = newspan(:)
            call union(tempspan, span_hint, newspan)
        end if
        
        datalength = slen(self%dataspan)
        newspan = allowed_span(newspan,ceiling(datalength*self%paddingfactor))

        if (any(newspan /= self%span)) then
            call resize( self%array,         newspan(1), slen(newspan) )
            call resize( self%array_tapered, newspan(1), slen(newspan) )
        end if
        self%span = newspan
        
        if (self%span(1) <= self%dataspan(1)-1) self%array(:self%dataspan(1)-1) = 0.
        self%array(self%dataspan(1):self%dataspan(2)) = strip%data(:) * factor
        if (self%dataspan(2)+1 <= self%span(2)) self%array(self%dataspan(2)+1:) = &
                                                self%array(self%dataspan(2))

        call dirtyfy_array( self )
                
    end subroutine
    
    subroutine probe_shift( self, ishift )
    
        type(t_probe), intent(inout) :: self
        integer, intent(in) :: ishift
        
        type(t_strip) :: strip
        integer, dimension(2) :: newdataspan

        if (.not. allocated(self%array)) return

        newdataspan(:) = self%dataspan(:) + ishift
        call strip_init( newdataspan, self%array(self%dataspan(1):self%dataspan(2)), strip )
        call probe_set_array( self, strip )
        call strip_destroy( strip )

    end subroutine

    
    subroutine probe_extend_span( self, span )
    
      ! extend span to include 'span'
    
        type(t_probe), intent(inout) :: self
        integer, dimension(2), intent(in) :: span
        
        integer, dimension(2) :: newspan
        type(t_strip) :: temp
        
        if (.not. allocated(self%array)) then 
            newspan = allowed_span( span, 0 )
            call resize( self%array, span(1), slen(span) )
            self%array(:) = 0.
            call resize( self%array_tapered, span(1), slen(span) )
            return
        end if

        call union( span, self%dataspan, newspan )
        newspan = allowed_span( newspan, 0 )
        if (all(self%span .eq. newspan)) return
        
      ! backup data, resize, reinsert data        
        call strip_init( self%dataspan, self%array(self%dataspan(1):self%dataspan(2)), temp )
        call resize( self%array, newspan(1), slen(newspan) )
        call resize( self%array_tapered, newspan(1), slen(newspan) )
        self%span = newspan

        
        if (self%span(1) <= self%dataspan(1)-1) &
            self%array(:self%dataspan(1)-1) = 0.
        self%array(self%dataspan(1):self%dataspan(2)) = temp%data(:)
        if (self%dataspan(2)+1 <= self%span(2)) &
            self%array(self%dataspan(2)+1:) = self%array(self%dataspan(2))
        
        call strip_destroy(temp)
        call dirtyfy_array( self )
        
        
    end subroutine
    

    subroutine probe_get_amp_spectrum( self, strip, df, which_processing )
    
        type(t_probe), intent(inout) :: self
        type(t_strip), intent(inout) :: strip
        real, intent(out) :: df
        integer, intent(in) :: which_processing
        
        call update_spectrum_filtered( self )
        
        if (.not. allocated( self%amp_spectrum )) then
            call strip_init( (/1,1/), (/0./), strip )
            return
        end if
        
        df = self%df
        if (plf_defined( self%filter ) .and. which_processing == FILTERED) then
            call strip_init( (/1,size(self%amp_spectrum_filtered)/), self%amp_spectrum_filtered, strip )
        else
            call strip_init( (/1,size(self%amp_spectrum)/), self%amp_spectrum, strip )
        end if
        
    end subroutine
    
    subroutine probe_get_plain( self, strip )

        type(t_probe), intent(inout) :: self
        type(t_strip), intent(inout) :: strip

        call update_array( self )
        
        if (.not. allocated( self%array )) then
            call strip_init( (/1,1/), (/0./), strip )
            return
        end if
        call strip_init( self%dataspan, self%array(self%dataspan(1):self%dataspan(2)), strip )

    end subroutine

    subroutine probe_get_tapered( self, strip )
    
        type(t_probe), intent(inout) :: self
        type(t_strip), intent(inout) :: strip
        
        integer, dimension(2) :: span

        call update_array_tapered( self )

        if (plf_defined( self%taper )) then
            if (.not. allocated( self%array_tapered )) then
                call strip_init( (/1,1/), (/0./), strip )
                return
            end if
            call intersection(discrete_plf_span(self%taper,self%dt), self%dataspan, span)
            if (span(1) > span(2)) span = self%dataspan
            call strip_init( span, self%array_tapered(span(1):span(2)), strip )
        else
            call probe_get_plain( self, strip )
        end if
        
    end subroutine

    subroutine probe_get_filtered( self, strip )

        type(t_probe), intent(inout) :: self
        type(t_strip), intent(inout) :: strip

        integer, dimension(2) :: span

        call update_array_filtered( self )

        if (plf_defined( self%filter )) then
            if (.not. allocated( self%array_filtered )) then
                call strip_init( (/1,1/), (/0./), strip )
                return
            end if
            if (plf_defined( self%taper )) then
                call intersection(discrete_plf_span(self%taper,self%dt), self%span, span)
                if (span(1) > span(2)) span = self%dataspan
            else
                span = self%dataspan
            end if
            call strip_init( span, self%array_filtered(span(1):span(2)), strip )
        else 
            call probe_get_tapered( self, strip )
        end if

    end subroutine

    subroutine probe_get( self, strip, which_processing )
        
        type(t_probe), intent(inout) :: self
        type(t_strip), intent(inout) :: strip
        integer, intent(in) :: which_processing

      ! according to `which_processing` call one of the probe_get_* functions

        if (which_processing == PLAIN) call probe_get_plain( self, strip )
        if (which_processing == TAPERED) call probe_get_tapered( self, strip )
        if (which_processing == FILTERED) call probe_get_filtered( self, strip )

    end subroutine
    
    subroutine probe_set_taper( self, plf )
    
        type(t_probe), intent(inout) :: self
        type(t_plf), intent(in) :: plf
        
        call plf_copy( plf, self%taper )
        call dirtyfy_array_tapered( self )
        
    end subroutine
    
    subroutine probe_set_filter( self, plf )
    
        type(t_probe), intent(inout) :: self
        type(t_plf), intent(in) :: plf
        
        call plf_copy( plf, self%filter )
        call dirtyfy_spectrum_filtered( self )
    
    end subroutine
    
    subroutine probe_set_factor( self, factor )

        type(t_probe), intent(inout) :: self
        real, intent(in) :: factor
        
        self%factor = factor

    end subroutine

    subroutine probes_adjust_spans( a, b )
    
        type(t_probe), intent(inout) :: a, b
        
        integer, dimension(2) :: newspan
        integer :: minlength
        
        if (a%dt /= b%dt) then
            call die("probes_adjust_spans(): both probes must have same dt")
        end if
        
        call union( a%dataspan, b%dataspan, newspan )
        minlength = max( ceiling(slen(a%dataspan)*a%paddingfactor), ceiling(slen(b%dataspan)*b%paddingfactor) )
        
        newspan = allowed_span( newspan, minlength )

        if (all(a%span == b%span) .and. slen(a%span) == slen(newspan) .and. &
             containing(a%span, b%dataspan) .and. containing(b%span, a%dataspan)) return
        
        call probe_extend_span( a, newspan )
        call probe_extend_span( b, newspan )
                
    end subroutine

    pure function scalar_product_2( a, b, dt, fa, fb ) result(prod)
        real, dimension(:), intent(in) :: a,b
        real, intent(in) :: dt, fa, fb
        real :: prod
        prod = dt * 0. ! (get rid of compiler warning)
        if (fa == 1. .and. fb == 1.) then
            prod = real(sum(real(a*b,8)))
        else
            prod = real(sum(real(a*fa*b*fb,8)))
        end if
    end function

    pure function l1norm_func( a, b, dt, fa, fb ) result(norm)
        real, dimension(:), intent(in) :: a,b
        real, intent(in) :: dt, fa, fb
        real :: norm
        if (fa == 1. .and. fb == 1.) then
            norm = real(dt*sum(real(abs(a-b),8)))
        else
            norm = real(dt*sum(real(abs(fa*a-fb*b),8)))
        end if
    end function

    pure function l2norm_func( a, b, dt, fa, fb ) result(norm)
        real, dimension(:), intent(in) :: a,b
        real, intent(in) :: dt, fa, fb
        real :: norm
        if (fa == 1. .and. fb == 1.) then
            norm = real(sqrt(dt*sum(real(a-b,8)**2)))
        else
            norm = real(sqrt(dt*sum(real(fa*a-fb*b,8)**2)))
        end if
    end function
    

    pure function scalar_product_1( a, dt, fa ) result(prod)
        real, dimension(:), intent(in) :: a
        real, intent(in) :: dt, fa
        real :: prod
        prod = dt * 0. ! (get rid of compiler warning)
        prod = fa**2 * real(sum(real(a*a,8)))
    end function

    pure function l1norm_func_1( a, dt, fa ) result(norm)
        real, dimension(:), intent(in) :: a
        real, intent(in) :: dt, fa
        real :: norm
        norm = fa * real(dt*sum(real(abs(a),8)))
    end function

    pure function l2norm_func_1( a, dt, fa ) result(norm)
        real, dimension(:), intent(in) :: a
        real, intent(in) :: dt, fa
        real :: norm
        norm = fa * real(sqrt(dt*sum(real(a,8)**2)))
    end function
    
    function probes_norm_timedomain( a, b, normfunction ) result(norm)

        interface
            pure function normfunction( a, b, dt, fa, fb ) result(norm)
                real, dimension(:), intent(in) :: a,b
                real, intent(in) :: dt, fa, fb
                real :: norm
            end function
        end interface

        type(t_probe), intent(inout) :: a, b
        real :: norm
        
        integer, dimension(2) :: span, a_temp_span, b_temp_span

        call probes_adjust_spans( a, b )
     
      ! when tapers are specified, restrict application to taper span
      ! otherwise, restrict to datarange
        if (plf_defined( a%taper ) .and. plf_defined( b%taper )) then
            call intersection(discrete_plf_span(a%taper,a%dt), a%span, a_temp_span)
            call intersection(discrete_plf_span(b%taper,b%dt), b%span, b_temp_span)
            if (a_temp_span(1) > a_temp_span(2)) then
                span = b_temp_span
            else if (b_temp_span(1) > b_temp_span(2)) then
                span = a_temp_span
            else
                call union( a_temp_span, b_temp_span, span )
            end if
        else
            call probes_adjust_spans( a, b )
            call union( a%dataspan,b%dataspan,span )      
        end if

        if (span(1) > span(2)) then
            call warn( 'comparator: applying timedomain norm to empty region.' )
            norm = 0.
            return
        end if
        
        if (plf_defined( a%filter ) .and. plf_defined( b%filter )) then
            call update_array_filtered( a )
            call update_array_filtered( b )
            norm = normfunction(a%array_filtered(span(1):span(2)), b%array_filtered(span(1):span(2)), a%dt, a%factor, b%factor )
        else if (plf_defined( a%taper ) .and. plf_defined( b%taper )) then
            call update_array_tapered( a )
            call update_array_tapered( b )
            norm = normfunction(a%array_tapered(span(1):span(2)), b%array_tapered(span(1):span(2)), a%dt, a%factor, b%factor)
        else
            norm = normfunction(a%array(span(1):span(2)), b%array(span(1):span(2)), a%dt, a%factor, b%factor)
        end if

    end function
    
    function probe_norm_timedomain( a, normfunction ) result(norm)

        interface
            pure function normfunction( a, dt, fa ) result(norm)
                real, dimension(:), intent(in) :: a
                real, intent(in) :: dt, fa
                real :: norm
            end function
        end interface

        type(t_probe), intent(inout) :: a
        real :: norm
        
        integer, dimension(2) :: span
        
        call intersection(discrete_plf_span(a%taper,a%dt), a%span, span)
        
      ! when taper is specified, restrict application to taper span
      ! otherwise, restrict to datarange
        if (plf_defined( a%taper )) then
            call intersection(discrete_plf_span(a%taper,a%dt), a%span, span)
        else 
            span = a%dataspan
        end if

        if (plf_defined( a%filter )) then
            call update_array_filtered( a )
            norm = normfunction(a%array_filtered(span(1):span(2)), a%dt, a%factor)
        else if (plf_defined( a%taper )) then
            call update_array_tapered( a )
            norm = normfunction(a%array_tapered(span(1):span(2)), a%dt, a%factor)
        else
            norm = normfunction(a%array(span(1):span(2)), a%dt, a%factor)
        end if

    end function

    function probes_norm_frequencydomain( a, b, normfunction ) result(norm)

        interface
            pure function normfunction( a, b, dt, fa, fb ) result(norm)
                real, dimension(:), intent(in) :: a,b
                real, intent(in) :: dt, fa, fb
                real :: norm
            end function
        end interface

        type(t_probe), intent(inout) :: a, b
        real :: norm

        call probes_adjust_spans( a, b )

        if (plf_defined( a%filter ) .and. plf_defined( b%filter)) then
            call update_spectrum_filtered( a )
            call update_spectrum_filtered( b )
            norm = normfunction( a%amp_spectrum_filtered, b%amp_spectrum_filtered, a%df, a%factor, b%factor)
        else
            call update_spectrum( a )
            call update_spectrum( b )
            norm = normfunction( a%amp_spectrum, b%amp_spectrum, a%df, a%factor, b%factor )
        end if

    end function
    
    function probe_norm_frequencydomain( a, normfunction ) result(norm)

        interface
            pure function normfunction( a, dt, fa ) result(norm)
                real, dimension(:), intent(in) :: a
                real, intent(in) :: dt, fa
                real :: norm
            end function
        end interface

        type(t_probe), intent(inout) :: a
        real :: norm

        if (plf_defined( a%filter )) then
            call update_spectrum_filtered( a )
            norm = normfunction( a%amp_spectrum_filtered, a%df, a%factor)
        else
            call update_spectrum( a )
            norm = normfunction( a%amp_spectrum, a%df, a%factor)
        end if

    end function
    
    function probes_norm( a, b, inmethod ) result(norm)
        
        type(t_probe), intent(inout) :: a, b
        integer, intent(in), optional :: inmethod
        
        real :: norm
        integer :: method
        
        method = L2NORM
        if (present( inmethod )) method = inmethod
        select case(method)
        
            case (L2NORM)
                norm = probes_norm_timedomain( a, b, l2norm_func )
                return
                
            case (L1NORM)
                norm = probes_norm_timedomain( a, b, l1norm_func )
                return

            case (SCALAR_PRODUCT)
                norm = probes_norm_timedomain( a, b, scalar_product_2 )
                return

            case (AMPSPEC_L2NORM)
                norm = probes_norm_frequencydomain( a, b, l2norm_func )
                return
                
            case (AMPSPEC_L1NORM)
                norm = probes_norm_frequencydomain( a, b, l1norm_func )
                return
            
            case default
                call die( "probes_norm(): unknown norm method" )
            
        end select
    
    end function
    
    function probe_norm( a, inmethod ) result(norm)
        
        type(t_probe), intent(inout) :: a
        integer, intent(in), optional :: inmethod
        
        real :: norm
        integer :: method
        
        method = L2NORM
        if (present( inmethod )) method = inmethod
        
        select case(method)
        
            case (L2NORM)
                norm = probe_norm_timedomain( a, l2norm_func_1 )
                return
                
            case (L1NORM)
                norm = probe_norm_timedomain( a, l1norm_func_1 )
                return
            
            case (SCALAR_PRODUCT)
                norm = probe_norm_timedomain( a, scalar_product_1 )
                return

            case (AMPSPEC_L2NORM)
                norm = probe_norm_frequencydomain( a, l2norm_func_1 )
                return
            
            case (AMPSPEC_L1NORM)
                norm = probe_norm_frequencydomain( a, l1norm_func_1 )
                return
            
            case default
                call die( "probe_norm(): unknown norm method" )
            
        end select
    
    end function

    function probes_scalar_product( a, b ) result(prod)
        type(t_probe), intent(inout) :: a, b
        real :: prod
        prod = probes_norm_timedomain( a, b, scalar_product_2 )
    end function
    
    subroutine probes_windowed_cross_corr( a, b, shiftrange, cross_corr )

        type(t_probe), intent(inout) :: a, b
        integer, intent(in), dimension(2) :: shiftrange
        real, intent(out), dimension(:) :: cross_corr

      ! Get cross correlation between a and b, which may be windowed.
      !
      ! If a and b are windowed, the window position is fixed at a and
      ! b is "pulled through it's window."
      ! Shifts in the range `shiftrange` are tested.
      ! The cross correlation is returned in `cross_corr`, which must be
      ! an array long enough to hold shiftrange(2)-shiftrange(1)+1 values
      ! The returned values should be interpreted as how good a and b 
      ! match when b is shifted by the corresponding number of elements.

        integer :: ishift, i

        ishift = shiftrange(1)
        
        do i=1,slen(shiftrange)
            call probe_shift( b, ishift )
            ishift = 1
            cross_corr(i) = probes_scalar_product( a,b )
        end do

      ! reset shift
        call probe_shift( b, -shiftrange(2) )

    end subroutine

    pure function allowed_span( span, minlength ) result(newspan)
    
        integer, dimension(2), intent(in) :: span
        integer, intent(in)               :: minlength
        integer, dimension(2)             :: newspan
        
        integer :: length, lengthp
        
        newspan(:) = span(:)
        length = slen(newspan)
        if (length < minlength) then
            length = minlength
        end if
        lengthp = next_power_of_two( length )
        newspan(1) = newspan(1) - floor((lengthp-slen(span))/2.)
        newspan(2) = newspan(1) + lengthp - 1
        
    end function
    
    pure function next_power_of_two(n) result(m)
    
        integer, intent(in) :: n
        integer :: m
        
        m = 2**ceiling(log(real(n))/log(2.))
    
    end function
    
    pure function containing( outer, inner )

        integer, dimension(2), intent(in) :: outer, inner
        logical :: containing
        integer, dimension(2) :: intersect
        
        call intersection(outer, inner, intersect)
        containing = all( intersect == inner )

    end function

    pure subroutine intersection( a, b, c )
    
        integer, dimension(2), intent(in) :: a, b
        integer, dimension(2), intent(out) :: c
        
        c(1) = max(a(1),b(1))
        c(2) = min(a(2),b(2))
        
    end subroutine

    pure subroutine union( a, b, c ) 
    
        integer, dimension(2), intent(in) :: a, b
        integer, dimension(2), intent(out) :: c
        
        c(1) = min(a(1),b(1))
        c(2) = max(a(2),b(2))
        
    end subroutine
    
    pure function slen( span )
        integer, dimension(2), intent(in) :: span
        integer :: slen
        slen = span(2) - span(1) + 1
    end function
    
    pure function discrete_plf_span( plf, dt )

        type(t_plf), intent(in) :: plf
        real, intent(in) :: dt
        integer, dimension(2) :: discrete_plf_span

        real, dimension(2) :: rspan

        rspan = plf_span( plf )
        discrete_plf_span(1) = ceiling(rspan(1)/dt)
        discrete_plf_span(2) = floor(rspan(2)/dt)
       
    end function

  ! calculation routines
    
    subroutine make_array_tapered( self )
            
        type(t_probe), intent(inout) :: self
    
        if (plf_defined( self%taper )) then
            self%array_tapered(:) = self%array(:)
            call plf_taper_array( self%taper, &
                                  self%array_tapered(self%dataspan(1):self%span(2)), &
                                  (/self%dataspan(1),self%span(2)/), self%dt, ip_cos )
        end if
            
    end subroutine
    
    subroutine make_spectrum( self )
        
        type(t_probe), intent(inout) :: self
        
        integer :: ntrans, e
     

        ntrans = size(self%array)

        e = int(ceiling(log(real(ntrans))/log(2.)))

        call resize( self%spectrum, 1, ntrans/2+1 )
        call resize( self%spectrum_filtered, 1, ntrans/2+1 )
        call resize( self%amp_spectrum, 1, ntrans/2+1 )
        call resize( self%amp_spectrum_filtered, 1, ntrans/2+1 )
        if (fftwplans(e) .eq. 0) then
            call sfftw_plan_dft_r2c_1d(fftwplans(e),ntrans,self%array, &
                                       self%spectrum,FFTW_ESTIMATE+FFTW_UNALIGNED)
        end if
        
        if (plf_defined( self%taper )) then
            call sfftw_execute_dft_r2c(fftwplans(e), self%array_tapered, self%spectrum)
        else 
            call sfftw_execute_dft_r2c(fftwplans(e), self%array, self%spectrum)
        end if
        
        self%amp_spectrum = abs(self%spectrum)
        
        self%df = 1./(ntrans * self%dt)
    
    end subroutine
    
    subroutine make_spectrum_filtered( self )
    
        type(t_probe), intent(inout) :: self
    
        if (plf_defined(self%filter)) then
            self%amp_spectrum_filtered(:) = self%amp_spectrum(:)
            self%spectrum_filtered(:) = self%spectrum(:)
            call plf_taper_array( self%filter, self%spectrum_filtered, &
                 (/0,size(self%spectrum_filtered)-1/), self%df, ip_cos)
            call plf_taper_array( self%filter, self%amp_spectrum_filtered, &
                 (/0,size(self%amp_spectrum_filtered)-1/), self%df, ip_cos )
        end if
    
    end subroutine
    
    subroutine make_array_filtered( self )

        type(t_probe), intent(inout) :: self

        integer :: ntrans, e

        if (plf_defined(self%filter)) then
            ntrans = size(self%array)
            call resize( self%array_filtered, lbound(self%array,1), ntrans )

            e = int(ceiling(log(real(ntrans))/log(2.)))
            if (fftwplans_inv(e) .eq. 0) then
                call sfftw_plan_dft_c2r_1d(fftwplans_inv(e),ntrans,self%spectrum_filtered, &
                                        self%array_filtered,FFTW_ESTIMATE+FFTW_UNALIGNED)
            end if

            call sfftw_execute_dft_c2r(fftwplans_inv(e), self%spectrum_filtered, self%array_filtered)
            
          ! normalize result
            self%array_filtered = self%array_filtered / ntrans

          ! make filtered trace zero, where taper is zero
            if (plf_defined(self%taper)) then
                call plf_taper_array( self%taper, &
                                  self%array_filtered(self%span(1):self%span(2)), &
                                  (/self%span(1),self%span(2)/), self%dt, ip_zero_one )
            end if

        end if

    end subroutine
    
  ! dataflow update routines
  
    subroutine update_array(self)
        type(t_probe), intent(inout) :: self
        self%array_dirty = .false.
    end subroutine
    
    subroutine update_array_tapered(self)
        type(t_probe), intent(inout) :: self
        call update_array(self)
        if (self%array_tapered_dirty) then
            call make_array_tapered(self)
        end if
        self%array_tapered_dirty = .false.
    end subroutine
    
    subroutine update_spectrum(self)
        type(t_probe), intent(inout) :: self
        call update_array_tapered(self)
        if (self%spectrum_dirty) then
            call make_spectrum(self)
        end if
        self%spectrum_dirty = .false.
    end subroutine

    subroutine update_spectrum_filtered(self)
        type(t_probe), intent(inout) :: self
        call update_spectrum(self)
        if (self%spectrum_filtered_dirty) then
            call make_spectrum_filtered(self)
        end if
        self%spectrum_filtered_dirty = .false.
    end subroutine
    
    subroutine update_array_filtered(self)
        type(t_probe), intent(inout) :: self
        call update_spectrum_filtered(self)
        if (self%array_filtered_dirty) then
            call make_array_filtered(self)
        end if
        self%array_filtered_dirty = .false.
    end subroutine

  ! dataflow dirtyfy routines
            
    subroutine dirtyfy_array(self)
        type(t_probe), intent(inout) :: self
        self%array_dirty = .true.
        call dirtyfy_array_tapered(self)
    end subroutine
    
    subroutine dirtyfy_array_tapered(self)
        type(t_probe), intent(inout) :: self
        self%array_tapered_dirty = .true.
        call dirtyfy_spectrum(self)
    end subroutine
    
    subroutine dirtyfy_spectrum(self)
        type(t_probe), intent(inout) :: self
        self%spectrum_dirty = .true.
        call dirtyfy_spectrum_filtered(self)
    end subroutine

    subroutine dirtyfy_spectrum_filtered(self)
        type(t_probe), intent(inout) :: self
        self%spectrum_filtered_dirty = .true.
        call dirtyfy_array_filtered(self)
    end subroutine
    
    subroutine dirtyfy_array_filtered(self)
        type(t_probe), intent(inout) :: self
        self%array_filtered_dirty = .true.
    end subroutine

end module

