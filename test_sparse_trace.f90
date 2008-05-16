! $Id: test_sparse_trace.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


program test_sparse_trace

    use util
    use sparse_trace

    implicit none

    type(t_strip) :: cont1, cont2, accu1, accu2 !, extender
    type(t_trace) :: sparse1, sparse2, joined, empty
    
    call test_begin( "test_sparse_trace" )
    
    call strip_init( (/21,40/), (/0.,0.,0.,1.,1.,  1.,0.,0.,0.,0.,  0.,0.,1.,0.,0.,  0.,0.,0.,1.,0. /), cont1 )
    call strip_init( (/51,51/), (/6./), cont2 )
    
    call trace_destroy( empty )
    
    call trace_pack( cont1, sparse1 )
    call trace_pack( cont2, sparse2 )
 
    if (sparse1%nstrips /= 2) &
        call test_fail("nstrips")
    
    if (any(strip_span(sparse1%strips(1)) .ne. (/24,27/) )) &
        call test_fail("span1")
    if (any(strip_span(sparse1%strips(2)) .ne. (/33,40/) )) &
        call test_fail("span2")
        
    call trace_join( sparse1, sparse2, joined )
         
    if (any( joined%span .ne. (/24,51/) )) &
        call test_fail("span3")
        
    call trace_unpack( joined, cont1 )
   
    if (any( strip_span(cont1) .ne. (/24,51/) )) &
        call test_fail("span4")
        
    if (cont1%data(ubound(cont1%data,1)) /= 6.) &
        call test_fail("last value")
    
    call trace_multiply_add( sparse1, accu1 )
    call trace_multiply_add( sparse2, accu1 )
    
    call trace_multiply_add( sparse2, accu2 )
    call trace_multiply_add( sparse1, accu2 )
    if (any( cont1%data  /= accu1%data) ) &
        call test_fail("muliply-add 1")
    
    if (any( cont1%data  /= accu2%data) ) &
        call test_fail("muliply-add 2")
    
    
    call strip_init( (/-2,2/), (/1.,0.,0.,0.,1./), cont1 )
    call strip_init( (/1,4/), (/3.,1.,1.,99./), cont2 )
    call trace_pack( cont2, sparse2 )
    if (any (strip_span(sparse2%strips(1)) .ne. (/1,4/) ) ) &
        call test_fail("pack (span)")
    if (any (sparse2%strips(1)%data .ne. (/3.,1.,1.,99./) ) ) &
        call test_fail("pack (data)")
    
    
    call strip_init( (/1,2/), (/1.,1./), cont1 )
    call strip_init( (/2,3/), (/1.,1./), cont2 )
    call trace_pack( cont2, sparse2 )
    call trace_multiply_add( sparse2, cont1,1.,itraceshift_=-1 )
    if (any( strip_span(cont1) .ne. (/1,2/)) .or. &
        any( cont1%data .ne. (/2.,2./))) &
        call test_fail("multiply-add 4")
        
    call strip_init( (/-2,2/), (/1.,0.,0.,0.,1./), cont1 )
    call trace_pack( cont1, sparse1 )
    call trace_join( empty, sparse1, sparse2 )
    call trace_unpack( sparse2, cont2 )
    if (any( cont2%data .ne. cont1%data )) &
        call test_fail("join w empty")
    
    !call strip_init( (/2,6/), (/1.,2.,3.,4.,5./), extender)
    !call strip_extend( extender, (/3,5/) )
    !if (any( extender%data .ne. (/2.,3.,4./))) &
    !    call test_fail("extend shrink data")
    
    !if (any( strip_span( extender) .ne. (/3,5/) )) &
    !    call test_fail("extend shrink span")
    
    call strip_init( (/1,1/), (/0./), cont1 )
    call strip_init( (/2,4/), (/1.,1.,0./), cont2 )
    call trace_pack( cont2, sparse2 )
    call trace_multiply_add( sparse2, cont1, 1., rtraceshift_=-0.25 )
    if (any( cont1%data .ne. (/0.25,1.,0.75,0./))) &
        call test_fail("multiply-add 6")
    

    call strip_init( (/-2,5/), (/0.,0.,1.,2.,2.,2.,2.,2./), cont1 )
    if (any( strip_dataspan(cont1) .ne. (/0,1/))) &
        call test_fail("strip dataspan 1")
    
    call strip_init( (/-2,5/), (/1.,1.,1.,2.,2.,2.,2.,3./), cont1 )
    if (any( strip_dataspan(cont1) .ne. (/-2,5/))) &
        call test_fail("strip dataspan 2")

    call strip_init( (/-2,0/), (/0.,0.,0./), cont1 )
    if (any( strip_dataspan(cont1) .ne. (/0,-2/))) &
        call test_fail("strip dataspan 3")


    call test_end()
    call cleanup()
end program
