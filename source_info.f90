! $Id: source_info.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


program source_info
  
  ! This tiny tool prints some information about the earthquake source models 
  ! implemented in seismosizer
  !
  ! usage: source_info [ source-type ]
  ! 
  ! Complete documentation is available on
  ! 
  !   http://kinherd.org/power/trac/wiki/SourceInfoTool
  !
 

    use util
    use source
    use better_varying_string
    use varying_string_getarg
    
    implicit none
    type(varying_string)        :: sourcename, paramname, paramunit
    integer                     :: nparams, sourcetype, isourcetype, iparam
    character, parameter             :: eol = char(10)
    real, allocatable,dimension(:)   :: mins, maxs, defaults

    
    g_pn = "source_info"

    g_usage = "usage: "//g_pn//" [ sourcename ]" // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/SourceInfoTool"
    
    if (iargc() .ne. 1) then
        write (*, "(a)", advance="no") "source types: "
        do isourcetype = 1,nsourcetypes
            call psm_get_source_name( sourcetypes(isourcetype), sourcename )
            write (*, "(a,a)", advance="no") char(sourcename), " "
        end do
        write (*,*)
        stop
    end if

    call vs_getarg( 1, sourcename )
    if (sourcename .eq. '--help') call usage()
    if (sourcename .eq. '-h')     call usage()

    call psm_get_source_id( sourcename, sourcetype )
    
    if (sourcetype .eq. 0) call die( "unknown sourcetype: " // sourcename )
    nparams = psm_get_n_source_params( sourcetype )

    write (*,"(a)") char("parameters: " // nparams)
    write (*,"(a)", advance="no") "parameter names: "
    do iparam = 1,nparams
        call psm_get_param_name( sourcetype, iparam, paramname )
        write (*, "(a,a)", advance="no") char(paramname), " "
    end do
    write(*,*)
    
    write (*,"(a)", advance="no") "parameter units: "
    do iparam = 1,nparams
        call psm_get_param_unit( sourcetype, iparam, paramunit )
        write (*, "(a,a)", advance="no") char(paramunit), " "
    end do
    write(*,*)
    
    call psm_get_param_limits( sourcetype, mins, maxs, .true. )
    
    write (*,"(a)", advance="no") "parameter hard min: "
    do iparam = 1,nparams
        write (*, "(es9.2,a)", advance="no") mins(iparam), " "
    end do
    write(*,*)
    
    write (*,"(a)", advance="no") "parameter hard max: "
    do iparam = 1,nparams
        write (*, "(es9.2,a)", advance="no") maxs(iparam), " "
    end do
    write(*,*)
    
    call psm_get_param_limits( sourcetype, mins, maxs, .false. )
    
    write (*,"(a)", advance="no") "parameter soft min: "
    do iparam = 1,nparams
        write (*, "(es9.2,a)", advance="no") mins(iparam), " "
    end do
    write(*,*)
    
    write (*,"(a)", advance="no") "parameter soft max: "
    do iparam = 1,nparams
        write (*, "(es9.2,a)", advance="no") maxs(iparam), " "
    end do
    write(*,*)
    
    call psm_get_param_defaults( sourcetype, defaults )
    write (*,"(a)", advance="no") "parameter defaults: "
    do iparam = 1,nparams
        write (*, "(es9.2,a)", advance="no") defaults(iparam), " "
    end do
    write(*,*)
    
    call psm_cleanup()
    call cleanup()
    
end program

