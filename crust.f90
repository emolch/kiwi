! $Id: gfdb_build.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


program crust
 
    use util
    use better_varying_string
    use varying_string_getarg
    use crust2x2
    use orthodrome

    ! use f90_unix_env
    
    implicit none
    
    character, parameter                     :: eol = char(10)
    type(t_geo_coords)                       :: location
    type(t_crust2x2_1d_profile)              :: profile
    real :: lat, lon
    type(varying_string)                     :: aux_path
    logical                        :: ok
   

    g_pn = 'crust'
    g_usage = 'usage: ' // g_pn // ' lat lon'
    
    if (iargc() /= 2) call usage()

    call vs_getenv( 'KIWI_HOME', aux_path )
    if (aux_path == '') then
        call die('Environment variable KIWI_HOME not set.' //eol // eol // &
                 'KIWI_HOME should contain the path to the directory' // eol // &
                 'containing the "aux" directory where minimizer looks for crustal' // eol // &
                 'models etc.' // eol // &
                 'An initial "aux" directory resides in the directory where' // eol // &
                 'you have built the kiwi tools, so you can set' // eol // &
                 'KIWI_HOME to your build directory.')
    end if
    aux_path = aux_path // '/aux'
    
    call crust2x2_load( char(aux_path // '/crust2x2'), ok )
    if (.not. ok) then
        call die()
    end if

    call real_getarg( 1, -90., 90., lat )
    call real_getarg( 2, -180., +180., lon )
    location%lat = lat
    location%lon = lon    

    call crust2x2_get_profile(location, profile)
    call crust2x2_describe_profile(profile, stdout)

    call cleanup()
    
end program
