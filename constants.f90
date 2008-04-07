! $Id: constants.f90 658 2007-08-03 12:48:49Z sebastian $ 
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

module constants

    implicit none

    real, parameter    :: pi          = 3.14159265358979
    real*8, parameter  :: pi_         = 3.14159265358979
    real, parameter    :: earthradius = 6371.*1000.
    real, parameter    :: earthradius_equator = 6378.14 * 1000.
    real*8, parameter  :: earth_oblateness = 1./298.257223563     ! WGS84

end module
