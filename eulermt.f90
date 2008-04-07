! $Id: eulermt.f90 658 2007-08-03 12:48:49Z sebastian $ 
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


program eulermt
    
    
    use euler
    use orthodrome
    
    implicit none
    
    real dip, strike, rake
    real, dimension(3,3) :: m_unrot = reshape((/0,0,-1,0,0,0,-1,0,0/),(/3,3/))
    real, dimension(3,3) :: m_rot, trotmat, rotmat, m_use
    real, dimension(3,3) :: ned2use = reshape((/0,0,-1,-1,0,0,0,1,0/),(/3,3/))
    real, dimension(3,3) :: use2ned
    integer :: iostat
    
    use2ned = transpose( ned2use )
    
    do  ! until end-of-file on stdin

        read (unit=*,fmt=*,iostat=iostat) strike, dip, rake
        if (iostat /= 0) exit
        
        call init_euler(d2r(dip),d2r(strike),d2r(-rake), rotmat)
        trotmat = transpose(rotmat)
        m_rot = matmul( rotmat, matmul( m_unrot, trotmat ) )
        
        write (*,'(6g15.5)') m_rot(1,1), m_rot(2,2), m_rot(3,3), m_rot(1,2), m_rot(1,3), m_rot(2,3)
        
        m_use = matmul( ned2use, matmul( m_rot, use2ned ) )
        
        write (*,'(6g15.5)') m_use(1,1), m_use(2,2), m_use(3,3), m_use(1,2), m_use(1,3), m_use(2,3)

    end do
        
end program

