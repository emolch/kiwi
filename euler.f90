! $Id: euler.f90 683 2007-12-01 12:09:15Z sebastian $ 
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


module euler

! this module only provides a single method to initialize a rotation matrix
! given eulerian angles

    implicit none

    public init_euler

  contains

    pure subroutine init_euler( alpha, beta, gamma, mat )
        
      ! given the euler angles alpha,beta,gamma, 
      ! make matrix mat a rotation matrix
        
      ! given coordinate system (x,y,z) and rotated system (xs,ys,zs)
      ! the line of nodes is the intersection between the x-y and the xs-ys
      ! planes.
      ! alpha is the angle between the z-axis and the zs-axis.
      ! beta is the angle between the x-axis and the line of nodes.
      ! gamma is the angle between the line of nodes and the xs-axis.
      
      ! Usage for moment tensors:
      !   real, dimension(3,3) :: m_unrot = reshape((/0,0,-1,0,0,0,-1,0,0/),(/3,3/))
      !   call init_euler(dip,strike,-rake, rotmat)
      !   m = matmul( rotmat, matmul( m_unrot, transpose(rotmat) ) )
        
        real, intent(in) :: alpha, beta, gamma
        real, intent(out), dimension(3,3) :: mat
        
        real :: ca,cb,cg,sa,sb,sg
        
        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sa = sin(alpha)
        sb = sin(beta)
        sg = sin(gamma)
        
        mat(1,1) = cb*cg-ca*sb*sg
        mat(2,1) = sb*cg+ca*cb*sg
        mat(3,1) = sa*sg
        mat(1,2) = -cb*sg-ca*sb*cg
        mat(2,2) = -sb*sg+ca*cb*cg
        mat(3,2) = sa*cg
        mat(1,3) = sa*sb
        mat(2,3) = -sa*cb
        mat(3,3) = ca
        
    end subroutine
    
end module
