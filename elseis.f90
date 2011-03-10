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

module elseis
      
!  S.Cesca & S. Heimann, Univ. Hamburg, October-November 2005
!  
!  This module contains routines to calculate elemenary seimograms for single forces
!  and moment tensors in a homogeneous, isotropic fullspace, including near field terms.
!
!  The equations are taken from Aki & Richards.
!
!  public methods are:
!
!        elseis    - to calculate all elementary seismograms at once
!        elseis_mt - to calculate elementary seismograms for single moment tensor components
!        elseis_sf - to calculate elementary seismograms for single single-force components
!        radpat_mt - to calculate radiation pattern coefficients for moment tensor Green function terms
!        radpat_sf - to calculate radiation pattern coefficients for single force Green function terms
!        make_direction_cosine - to calculate the direction cosines needed by the elseis routines

    use differentiation
    use integration

    implicit none
    
    private
    public :: elseis_all, elseis_mt, elseis_sf, radpat_mt, radpat_sf, &
              material_factors_mt, material_factors_sf, factors_mt, factors_sf, &
              make_direction_cosine, make_istfs

    real, parameter :: PI = 3.14159265358979
    integer, parameter, public :: nelseism = 36
    integer, parameter, dimension(3,3) :: delta = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))

  contains

!********************************************************
!  ELEMENTARY SEISMOGRAMS
!********************************************************

     pure subroutine elseis_all(  rho,alpha,beta,                             &
                        toffset,dt,                                 &
                        stf,dstf,istf,istftau,                      &
                        coord,                                      &
                        nfflag,ffflag,mtflag,sfflag,                &
                        all_elseism )
 
        real,    intent(in)                   :: rho, alpha, beta, toffset, dt
        real,    intent(in), dimension(:)     :: stf,dstf,istf,istftau
        real,    intent(in), dimension(3)     :: coord
        logical, intent(in)                   :: nfflag,ffflag,mtflag,sfflag
        real,    intent(out), dimension(:,:)  :: all_elseism
        
        ! calculate elementary seismograms for every moment tensor and single force component
        ! the subroutine parameters are (units are assumed to be in SI)
        ! (in):
        !   rho, alpha, beta   : density, p-wave-velocity, s-wave-velocity
        !   toffset            : temporal offset of first element of output seismograms
        !   dt                 : sampling distance of source time function(s) and output
        !   stf                : source time function, stf(1) must be zero
        !   dstf               : differentiated source time function
        !   istf               : integrated source time function
        !   istftau            : integration of source time function multiplied with time
        !   coord              : coordinates of station relative to source
        !   nfflag             : turn on/off near field calculations
        !   ffflag             : turn on/off far field calculations
        !   mtflag             : turn on/off calculation of el. seismograms from moment tensor
        !   sfflag             : turn on/off calculation of el. seismograms from single force
        ! (out):
        !   all_elseism        : calculated elementary seismograms, first dim must be 36
        !                        second dim is the number of samples in each trace
        
        real, dimension(3)     :: gamma
        real, dimension(5)     :: matfac_mt, radpat, factor_mt
        real, dimension(3)     :: matfac_sf, factor_sf
        real                   :: r
        integer                :: n,p,q,i
       
        call make_direction_cosine( coord, gamma, r )
 
        call material_factors_sf( rho,alpha,beta, matfac_sf)
        call material_factors_mt( rho,alpha,beta, matfac_mt)
       
        if (mtflag) then
            i=1
            do n=1,3
                do p=1,3
                    do q=1,3
                        call radpat_mt( gamma,n,p,q,radpat )
                        call factors_mt( matfac_mt, radpat, r, factor_mt )
                        call elseis_mt( factor_mt, r, alpha, beta, toffset,dt,stf,dstf,istf,istftau, &
                                        nfflag,ffflag,all_elseism(i,:) )
                        i = i+1
                    end do
                end do
            end do
        end if
        
        if (sfflag) then
            i=28
            do n=1,3
                do p=1,3
                    call radpat_sf( gamma,n,p,radpat )
                    call factors_sf( matfac_sf, radpat, r, factor_sf )
                    call elseis_sf( factor_sf, r, alpha, beta, toffset,dt,stf,istf,istftau, &
                                    nfflag,ffflag,all_elseism(i,:) )
                    i = i+1
                end do
            end do
        end if
        
    end subroutine elseis_all
    
          
!********************************************************
!  ELEMENTARY SEISMOGRAMS: MOMENT TENSOR
!********************************************************

     pure subroutine elseis_mt( factors,                &
                          r, alpha, beta,              &
                          toffset,dt,                  &
                          stf,dstf,istf,istftau,       &
                          nfflag,ffflag,               &
                          elseism,addweight )

        real,    intent(in), dimension(5)     :: factors
        real,    intent(in)                   :: r, alpha, beta
        real,    intent(in)                   :: toffset, dt
        real,    intent(in), dimension(:)     :: stf,dstf,istf,istftau
        logical, intent(in)                   :: nfflag,ffflag
        real,    intent(inout), dimension(:)  :: elseism
        real,    intent(in), optional         :: addweight ! if present, add calculated seismogram to elseism, weighted by this

        integer ::    lstf,npt
        real    ::    term
        integer ::    ita,itb,it,ita_delta, itb_delta
        real    ::    t, ta, tb, ta_delta, tb_delta, integral_term
        
        npt = size(elseism)
        lstf = size(stf)
        
        ita_delta = nint(toffset/dt - r/alpha/dt)
        itb_delta = nint(toffset/dt - r/beta/dt)
        do it=1,npt
        
            t   = toffset + (it-1)*dt
            ta = t-r/alpha
            tb = t-r/beta
            ita = ita_delta + (it-1)
            itb = itb_delta + (it-1)
            call to_bounds( 0, lstf-1, ita )
            call to_bounds( 0, lstf-1, itb )
            
            if (nfflag) then  ! only needed for second order terms in near field 
                ta_delta = ta - ita*dt
                tb_delta = tb - itb*dt
            end if
            ita = ita+1
            itb = itb+1
            
            term = 0.0
            
            ! near field...
            if (nfflag) then
                ! this includes second order terms, because, when r is small ita=itb
                ! first order would be:
                !  integral_term = t * ( istf(ita)-istf(itb) ) &
                !         - ( istftau(ita)-istftau(itb) )
                ! this is also needed at the end, when (stf(lstf) /= 0)
                integral_term =   t * ( istf(ita)    - istf(itb) + ta_delta*stf(ita)    - tb_delta*stf(itb) )    &
                                  - (   istftau(ita) + ta_delta*stf(ita)*(ita-1)*dt  + 0.5*stf(ita)*ta_delta**2  &
                                      - istftau(itb) - tb_delta*stf(itb)*(itb-1)*dt  - 0.5*stf(itb)*tb_delta**2 )
                
                
                term = term + factors(1) * integral_term
                term = term + factors(2) * stf(ita)
                term = term + factors(3) * stf(itb)
                
            end if
           
            ! far field...
            if (ffflag) then
                term = term + factors(4) * dstf(ita)
                term = term + factors(5) * dstf(itb)
            end if
            
            if (present(addweight)) then
                elseism(it) = elseism(it) + term*addweight
            else
                elseism(it) = term
            end if
            
        end do

    end subroutine elseis_mt
      
      
!********************************************************
!  ELEMENTARY SEISMOGRAMS: SINGLE FORCES
!********************************************************
      
    pure subroutine elseis_sf( factors,                &
                          r, alpha, beta,              &
                          toffset,dt,                  &
                          stf,istf,istftau,            &
                          nfflag,ffflag,               &
                          elseism,addweight )

        real,    intent(in), dimension(3)     :: factors
        real,    intent(in)                   :: r, alpha, beta
        real,    intent(in)                   :: toffset, dt
        real,    intent(in), dimension(:)     :: stf,istf,istftau
        logical, intent(in)                   :: nfflag,ffflag
        real,    intent(inout), dimension(:)  :: elseism
        real,    intent(in), optional         :: addweight ! if present, add calculated seismogram to elseism, weighted by this

        integer ::    lstf,npt
        real    ::    term
        integer ::    ita,itb,it,ita_delta, itb_delta
        real    ::    t
        real    ::    ta, tb, ta_delta, tb_delta, integral_term

        npt = size(elseism)
        lstf = size(stf)
        
        ita_delta = nint(toffset/dt - r/alpha/dt)
        itb_delta = nint(toffset/dt - r/beta/dt)
        do it=1,npt
        
            t   = toffset + (it-1)*dt
            ta = t-r/alpha
            tb = t-r/beta
            ita = ita_delta + (it-1)
            itb = itb_delta + (it-1)
            call to_bounds( 0, lstf-1, ita )
            call to_bounds( 0, lstf-1, itb )
            
            if (nfflag) then  ! only needed for second order terms in near field 
                ta_delta = ta - ita*dt
                tb_delta = tb - itb*dt
            end if
            ita = ita+1
            itb = itb+1
           
            term = 0.0
            
            ! near field...
            if (nfflag) then
                ! this includes second order terms, because, when r is small ita=itb
                ! first order would be:
                !  integral_term = t * ( istf(ita)-istf(itb) ) &
                !         - ( istftau(ita)-istftau(itb) )
                integral_term =   t * ( istf(ita)    - istf(itb) + ta_delta*stf(ita)    - tb_delta*stf(itb) )    &
                                  - (   istftau(ita) + ta_delta*stf(ita)*(ita-1)*dt  + 0.5*stf(ita)*ta_delta**2  &
                                      - istftau(itb) - tb_delta*stf(itb)*(itb-1)*dt  - 0.5*stf(itb)*tb_delta**2 )
                
                term = term + factors(1)*integral_term
            end if
            
            ! far field...
            if (ffflag) then
                term = term + factors(2)*stf(ita)
                term = term + factors(3)*stf(itb)
            end if 

            if (present(addweight)) then
                elseism(it) = elseism(it) + term*addweight
            else
                elseism(it) = term
            end if
            
        end do
    
    end subroutine elseis_sf

!********************************************************


    pure subroutine factors_mt( matfac, radpat, r, factors )
        
        real, intent(in), dimension(5)  :: matfac, radpat
        real, intent(in)                :: r
        real, intent(out), dimension(5) :: factors
        
        factors(1) = matfac(1)*radpat(1)/r**4
        factors(2) = matfac(2)*radpat(2)/r**2
        factors(3) = matfac(3)*radpat(3)/r**2
        factors(4) = matfac(4)*radpat(4)/r
        factors(5) = matfac(5)*radpat(5)/r
        
    end subroutine factors_mt

    pure subroutine factors_sf( matfac, radpat, r, factors )
    
        real, intent(in), dimension(3)  :: matfac, radpat
        real, intent(in)                :: r
        real, intent(out), dimension(3) :: factors
        
        factors(1) = matfac(1)*radpat(1)/r**3
        factors(2) = matfac(2)*radpat(2)/r
        factors(3) = matfac(3)*radpat(3)/r
        
    end subroutine factors_sf

!********************************************************
      
    pure subroutine radpat_mt(gamma,n,p,q,rpc)
      
        real, intent(in)    :: gamma(3)
        integer, intent(in) :: n,p,q
        real, intent(out)   :: rpc(5)
      
        ! Calculation of radiation pattern for Green functions
        ! depending on Moment Tensor. 

        ! for near field terms
        ! is in range [-6,6]
        rpc(1)=(15*gamma(n)*gamma(p)*gamma(q)) &
                  -(3*gamma(n)*delta(p,q)) &
                  -(3*gamma(p)*delta(n,q)) &
                  -(3*gamma(q)*delta(n,p))

        ! for intermediate-field terms
        ! [-3,3]
        rpc(2)=(6*gamma(n)*gamma(p)*gamma(q)) &
                  -(gamma(n)*delta(p,q)) &
                  -(gamma(p)*delta(n,q)) &
                  -(gamma(q)*delta(n,p))
        
        ! [-2,2]
        rpc(3)=-((6*gamma(n)*gamma(p)*gamma(q)) &
                  -(gamma(n)*delta(p,q)) &
                  -(gamma(p)*delta(n,q)) &
                  -(2*gamma(q)*delta(n,p)))

        ! for far field terms
        ! [-1,1]
        rpc(4)=gamma(n)*gamma(p)*gamma(q)
        
        ! [-1,1]
        rpc(5)=-(gamma(n)*gamma(p)-delta(n,p))*gamma(q)
      
    end subroutine radpat_mt
            
!********************************************************
      
    pure subroutine radpat_sf(gamma,n,p,rpfc)
      
        real, intent(in)    :: gamma(3)
        integer, intent(in) :: n,p
        real, intent(out)   :: rpfc(3)

        ! calculation of radiation pattern for the Greens functions
        ! depending on single forces. 
        ! [-<2,2]
        rpfc(1)=3.0*gamma(n)*gamma(p)-delta(n,p)
        
        ! [-<1,1]
        rpfc(2)=gamma(n)*gamma(p)
        
        ! [-<1,1]
        rpfc(3)=delta(n,p)-gamma(n)*gamma(p)

    end subroutine radpat_sf
    
!********************************************************
      
    pure subroutine material_factors_mt( rho, alpha, beta, matfac )
        
        real, intent(in) :: rho, alpha, beta
        real, intent(out) :: matfac(5)
   
        ! calculates the material-dependent factors for the Greens function terms
        ! for moment tensors
   
        matfac(1) = 1.0/(4.0*PI*rho)
        matfac(2) = 1.0/(4.0*PI*rho*alpha**2)
        matfac(3) = 1.0/(4.0*PI*rho*beta**2)
        matfac(4) = 1.0/(4.0*PI*rho*alpha**3)
        matfac(5) = 1.0/(4.0*PI*rho*beta**3)
    
    end subroutine material_factors_mt

 !********************************************************      
     
    pure subroutine material_factors_sf( rho, alpha, beta, matfac )
        
        real, intent(in) :: rho, alpha, beta
        real, intent(out) :: matfac(3)
   
        ! calculates the material-dependent factors for the Greens function terms
        ! for single forces
        
        matfac(1) = 1.0/(4.0*PI*rho)
        matfac(2) = 1.0/(4.0*PI*rho*alpha**2)
        matfac(3) = 1.0/(4.0*PI*rho*beta**2)
        
    end subroutine material_factors_sf
    
!********************************************************
      
    pure subroutine make_direction_cosine( coord, gamma, r)
            
        real, intent(in), dimension(3)    :: coord
        real, intent(out), dimension(3)   :: gamma
        real, intent(out)                 :: r

        ! make direction cosines and r
        
        ! N,E,Z(down)=coord(1),coord(2),coord(3) 
        r=sqrt(coord(1)**2+coord(2)**2+coord(3)**2)

        ! gamma(i)=dr/dx(i)=cos(r,x(i))=coord(i)/r
        gamma(:)=coord(:)/r
        
    end subroutine make_direction_cosine
    
!********************************************************   
      
    pure subroutine make_istfs( dt, stf, istf, istftau )
    
        real, intent(in)                :: dt
        real, intent(in), dimension(:)  :: stf
        real, intent(out), dimension(:) :: istf, istftau
        
        ! calculate istf and istftau,
        ! used for calculation of near field terms in the elseis routines
        
        real, dimension(size(stf)) :: stftau
        integer :: i
        
        do i=1,size(stf)
            stftau(i) = stf(i) * (i-1)*dt
        end do
        call antiderivate( dt, stf, istf )
        call antiderivate( dt, stftau, istftau )
        
    end subroutine make_istfs
    
!********************************************************   
      
    pure subroutine to_bounds( lb, ub, i )

        integer, intent(in) :: lb, ub
        integer, intent(inout) :: i

        ! assure i is in [lb,ub] 

        if (i < lb) i = lb
        if (ub < i) i = ub

    end subroutine to_bounds

end module elseis
