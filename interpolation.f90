! 
! Copyright 2007 Francesco Pacchiani & Sebastian Heimann
! 

module interpolation

    use constants
    use util
    implicit none
    include 'fftw3.f'
    
    public gulunay2d, gulunay3d

  contains
    
    subroutine gulunay2d(A,Inter, ntmargin, nxmargin)
        
      ! Interpolate Green's function using the Generalized f-k interpolation method,
      ! also known as Gulunay's [2003] method.
        
        real, intent(inout), dimension(:,:)   :: A       ! (t,s)
        real, intent(out), dimension(:,:)     :: Inter   ! (t,kk)
        integer, intent(in)                   :: ntmargin, nxmargin
        
        integer                               :: t,s, kk, ff, fny
        integer                               :: x
        integer*8                             :: plan
        real                                  :: m
        integer                               :: l, il
        
        real,    dimension(size(A,1),                           size(Inter,2))  :: B
        real,    dimension(size(A,1)*(size(Inter,2)/size(A,2)), size(Inter,2))  :: C,D
        complex, dimension(size(A,1)/2+1,                       size(Inter,2))  :: fB, fInter, Operator   ! (t/2+1,kk)
        complex, dimension(size(C,1)/2+1,                       size(Inter,2))  :: fC, fD                 ! (ff/2+1,kk)
        
      ! check sizes
      
        if (size(Inter,1) /= size(A,1)) then
            call die( "gulunay2d(): time length of input and output arrays do not match")
        end if
        if (mod(size(Inter,2),size(A,2)) /= 0) then
            call die( "gulunay2d(): number of traces in ouput array is not a multiple of size of input array")
        end if
        
        l = size(Inter,2)/size(A,2)
        
        t=size(A,1)
        s=size(A,2)
        
        kk=size(Inter,2)
        ff=l*t
             
        ! --- Taper ---
        !
        do x=1,nxmargin/l
            A(:,x) = A(:,x) * (1. - cos(2.*pi*((x-1)/(2.*nxmargin/l)))) / 2.
        end do
        
        do x=s-nxmargin/l+1,s
            A(:,x) = A(:,x) * (1. - cos(2.*pi*((s-x)/(2.*nxmargin/l)))) / 2.
        end do
        
        do x=1,ntmargin/l
            A(x,:) = A(x,:) * (1. - cos(2.*pi*((x-1)/(2.*ntmargin/l)))) / 2.
        end do
        
        do x=t-ntmargin/l+1,t
            A(x,:) = A(x,:) * (1. - cos(2.*pi*((t-x)/(2.*ntmargin/l)))) / 2.
        end do
        
        !
        ! --- Insert zero traces in original data + 2-D FFT ---
        ! 
        B(:,1:kk:l) = A(:,:)
        do il=2,l
            B(:,il:kk:l) = 0.
        end do
        
        call sfftw_plan_dft_r2c_2d(plan,t,kk,B,fB,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)
        
        !
        ! --- Zero pad the original data in all dimensions + 2-D FFT ---
        !
        C(:,:) = 0.
        C(1:t,1:s) = A(:,:)
        
        call sfftw_plan_dft_r2c_2d(plan,ff,kk,C,fC,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)
        
        !
        ! --- Replace even traces in C with zeros + 2-D FFT ---
        !
        D(:,1:s:l) = C(:,1:s:l)
        do il=2,l
            D(:,il:s:l) = 0.
        end do
        
        call sfftw_plan_dft_r2c_2d(plan,ff,kk,D,fD,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)    
        
        !
        ! --- Add white noise to spectrum ---
        !
        fny = t/2+1
        m = 0.01*maxval( abs(fD(fny,:)) )
        where (abs(fD(1:fny,:)) .lt. m/1000.) 
            fD(1:fny,:) = cmplx(m,aimag(fD(1:fny,:)))
        end where
            
        where (abs(fD(1:fny,:)) .lt. m)
            fD(1:fny,:) = m/abs(fD(1:fny,:)) * fD(1:fny,:)
        end where
        
        !
        ! --- Calculate the interpolator ---
        !
        Operator(:,:) = fC(1:fny,:) / fD(1:fny,:)
        
        !
        ! --- Clip the iterpolation operator ---
        !        
        where (abs(Operator) .gt. l)
            Operator = l/abs(Operator) * Operator ! cmplx(l,aimag(Operator))
        end where
    
        where (abs(Operator) .lt. l*0.5)
            Operator = 0. ! cmplx(0.,aimag(Operator))
        end where
        
        !
        ! --- Compute interpolated traces ---
        !
        fInter(:,:) = fB(:,:) * Operator(:,:) / (t*kk)
        
        
        !
        ! --- Return in the time domain ---
        !
        call sfftw_plan_dft_c2r_2d(plan,t,kk,fInter,Inter,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)
            
    end subroutine
    
    subroutine gulunay3d(A,Inter, ntmargin, nxmargin, nzmargin)
       
        !  Interpolate Green's function using the Generalized f-k interpolation method,
        !  also known as Gulunay's [2003] method.
        
        real, intent(inout), dimension(:,:,:)   :: A       ! (t,sz,sx)
        real, intent(out), dimension(:,:,:)     :: Inter   ! (t,kkz,kkx)
        integer, intent(in)                   :: ntmargin, nxmargin, nzmargin
        
        integer                               :: t,sx,sz, kkx,kkz, ff, fny
        integer                               :: x
        integer*8                             :: plan
        integer                               :: l, il
        real                                  :: m, ls, lowcut    

        real,    dimension(size(A,1),                           size(Inter,2), size(Inter,3))  :: B
        real,    dimension(size(A,1)*(size(Inter,2)/size(A,2)), size(Inter,2), size(Inter,3))  :: C,D
        complex, dimension(size(A,1)/2+1,                       size(Inter,2), size(Inter,3))  :: fB, fInter, Operator   ! (t/2+1,kk)
        complex, dimension(size(C,1)/2+1,                       size(Inter,2), size(Inter,3))  :: fC, fD                 ! (ff/2+1,kk)
       
        if (size(Inter,1) /= size(A,1)) then
            call die( "gulunay3d(): time length of input and output arrays do not match")
        end if
        
        if (size(Inter,2)/size(A,2) /= size(Inter,3)/size(A,3)) then
            call die( "gulunay3d(): proportions of dimensions 2 and 3 between input and output arrays do not match")
        end if
        
        if (mod(size(Inter,2),size(A,2)) /= 0 .or. mod(size(Inter,3),size(A,3)) /= 0 ) then
            call die( "gulunay3d(): number of traces in ouput array is not a multiple of size of input array")
        end if
        
        l = size(Inter,2)/size(A,2)
        
        t=size(A,1)
        sz=size(A,2)
        sx=size(A,3)
           
        kkx=l*sx
        kkz=l*sz
        
        ff=l*t
                
        !
        ! --- Taper ---
        !
        do x=1,nxmargin/l
            A(:,:,x) = A(:,:,x) * (1. - cos(2.*pi*((x-1)/(2.*nxmargin/l)))) / 2.
        end do
        
        do x=sx-nxmargin/l+1,sx
            A(:,:,x) = A(:,:,x) * (1. - cos(2.*pi*((sx-x)/(2.*nxmargin/l)))) / 2.
        end do
        
        do x=1,nzmargin/l
            A(:,x,:) = A(:,x,:) * (1. - cos(2.*pi*((x-1)/(2.*nzmargin/l)))) / 2.
        end do
        
        do x=sz-nzmargin/l+1,sz
            A(:,x,:) = A(:,x,:) * (1. - cos(2.*pi*((sz-x)/(2.*nzmargin/l)))) / 2.
        end do
        
        do x=1,ntmargin/l
            A(x,:,:) = A(x,:,:) * (1. - cos(2.*pi*((x-1)/(2.*ntmargin/l)))) / 2.
        end do
        
        do x=t-ntmargin/l+1,t
            A(x,:,:) = A(x,:,:) * (1. - cos(2.*pi*((t-x)/(2.*ntmargin/l)))) / 2.
        end do
        
        !
        ! --- Insert zero traces in original data + 2-D FFT ---
        !
        B(:,1:kkz:l,1:kkx:l) = A(:,:,:)
        do il=2,l
            B(:,il:kkz:l,il:kkx:l) = 0.
        end do
        
        call sfftw_plan_dft_r2c_3d(plan,t,kkz,kkx,B,fB,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)
        
        !
        ! --- Zero pad the original data in all dimensions + 2-D FFT ---
        !
        C(:,:,:) = 0.
        C(1:t,1:sz,1:sx) = A(:,:,:)
        
        call sfftw_plan_dft_r2c_3d(plan,ff,kkz,kkx,C,fC,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)
        
        !
        ! --- Replace even traces in C with zeros + 2-D FFT ---
        !
        D(:,1:sz:l,1:sx:l) = C(:,1:sz:l,1:sx:l)
        do il=2,l
            D(:,il:sz:l,il:sx:l) = 0.
        end do
        
        call sfftw_plan_dft_r2c_3d(plan,ff,kkz,kkx,D,fD,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)    
        
        !
        ! --- Add white noise to spectrum ---
        !
        fny = t/2+1
        m = 0.01*maxval( abs(fD(fny,:,:)) )
        where (abs(fD(1:fny,:,:)) .lt. m/1000.) 
            fD(1:fny,:,:) = cmplx(m,aimag(fD(1:fny,:,:)))
        end where
            
        where (abs(fD(1:fny,:,:)) .lt. m)
            fD(1:fny,:,:) = m/abs(fD(1:fny,:,:)) * fD(1:fny,:,:)
        end where
        
        !
        ! --- Calculate the interpolator ---
        !
        Operator(:,:,:) = fC(1:fny,:,:) / fD(1:fny,:,:)
            
        !
        ! --- Clip the iterpolation operator ---
        !
        
        lowcut = 0.5*l**2
        ls = l**2
        
        where (abs(Operator) .gt. ls)
            Operator = ls/abs(Operator) * Operator
        end where
        
        where (abs(Operator) .lt. lowcut)
            Operator = 0.
        end where        
        
        !
        ! --- Compute interpolated traces ---
        !
        fInter(:,:,:) = fB(:,:,:) * Operator(:,:,:) / (t*kkx*kkz)
            
        !
        ! --- Return in the time domain ---
        !
        call sfftw_plan_dft_c2r_3d(plan,t,kkz,kkx,fInter,Inter,FFTW_ESTIMATE)
        call sfftw_execute(plan)
        call sfftw_destroy_plan(plan)
            
    end subroutine   
   
end module
