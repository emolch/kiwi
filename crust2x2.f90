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

module crust2x2

  ! This module allows access to the the global crustal model at 2x2 Degrees by 
  ! Gabi Laske, Guy Masters and Christine Reif.
  !
  ! The model is available at 
  !    http://mahi.ucsd.edu/Gabi/rem.dir/crust/crust2.html
  
    
    use better_varying_string
    use util
    use unit
    use orthodrome
    
    implicit none
    
    private
    
    public crust2x2_load
    public crust2x2_get_profile
    public crust2x2_get_profile_averages
    public crust2x2_describe_profile
    public crust2x2_get_at_depth

    integer, parameter :: nlayers = 7
    integer, parameter :: ntypes = 360
    integer, parameter :: nla  = 90
    integer, parameter :: nlo = 180
    
    type, public :: t_crust2x2_1d_profile
        character(len=2)   :: id
        real, dimension(nlayers+1) :: vp, vs, rho
        real, dimension(nlayers)   :: thickness
        real                       :: elevation
    end type
    
    character(len=*), parameter :: fn_keys =      'CNtype2_key.txt'
    character(len=*), parameter :: fn_map =       'CNtype2.txt'
    character(len=*), parameter :: fn_elevation = 'CNelevatio2.txt'
    
    character(len=12) :: names(nlayers)
    logical :: found
    
    data names/'water','ice','soft sed.','hard sed.','upper crust','middle crust','lower crust'/
    
    integer, parameter :: LWATER = 1
    integer, parameter :: LICE = 2
    integer, parameter :: LSOFTSED = 3
    integer, parameter :: LHARDSED = 4
    integer, parameter :: LUPPERCRUST = 5
    integer, parameter :: LMIDDLECRUST = 6
    integer, parameter :: LLOWERCRUST = 7
    integer, parameter :: LBELOWCRUST = 8
    
    logical, public :: crust2x2_loaded = .false.
    type(t_crust2x2_1d_profile), dimension(nlo,nla) :: model
    
  contains
  
    subroutine crust2x2_load( directory, ok )
    
      ! Load and initialize the crust2x2 model. 
      !
      ! It looks into `directory` for the model files.
    
        character(len=*) :: directory
        logical, intent(out) :: ok
        
        call load_crustal_model( directory,  model, ok )
        if (ok) crust2x2_loaded = .true.
        
    end subroutine
    
    subroutine crust2x2_get_profile( location, profile )
    
      ! Get crustal profile at a specific location.
      
        type(t_geo_coords), intent(in)          :: location
        type(t_crust2x2_1d_profile), intent(out) :: profile
        integer :: ilat, ilon
        
        if (.not. crust2x2_loaded) then
            call die("crust2x2 model not loaded")
        end if 
        
        call latlon2indices( real(location%lat), real(location%lon), ilat, ilon )
        profile = model(ilon,ilat)
       
    end subroutine
    
    subroutine crust2x2_describe_profile( profile, u )
    
      ! Print information about a specific crustal profile.
      !
      ! The information is printed to unit `u`.
    
        type(t_crust2x2_1d_profile), intent(in) :: profile
        integer, intent(in) :: u
        real :: vvp, vvs, vrho, vthi
        integer :: i
        
        call crust2x2_get_profile_averages( profile, vvp, vvs, vrho, vthi )
            
        write (u,"(a, a, g15.5)") 'type, elevation: ', profile%id, profile%elevation
        write (u,"(a, 4g15.5)") 'crustal thickness, ave. vp, vs, rho:', vthi,vvp,vvs,vrho                    
        write (u,"(a, a15, 3g15.5)") 'Mantle below Moho: ave. vp, vs, rho:', ' ',&
                profile%vp(LBELOWCRUST), &
                profile%vs(LBELOWCRUST), &
                profile%rho(LBELOWCRUST)
    
    
        write (u,*)
        write (u,'(i3,a)') nlayers,'-layer crustal profile (thickness, vp,vs,rho)' 
        
        do i=1,nlayers
            write (u,"(4g15.5,a)") profile%thickness(i), profile%vp(i), &
                profile%vs(i),profile%rho(i),names(i)
        end do
    
    end subroutine
    
    subroutine crust2x2_get_profile_averages( profile, vvp, vvs, vrho, vthi )
    
      ! Get crustal averages for vp, vs and density and total crustal thickness,
      ! for a given profile.
      !
      ! Takes into account ice layer.
      ! Does not take into account water layer.
      
        type(t_crust2x2_1d_profile), intent(in) :: profile
        real, intent(out) :: vvp, vvs, vrho, vthi
        integer :: i
        
        vthi=0.
        vvp=0.
        vvs=0.
        vrho=0.
        
        do i=2,nlayers
            vthi = vthi + profile%thickness(i)
            vvp  = vvp  + profile%thickness(i) / profile%vp(i)
            vvs  = vvs  + profile%thickness(i) / profile%vs(i)
            vrho = vrho + profile%thickness(i) * profile%rho(i)
        end do
        
        vvp  = vthi / vvp
        vvs  = vthi / vvs
        vrho = vrho / vthi
        
    end subroutine
    
    subroutine crust2x2_get_at_depth( profile, depth, vp, vs, rho )
    
      ! Retrieve vp, vs and density at given depth in profile.
    
        type(t_crust2x2_1d_profile), intent(in) :: profile
        real, intent(in)  :: depth
        real, intent(out) :: vp, vs, rho
        integer :: i
        real :: d
        
        d = 0
        do i=3,nlayers
            d = d + profile%thickness(i)
            if (d >= depth) then
                vp = profile%vp(i)
                vs = profile%vs(i)
                rho = profile%rho(i)
                return
            end if
        end do
        
        vp = profile%vp(LBELOWCRUST)
        vs = profile%vs(LBELOWCRUST)
        rho = profile%rho(LBELOWCRUST)
        
    end subroutine
    
  ! internal stuff below here
    
    subroutine latlon2indices( lat, lon, ilat, ilon )
        real, intent(in) :: lat, lon
        integer, intent(out) :: ilat, ilon
        
        real :: flat, flon, cola, dx
        
        flat = lat
        flon = lon
        flat = clip( flat,-90.,90. )
        flon = wrap( flon,-180.,180. )
        
        dx = 360./nlo
        cola=90.-flat
        ilat=int(cola/dx)+1
        ilon=int((flon+180.)/dx)+1
        
    end subroutine
  
    subroutine swap(a,b)
        real, intent(inout) :: a,b
        real :: aux
        aux = a
        a = b
        b = aux
    end subroutine
    
    elemental function wrap( x, mi, ma )
        real, intent(in) :: x, mi, ma
        real :: wrap
        wrap = x
        
        if (mi .le. x .and. x .le. ma) return
        wrap = x - floor((x-mi)/(ma-mi)) * (ma-mi)
    end function
    
    elemental function clip( x, mi, ma )
        real, intent(in) :: x, mi, ma
        real :: clip

        clip = min(max(mi,x),ma)
        
    end function
    
    subroutine load_crustal_model( directory,  amap, ok  )
        character(len=*), intent(in) :: directory
        type(t_crust2x2_1d_profile), dimension(nlo,nla), intent(inout) :: amap
        logical, intent(out) :: ok
        
        type(t_crust2x2_1d_profile), dimension(ntypes) :: ctypes
        
        character(len=2) :: ctype_ids(nlo)
        integer :: i,j,l, ilat
        integer :: u, iostat
        
        call claim_unit( u )
        
      ! read in key crust types
      
      ! skip header
        ok = .false.
        open(u,file=directory//'/'//fn_keys, status='old', iostat=iostat)
        if (iostat /= 0) then
            call error( "can't open file: "// directory//'/'//fn_keys )
            return
        end if
        
        do i=1,5
            read(u,*)
        end do
        
        do i=1,ntypes
            read(u,*) ctypes(i)%id
            read(u,*) ( ctypes(i)%vp(l), l=1,nlayers+1 )
            read(u,*) ( ctypes(i)%vs(l), l=1,nlayers+1 )
            read(u,*) ( ctypes(i)%rho(l), l=1,nlayers+1 )
            read(u,*) ( ctypes(i)%thickness(l), l=1,nlayers )
            
          ! want SI units
            do l=1,nlayers+1
                ctypes(i)%vp(l) = ctypes(i)%vp(l)*1000.
                ctypes(i)%vs(l) = ctypes(i)%vs(l)*1000.
                ctypes(i)%rho(l) = ctypes(i)%rho(l)*1000.
            end do
            do l=1,nlayers
                ctypes(i)%thickness(l) = ctypes(i)%thickness(l)*1000.
            end do
            
          ! flip ice and water layers
            call swap( ctypes(i)%vp(1), ctypes(i)%vp(2) )
            call swap( ctypes(i)%vs(1), ctypes(i)%vs(2) )
            call swap( ctypes(i)%rho(1), ctypes(i)%rho(2) )
            call swap( ctypes(i)%thickness(1), ctypes(i)%thickness(2) )
        end do
        close(u) 
        
      ! read CNtype file
        open(u,file=directory//'/'//fn_map, status='old', iostat=iostat)
        if (iostat /= 0) then
            call error( "can't open file: "// directory//'/'//fn_map )
            return
        end if
        read (u,*)

        do j=1,nla
            read(u,*) ilat, ctype_ids
            do i=1,nlo
                found = .false.
                type_loop: do l=1,ntypes
                    if(ctype_ids(i).eq.ctypes(l)%id )then
                        amap(i,j) = ctypes(l)
                        found = .true.
                        exit type_loop
                    endif
                end do type_loop
                if (.not. found) then
                    call die('crust type code not found: ' // ctype_ids(i) //  &
                              ' latitude: ' // ilat // ' long index: ' // i )
                end if
            end do
        end do
        close(u)
        
        open(u,file=directory//'/'//fn_elevation, status='old', iostat=iostat)
        if (iostat /= 0) then
            call error( "can't open file: "// directory//'/'//fn_elevation )
            return
        end if
        
        read (u,*)
        do j=1,nla
            read(u,*) ilat, (amap(i,j)%elevation,i=1,nlo) 
            
            ! replace water thickness with more accurate value from evelation map
            do i=1,nlo
                if (amap(i,j)%elevation .lt. 0. .and. amap(i,j)%thickness(LWATER) .ne. 0.) then
                    amap(i,j)%thickness(LWATER) = - amap(i,j)%elevation
                end if
            end do
            
        end do
        close(u)
        
        call release_unit( u )
        ok = .true.
        
    end subroutine
    
    
end module