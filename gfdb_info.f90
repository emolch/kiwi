! $Id: gfdb_info.f90 702 2008-03-31 07:07:40Z sebastian $ 
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


program gfdb_info

  ! This program is used to get some information about a Greens function database 
  ! built with gfdb_build.
  ! 
  ! usage: gfdb_info database
  !
  ! Complete documentation is available on
  ! 
  !   http://kinherd.org/power/trac/wiki/GFDBInfoTool
  !
  
    use util
    use gfdb
    use better_varying_string
    use varying_string_getarg

    ! use f90_unix_env

    implicit none    

    type(varying_string)  :: basefn, command
    type(t_gfdb)          :: db
    integer               :: indexmemory_per_chunk
    integer, allocatable, dimension(:,:) :: ntraces_in_chunks
    integer :: i, ntraces, ntraces_used
    character, parameter :: eol = char(10)

    
    g_pn = 'gfdb_info'
    g_usage = 'usage: '//g_pn//' database' // eol // eol // &
              "documentation:  " // &
              "http://kinherd.org/power/trac/wiki/GFDBInfoTool"

    if (iargc() /= 1 .and. iargc() /= 2) call usage()

    call vs_getarg( 1, basefn )

    call gfdb_init(db, basefn)
    
    if (iargc() == 2) then
        call vs_getarg( 2, command )
        if (command == 'infomap') then
            call gfdb_dump_infomap( db, stdout )
        else if (command == 'contents') then
            call gfdb_dump_contents( db, stdout )
        else if (command == 'missing') then
            call gfdb_dump_missing( db, stdout )
        else if (command == 'stats') then
            call gfdb_dump_stats( db, stdout )
        else if (command == 'chunkstats') then
            call gfdb_infos( db, indexmemory_per_chunk, ntraces_in_chunks )
            write (*,'(a)') char("indexmemory_per_chunk="// indexmemory_per_chunk)
            ntraces = 0
            ntraces_used = 0
            do i=1,size(ntraces_in_chunks,2)
                write (*,'(a)') char("chunk"//i//"_traces="// ntraces_in_chunks(1,i) // "/" //&
                                                        ntraces_in_chunks(2,i) )
                ntraces_used = ntraces_used+ntraces_in_chunks(1,i)
                ntraces = ntraces+ntraces_in_chunks(2,i)
            end do
            write (*,'(a)') char("total_traces="// ntraces_used // "/" //&
                                                        ntraces )
        end if
    else
    
        write (*,'(a)') char("dt="// db%dt)
        write (*,'(a)') char("dx="// db%dx)
        write (*,'(a)') char("dz="// db%dz)
        write (*,'(a)') char("firstx="// db%firstx)
        write (*,'(a)') char("firstz="// db%firstz)
        write (*,'(a)') char("nx="// db%nx)
        write (*,'(a)') char("nz="// db%nz)
        write (*,'(a)') char("ng="// db%ng)
        write (*,'(a)') char("nchunks="// db%nchunks)
        write (*,'(a)') char("nxc="// db%nxc)
        
    end if
    
    call gfdb_destroy( db )

    call cleanup()

end program
