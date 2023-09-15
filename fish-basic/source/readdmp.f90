!=======================================================================
!
!     readdmp
!
!=======================================================================

      subroutine readdmp(iter,time,s,v,nx,ny,nz)

      use comm_module
!cscs: iflport for ifort compiler
    use iflport, dummy1=>time,dummy2=>date !required for ifort

      implicit none

      integer, intent(out) :: iter
      type(comm_wld), intent(in) :: nx,ny,nz
      real, intent(out) :: time
      real(4), dimension(ns,nx%l,ny%l,nz%l), intent(out) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l), intent(out) :: v

!-----------------------------------------------------------------------

      logical :: found
      character(8) :: date
      character(11) :: suffix
      character(12) :: num
      character(24) :: host
      character(57) :: path
      character(72) :: line
      integer(4) :: result
      integer :: i,j,k,m,status,nvf,nsf,i1,i2,i3,i4,i5,i6,i7,i8
      type(comm_wld) :: nxf,nyf,nzf
      real :: dxf

      call prof_enter(20,'readdmp')

!.....prepare suffix for filename of logfile............................
      write(suffix,22) mr
22    format('_',i4)
      do m=1,5
        if (suffix(m:m).eq.' ') suffix(m:m) = '0'
      enddo

!-----start from scratch------------------------------------------------
      if (restart.eq.0) then

!.....create directories................................................
        result = system('mkdir -p '//trim(glbpath)// &
     &    'proc'//suffix(2:5)//'/')
        result = system('mkdir -p '//trim(binpath)// &
     &    'proc'//suffix(2:5)//'/')

!.....is there an obsolete logfile?.....................................
        inquire(file=trim(txtpath)//'delete',exist=found)
        if (found) inquire(file=trim(txtpath)//trim(name)//suffix(1:5) &
     &    //'.dmp',exist=found)
        if (found) then
          open(1,file=trim(txtpath)//trim(name)//suffix(1:5)//'.dmp', &
     &      status='old')

!.....go through list in logfile........................................
          do
            read(1,66,end=100) line
66          format(a72)
            host = adjustl(line(25:48))
            path = trim(glbpath)//'proc'//suffix(2:5)//'/'

!.....if the restart file still exists, delete it!.....................
            write(suffix,77) line(8:12),mr
77          format('_',a5,'_',i4)
            do m=1,11
              if (suffix(m:m).eq.' ') suffix(m:m) = '0'
            enddo
            inquire(file=trim(path)//trim(name)//suffix//'.dmp', &
     &        exist=found)
            if (found) then
              if (mr.eq.0) write(6,*) &
     &          'Info: deleting ',trim(name)//suffix(1:6)//'_xxxx.dmp'
              open(2,file=trim(path)//trim(name)//suffix//'.dmp', &
     &          status='old',form='unformatted')
              close(2,status='delete')
            endif
          enddo
100       continue
          close(1,status='delete')
        endif

!.....now generate initial state for new run............................
        call init(s,v,nx,ny,nz)
        time = 0.
        call output(0,time,s,v,nx,ny,nz)

!-----read in from existing file----------------------------------------
      else

!.....construct path to restart file for this process...................
        path = trim(glbpath)//'proc'//suffix(2:5)//'/'

!.....read general information..........................................
        write(suffix,11) restart,mr
11      format('_',i5,'_',i4)
        do m=1,11
          if (suffix(m:m).eq.' ') suffix(m:m) = '0'
        enddo
        if (mr.eq.0) write(6,*) 'Info: reading ',trim(name)//suffix(1:6)
        open(1,file=trim(path)//trim(name)//suffix//'.dmp', &
     &    status='old',form='unformatted')
        read(1) date
        read(1) time,dxf
        read(1) nvf,nsf,nxf,nyf,nzf
        if (mr.eq.0) then
          if (dxf.ne.dx) &
     &      write(6,*) 'Warning: mismatch between dx and dxf: ',dx,dxf
          if (nvf.ne.nv) &
     &      write(6,*) 'Warning: mismatch between nv and nvf: ',nv,nvf
          if (nsf.ne.ns-ndbg) &
     &      write(6,*) 'Warning: mismatch between ns and nsf: ', &
     &        ns-ndbg,nsf
          if ((nx%g.ne.nxf%g).or.(ny%g.ne.nyf%g).or.(nz%g.ne.nzf%g)) &
     &      write(6,*) 'Warning: mismatch in resolution: (', &
     &        nx%g,ny%g,nz%g,') (',nxf%g,nyf%g,nzf%g,')'
          if ((nx%n-nx%m.ne.nxf%n-nxf%m).or. &
     &        (ny%n-ny%m.ne.nyf%n-nyf%m).or. &
     &        (nz%n-nz%m.ne.nzf%n-nzf%m)) &
     &      write(6,*) 'Warning: mismatch in decomposition: (', &
     &        nx%n-nx%m+1,ny%n-ny%m+1,nz%n-nz%m+1,') (', &
     &        nxf%n-nxf%m+1,nyf%n-nyf%m+1,nzf%n-nzf%m+1,')'
        endif

!.....read data.........................................................
        read(1) v(:,:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
        read(1) s(1:ns-ndbg,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
        close(1)
      endif
      iter = restart*skip

!.....make title entry in logfile.......................................
      write(suffix,22) mr
      do m=1,5
        if (suffix(m:m).eq.' ') suffix(m:m) = '0'
      enddo
      open(1,file=trim(txtpath)//trim(name)//suffix(1:5)//'.dmp', &
     &  position='append',status='unknown')
      write(1,44)
      write(1,55) 'File','Time','Host','Date'
      write(1,44)
44    format(60('-'))
55    format(2a12,a24,a12)
      close(1)

      call prof_exit(20)

      end subroutine readdmp

!=======================================================================
