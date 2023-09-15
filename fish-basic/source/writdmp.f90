!=======================================================================
!
!     writdmp
!
!=======================================================================

      subroutine writdmp(num,time,s,v,nx,ny,nz)

      use comm_module
!cscs: iflport for ifort compiler
    use iflport, dummy1=>time,dummy2=>date !required for ifort

      implicit none

      integer, intent(in) :: num
      type(comm_wld), intent(in) :: nx,ny,nz
      real, intent(in) :: time
      real(4), dimension(ns,nx%l,ny%l,nz%l), intent(in) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l), intent(in) :: v

!-----------------------------------------------------------------------

      logical :: found
      character(8) :: date
      character(11) :: suffix
      character(24) :: host
      integer :: i,j,k,m,status

      call prof_enter(19,'writdmp')

!.....dump general information..........................................
      write(suffix,11) num,mr
11    format('_',i5,'_',i4)
      do m=1,11
        if (suffix(m:m).eq.' ') suffix(m:m) = '0'
      enddo
      if (mr.eq.0) write(6,*) 'Info: writing '//trim(name)//suffix(1:6)
      open(1,file=trim(glbpath)//'proc'//suffix(8:11)//'/' &
     &  //trim(name)//suffix//'.dmp', &
     &  status='unknown',form='unformatted')
      write(1) date
      write(1) time,dx

!.....dump data size and content........................................
      write(1) nv,ns-ndbg,nx,ny,nz
      write(1) v(:,:,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      write(1) s(1:ns-ndbg,nx%m:nx%n,ny%m:ny%n,nz%m:nz%n)
      close(1)

!.....record dump in logfile............................................
      write(suffix,22) mr
22    format('_',i4)
      do m=1,5
        if (suffix(m:m).eq.' ') suffix(m:m) = '0'
      enddo
      open(1,file=trim(txtpath)//trim(name)//suffix(1:5)//'.dmp', &
     &  position='append',status='unknown')
!cscs: different non-standard commands for ifort and pathscale
!     status = hostnam(host)	!ifort
      status = hostnm(host)	!pathscale
      host = adjustr(host)
      call date_and_time(date)
      write(1,33) num,time,host,date
33    format(i12,g12.4,a24,a12)
      close(1)
      if (mr.eq.0) call write_restart_loc(num)

      call prof_exit(19)

      end subroutine writdmp

!=======================================================================
