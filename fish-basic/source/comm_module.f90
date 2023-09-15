!***********************************************************************
!
!     comm_module
!
!***********************************************************************

      module comm_module

      use mhd_parameter_module

      implicit none

!.....new type for dimensions, work load distribution, comm handles.....
      type comm_wld
        integer :: d !global label of direction {x=1,y=2,z=3}
        integer :: g !global number of zones in that direction
        integer :: r !global index of index 0 in local array
        integer :: m,n !start and end of local section without buffers
        integer :: l !dimension of array inclusive buffers
        integer, dimension(4) :: requests !communication handles
      end type comm_wld

!-----------------------------------------------------------------------
!
!     Indexing conventions, for example in x-direction:
!
!     The total computational domain is split into nx%g zones. nx%g is
!     set in this module to be nx%g=nxg (nxg is specified in the
!     parameter_module.f).
!
!     The number of available processes is compared to the total number
!     of zones in three dimensions. The algorithm below tries to divide
!     the zones into equally sized rectangular subdomains, each assigned
!     to one process. Each subdomain has nx%n-nx%m+1 zones in x-direction
!     and carries a halo of buffer zones that mirror the boundary of
!     neighboring subdomains. The length of the subdomain including the
!     buffer is nx%l:
!
!     1   nx%m             nx%n  nx%l (local indices in process 1)
!     v     v                v     v
!     buffer    subdomain1    buffer
!     ......------------------......
!                       ......------------------......
!                       buffer    subdomain2    buffer
!                       ^     ^                ^     ^
!                       1   nx%m             nx%n  nx%l (in process 2)
!         
!     Usually local indices running from 1 to nx%l are used to describe
!     the subdomain. The global index is obtained from the local one
!     by addition of the (local) reference nx%r. The y- and z-direction
!     are described in the same way. Thus, a typical loop assigning to
!     each (non-buffer) zone in a subdomain the value of a function that
!     depends on the global position with respect to an origin (rx,ry,rz)
!     would read:
!
!     do k=nz%m,nz%n
!       do j=ny%m,ny%n
!         do i=nx%m,nx%n
!           field(i,j,k)=global_function(nx%r+i-rx,ny%r+j-ry,nz%r+k-rz)
!         enddo
!       enddo
!     enddo
!
!-----------------------------------------------------------------------

!.....standard variables................................................
      include 'mpif.h'
      integer :: ncpu,mr,ierr
      integer, parameter :: nr=4 !number of requests per direction

!.....message tags......................................................
      integer :: xtagdw=1,xtagup=2,ytagdw=3,ytagup=4,ztagdw=5,ztagup=6

!.....storage for communication buffers.................................
      real(4), dimension((ns+3*nv)*buf &
     &  *(2*yzbuf+nmax)*(2*yzbuf+nmax)) :: send0,send1,recv0,recv1

      contains

!=======================================================================
!
!     comm_mpinit
!
!=======================================================================

      subroutine comm_mpinit(nx,ny,nz)
      
      implicit none    

      type(comm_wld), intent(out) :: nx,ny,nz

!-----------------------------------------------------------------------
!
!     - initialize message passing interface
!     - set workload distribution parameters for type(comm_wld)
!     - open half-channels for persistent communication between cubes
!
!-----------------------------------------------------------------------

      integer :: npercpu,nxmn,nymn,nzmn,i,j,k,imax,jmax,kmax
      integer :: mr0,mr1

!.....initialization for message passing................................
      call mpi_init(ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_comm_size(mpi_comm_world, ncpu, ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_comm_rank(mpi_comm_world, mr, ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      if (mr.eq.0) then
        write(6,*) 'MPI mhd code initialized'
        write(6,22) nxg,nyg,nzg,ncpu
22      format(' Cube of sides = ',i4,'*',i4,'*',i4,' over ', &
     &    i6,' processes')
      endif

!.....determine decomposition in each direction.........................
      call decompose(ncpu,nxg,nyg,nzg,nxmn,nymn,nzmn)
      if (nxmn.le.0) then
        write(6,11) nxg,nyg,nzg,ncpu
11      format('Error: ',i5,'*',i5,'*',i5, &
     &    ' does not fit on ',i8,' cpus!')
        call comm_mpistop('mpinit')
      endif
        
!.....3D index of mr (max is added to avoid negative indices later).....
      imax = nxg/nxmn
      jmax = nyg/nymn
      kmax = nzg/nzmn
      i = imax+mod(mr,imax)+1
      j = jmax+mod(mr/imax,jmax)+1
      k = kmax+mod(mr/jmax/imax,kmax)+1

!.....assign comm_wld for x-direction...................................
      nx%d = 1
      nx%g = nxg
      nx%r = mod(i-1,imax)*nxmn-yzbuf
      nx%m = yzbuf+1
      nx%n = yzbuf+nxmn
      nx%l = 2*yzbuf+nxmn
      if (mr.eq.0) then
        write(*,*) 'X axis broken into ',nx%n-nx%m+1,'cell-thick slabs'
        if (nx%n-nx%m+1.gt.nmax) then
          write(6,*) 'Error: parameter nmax=',nmax,'is too small!'
          call comm_mpistop('mpinit')
        endif
      endif

!.....assign comm_wld for y-direction...................................
      ny%d = 2
      ny%g = nyg
      ny%r = mod(j-1,jmax)*nymn-yzbuf
      ny%m = yzbuf+1
      ny%n = yzbuf+nymn
      ny%l = 2*yzbuf+nymn
      if (mr.eq.0) then
        write(*,*) 'Y axis broken into ',ny%n-ny%m+1,'cell-thick slabs'
        if (ny%n-ny%m+1.gt.nmax) then
          write(6,*) 'Error: parameter nmax=',nmax,'is too small!'
          call comm_mpistop('mpinit')
        endif
      endif
 
!.....assign comm_wld for z-direction...................................
      nz%d = 3
      nz%g = nzg
      nz%r = mod(k-1,kmax)*nzmn-yzbuf
      nz%m = yzbuf+1
      nz%n = yzbuf+nzmn
      nz%l = 2*yzbuf+nzmn
      if (mr.eq.0) then
        write(*,*) 'Z axis broken into ',nz%n-nz%m+1,'cell-thick slabs'
        if (nz%n-nz%m+1.gt.nmax) then
          write(6,*) 'Error: parameter nmax=',nmax,'is too small!'
          call comm_mpistop('mpinit')
        endif
      endif

!.....install persistent communications in x-direction..................
      mr0 = (mod(k-1,kmax)*jmax + mod(j-1,jmax))*imax + mod(i-1-1,imax)
      mr1 = (mod(k-1,kmax)*jmax + mod(j-1,jmax))*imax + mod(i-1+1,imax)
      call mpi_send_init(send0,(ns+3*nv)*buf*ny%l*nz%l, &
     &  MPI_REAL,mr0,xtagdw,MPI_COMM_WORLD, &
     &  nx%requests(1),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_send_init(send1,(ns+3*nv)*buf*ny%l*nz%l, &
     &  MPI_REAL,mr1,xtagup,MPI_COMM_WORLD, &
     &  nx%requests(2),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_recv_init(recv1,(ns+3*nv)*buf*ny%l*nz%l, &
     &  MPI_REAL,mr1,xtagdw,MPI_COMM_WORLD, &
     &  nx%requests(3),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_recv_init(recv0,(ns+3*nv)*buf*ny%l*nz%l, &
     &  MPI_REAL,mr0,xtagup,MPI_COMM_WORLD, &
     &  nx%requests(4),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')

!.....install persistent communications in y-direction..................
      mr0 = (mod(k-1,kmax)*jmax + mod(j-1-1,jmax))*imax + mod(i-1,imax)
      mr1 = (mod(k-1,kmax)*jmax + mod(j-1+1,jmax))*imax + mod(i-1,imax)
      call mpi_send_init(send0,(ns+3*nv)*buf*nz%l*nx%l, &
     &  MPI_REAL,mr0,ytagdw,MPI_COMM_WORLD, &
     &  ny%requests(1),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_send_init(send1,(ns+3*nv)*buf*nz%l*nx%l, &
     &  MPI_REAL,mr1,ytagup,MPI_COMM_WORLD, &
     &  ny%requests(2),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_recv_init(recv1,(ns+3*nv)*buf*nz%l*nx%l, &
     &  MPI_REAL,mr1,ytagdw,MPI_COMM_WORLD, &
     &  ny%requests(3),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_recv_init(recv0,(ns+3*nv)*buf*nz%l*nx%l, &
     &  MPI_REAL,mr0,ytagup,MPI_COMM_WORLD, &
     &  ny%requests(4),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')

!.....install persistent communications in z-direction..................
      mr0 = (mod(k-1-1,kmax)*jmax + mod(j-1,jmax))*imax + mod(i-1,imax)
      mr1 = (mod(k-1+1,kmax)*jmax + mod(j-1,jmax))*imax + mod(i-1,imax)
      call mpi_send_init(send0,(ns+3*nv)*buf*nx%l*ny%l, &
     &  MPI_REAL,mr0,ztagdw,MPI_COMM_WORLD, &
     &  nz%requests(1),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_send_init(send1,(ns+3*nv)*buf*nx%l*ny%l, &
     &  MPI_REAL,mr1,ztagup,MPI_COMM_WORLD, &
     &  nz%requests(2),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_recv_init(recv1,(ns+3*nv)*buf*nx%l*ny%l, &
     &  MPI_REAL,mr1,ztagdw,MPI_COMM_WORLD, &
     &  nz%requests(3),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')
      call mpi_recv_init(recv0,(ns+3*nv)*buf*nx%l*ny%l, &
     &  MPI_REAL,mr0,ztagup,MPI_COMM_WORLD, &
     &  nz%requests(4),ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpinit')

      end subroutine comm_mpinit

!=======================================================================
!
!     comm_mpistartall
!
!=======================================================================

      subroutine comm_mpistartall(requests)

      integer, dimension(nr), intent(in) :: requests

      call mpi_startall(nr,requests,ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpistartall')

      end subroutine comm_mpistartall
      
!=======================================================================
!
!     comm_mpibarrier
!
!=======================================================================

      subroutine comm_mpibarrier

      call mpi_barrier(mpi_comm_world,ierr)

      end subroutine comm_mpibarrier

!=======================================================================
!
!     comm_mpiwaitall
!
!=======================================================================

      subroutine comm_mpiwaitall(requests)

      integer, dimension(nr), intent(in) :: requests

!-----------------------------------------------------------------------

      integer, dimension(MPI_STATUS_SIZE,nr) :: status

      call mpi_waitall(nr,requests,status,ierr)
      if (ierr.ne.mpi_success) call comm_mpistop('mpiwaitall')

      end subroutine comm_mpiwaitall
      
!=======================================================================
!
!     comm_bufferupdate
!
!=======================================================================

      subroutine comm_bufferupdate(s,v,nx,ny,nz)

      type(comm_wld), intent(in) :: nx,ny,nz
      real(4), dimension(ns,nx%l,ny%l,nz%l) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l) :: v

!-----------------------------------------------------------------------
!
!     this routine updates all buffers. It is intended for
!     initialization and debugging. Buffer updates during regular
!     execution are handled in the subroutine sweep.
!
!-----------------------------------------------------------------------

      integer, parameter :: b1=yzbuf
      integer :: srs,srv
      integer, dimension(4) :: sshape
      integer, dimension(5) :: vshape

!.....copy to send buffer...............................................
      srs = ns*nx%l*ny%l*b1
      srv = 3*nv*nx%l*ny%l*b1
      send0(1:srs)         = reshape(s(:,:,:,nz%m:nz%m+b1-1),(/ srs /))
      send0(srs+1:srs+srv)=reshape(v(:,:,:,:,nz%m:nz%m+b1-1),(/ srv /))
      send1(1:srs)         = reshape(s(:,:,:,nz%n-b1+1:nz%n),(/ srs /))
      send1(srs+1:srs+srv)=reshape(v(:,:,:,:,nz%n-b1+1:nz%n),(/ srv /))

!.....do communication..................................................
      call comm_mpistartall(nz%requests)
      call comm_mpiwaitall(nz%requests)

!.....fetch from receive buffer.........................................
      sshape = (/ ns,nx%l,ny%l,b1 /)
      vshape = (/ 3,nv,nx%l,ny%l,b1 /)
      s(:,:,:,nz%m-b1:nz%m-1) = reshape(recv0(1:srs),sshape)
      v(:,:,:,:,nz%m-b1:nz%m-1)=reshape(recv0(srs+1:srs+srv),vshape)
      s(:,:,:,nz%n+1:nz%n+b1) = reshape(recv1(1:srs),sshape)
      v(:,:,:,:,nz%n+1:nz%n+b1)=reshape(recv1(srs+1:srs+srv),vshape)

      end subroutine comm_bufferupdate

!=======================================================================
!
!     comm_mpistop
!
!=======================================================================

      subroutine comm_mpistop(msg)

      implicit none     

      character(*) :: msg

      write(*,*) 'Process ',mr,' stopped with message ',msg, &
     &  ' by call to comm_mpistop'
      
      call mpi_abort(mpi_comm_world,ierr)
      stop

      end subroutine comm_mpistop

!=======================================================================

      end module comm_module

!***********************************************************************
