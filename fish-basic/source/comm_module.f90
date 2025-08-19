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

!.....standard variables (non-MPI version)..............................
      integer :: ncpu,mr,ierr
      integer, parameter :: nr=4 !number of requests per direction (dummy)

!.....message tags (dummy)..............................................
      integer :: xtagdw=1,xtagup=2,ytagdw=3,ytagup=4,ztagdw=5,ztagup=6

!.....storage for communication buffers (dummy)........................
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
!     - initialize single CPU version (no MPI)
!     - set workload distribution parameters for type(comm_wld)
!     - setup for single process execution
!
!-----------------------------------------------------------------------

      integer :: npercpu,nxmn,nymn,nzmn,i,j,k,imax,jmax,kmax
      integer :: mr0,mr1

!.....initialization for single CPU.....................................
      ncpu = 1
      mr = 0
      ierr = 0
      write(6,*) 'Single CPU mhd code initialized'
      write(6,22) nxg,nyg,nzg,ncpu
22    format(' Cube of sides = ',i4,'*',i4,'*',i4,' over ', &
     &    i6,' processes')

!.....single CPU decomposition (no decomposition needed)................
      nxmn = nxg
      nymn = nyg
      nzmn = nzg
      
!.....3D index of mr (single process)...................................
      imax = 1
      jmax = 1
      kmax = 1
      i = 1
      j = 1
      k = 1

!.....assign comm_wld for x-direction...................................
      nx%d = 1
      nx%g = nxg
      nx%r = -yzbuf
      nx%m = yzbuf+1
      nx%n = yzbuf+nxmn
      nx%l = 2*yzbuf+nxmn
      write(*,*) 'X axis: full domain ',nx%n-nx%m+1,' cells'
      if (nx%n-nx%m+1.gt.nmax) then
        write(6,*) 'Error: parameter nmax=',nmax,'is too small!'
        stop
      endif

!.....assign comm_wld for y-direction...................................
      ny%d = 2
      ny%g = nyg
      ny%r = -yzbuf
      ny%m = yzbuf+1
      ny%n = yzbuf+nymn
      ny%l = 2*yzbuf+nymn
      write(*,*) 'Y axis: full domain ',ny%n-ny%m+1,' cells'
      if (ny%n-ny%m+1.gt.nmax) then
        write(6,*) 'Error: parameter nmax=',nmax,'is too small!'
        stop
      endif
 
!.....assign comm_wld for z-direction...................................
      nz%d = 3
      nz%g = nzg
      nz%r = -yzbuf
      nz%m = yzbuf+1
      nz%n = yzbuf+nzmn
      nz%l = 2*yzbuf+nzmn
      write(*,*) 'Z axis: full domain ',nz%n-nz%m+1,' cells'
      if (nz%n-nz%m+1.gt.nmax) then
        write(6,*) 'Error: parameter nmax=',nmax,'is too small!'
        stop
      endif

!.....no communication setup needed for single CPU.....................
      nx%requests = 0
      ny%requests = 0
      nz%requests = 0

      end subroutine comm_mpinit

!=======================================================================
!
!     comm_mpistartall
!
!=======================================================================

      subroutine comm_mpistartall(requests)

      integer, dimension(nr), intent(in) :: requests

!.....no-op for single CPU version......................................
      return

      end subroutine comm_mpistartall
      
!=======================================================================
!
!     comm_mpibarrier
!
!=======================================================================

      subroutine comm_mpibarrier

!.....no-op for single CPU version......................................
      return

      end subroutine comm_mpibarrier

!=======================================================================
!
!     comm_mpiwaitall
!
!=======================================================================

      subroutine comm_mpiwaitall(requests)

      integer, dimension(nr), intent(in) :: requests

!-----------------------------------------------------------------------

!.....no-op for single CPU version......................................
      return

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
!     this routine updates all buffers using periodic boundary 
!     conditions for single CPU version. No MPI communication needed.
!
!-----------------------------------------------------------------------

      integer, parameter :: b1=yzbuf

!.....apply periodic boundary conditions in z-direction................
!     Copy from top of domain to bottom buffer
      s(:,:,:,nz%m-b1:nz%m-1) = s(:,:,:,nz%n-b1+1:nz%n)
      v(:,:,:,:,nz%m-b1:nz%m-1) = v(:,:,:,:,nz%n-b1+1:nz%n)
!     Copy from bottom of domain to top buffer
      s(:,:,:,nz%n+1:nz%n+b1) = s(:,:,:,nz%m:nz%m+b1-1)
      v(:,:,:,:,nz%n+1:nz%n+b1) = v(:,:,:,:,nz%m:nz%m+b1-1)

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
      
      stop

      end subroutine comm_mpistop

!=======================================================================

      end module comm_module

!***********************************************************************
