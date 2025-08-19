!=======================================================================
!
!     sweep
!
!=======================================================================

      subroutine sweep(forward,update,yzbin,xb,yzbout,s,v,nx,ny,nz,dt)

      use comm_module

      implicit none

      external update

      type(comm_wld), intent(in) :: nx,ny,nz
      logical, intent(in) :: forward
      integer, intent(in) :: yzbin,xb,yzbout
      real, intent(in) :: dt
      real(4), dimension(ns,nx%l,ny%l,nz%l) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l) :: v

!-----------------------------------------------------------------------
!
!     Input/Output:
!     forward  ... direction of sweep (is passed to the update function)
!     update   ... function that updates s,v
!     yzbin    ... number of buffer zones with valid data on input
!     xb       ... number of border zones that get lost in update
!     yzbout   ... desired number of valid buffer zones on output
!     s,v      ... data arrays
!     nx,ny,nz ... space dimensions of data arrays
!     dt       ... time step
!
!     Others:
!     nxm = nx%m-yzbin   ... first valid zone on input
!     yzb = yzbout-yzbin ... length to be added to yzbuffer
!
!     actions at left boundary (right boundary is symmetric):
!
!                                 ___ = array section with original data
!                        nxm      --- = array section with updated data
!                         |       ... = array section with corrupt data
!                         V
!                        |yzbin|yzbin|
!                        |    xb    |   xb     |
!     |    xb    |  yzb  |    xb    |  yzb  |
!
!             |..........|____________________________ s _______________
!                                    |____ send buf ____|
!     initiate communications
!     |..................|___ sl ______________|
!                        |..........|----------------- s ---------------
!     continue when communications complete
!     |__ receive buf ___|
!     |______________________ sl ______________|
!     |..........|----------- sl ---|..........|
!             |..|------------------------------------ s ---------------
!
!                |    yzbout   |
!             |      yzbuf     |
!                               ^
!                               |
!                             nx%m
!
!-----------------------------------------------------------------------

      type(comm_wld) :: nxl,nxr
      integer :: yzb,b1,b2,nxm,nxtm,nxn,nxtn,srs,srv
      integer, dimension(4) :: sshape
      integer, dimension(5) :: vshape
      real(4), dimension(:,:,:,:), allocatable :: sl,sr
      real(4), dimension(:,:,:,:,:), allocatable :: vl,vr

      call prof_enter(11,'sweep')

!.....figure out dimensions.............................................
      yzb = yzbout-yzbin !length to be added to yzbuffer
      b1 = xb+yzb        !length of send and receive data
      b2 = 2*xb          !overlap of work arrays sl,sr,vl,vr with s,v
      nxm = nx%m-yzbin   !first valid zone on input
      nxtm = xb+yzb+1    !corresponding element in left work array
      nxn = nx%n+yzbin   !last valid zone on input
      nxtn = 2*xb        !corresponding element in right work array
      nxl%d = nx%d        !dimension of left work arrays sl,vl
      nxl%g = nx%g
      nxl%r = nx%r+nxm-nxtm
      nxl%m = nxtm-b1     !first boundary zone
      nxl%n = nxtm-1      !last boundary zone
      nxl%l = 3*xb+yzb
      nxr%d = nx%d        !dimension of right work arrays sr,vr
      nxr%g = nx%g
      nxr%r = nx%r+nxn-nxtn
      nxr%m = nxtn+1      !first boundary zone
      nxr%n = nxtn+b1     !last boundary zone
      nxr%l = 3*xb+yzb
      if (yzbout.gt.yzbuf) then
        write(6,*) 'Error: yzbout>yzbuf in sweep.f90!'
        call comm_mpistop('sweep')
      endif
      if (b1.gt.buf) then
        write(6,*) 'Error: xb+yzbout-yzbin>buf in sweep.f90!'
        call comm_mpistop('sweep')
      endif
      if (nx%n-nx%m.lt.max(b1,b2)) then
        write(6,*) 'Error: xb+max(yzbout-yzbin,xb)>nx%n-nx%m'
        write(6,*) 'in sweep.f90!'
        call comm_mpistop('sweep')
      endif
      if (xb.eq.0) then
        write(6,*) 'Error: xb=0 does not work with boundary()!'
        call comm_mpistop('sweep')
      endif
      allocate(sl(ns,nxl%l,ny%l,nz%l),sr(ns,nxr%l,ny%l,nz%l), &
     &  vl(3,nv,nxl%l,ny%l,nz%l),vr(3,nv,nxr%l,ny%l,nz%l))

!.....write boundaries to send buffer (single CPU version)...............
      call prof_enter(2,'copy and periodic boundaries')
      srs = ns*b1*ny%l*nz%l
      srv = 3*nv*b1*ny%l*nz%l
      ! Apply periodic boundaries directly for single CPU
      send0(1:srs) = &
     &  reshape(s(:,nx%n-yzbin-b1+1:nx%n-yzbin,:,:),(/ srs /))
      send0(srs+1:srs+srv) = &
     &  reshape(v(:,:,nx%n-yzbin-b1+1:nx%n-yzbin,:,:),(/ srv /))
      send1(1:srs) = &
     &  reshape(s(:,nx%m+yzbin:nx%m+yzbin+b1-1,:,:),(/ srs /))
      send1(srs+1:srs+srv) = &
     &  reshape(v(:,:,nx%m+yzbin:nx%m+yzbin+b1-1,:,:),(/ srv /))
      ! Copy to receive buffers immediately (no communication needed)
      recv0(1:srs+srv) = send0(1:srs+srv)
      recv1(1:srs+srv) = send1(1:srs+srv)

!.....no communication needed for single CPU............................
      ! call comm_mpistartall(nx%requests) - not needed

!.....copy boundaries of s,v into work arrays sl,vl,sr,vr...............
      sl(:,nxtm:nxtm+b2-1,:,:) = s(:,nxm:nxm+b2-1,:,:)       !left
      vl(:,:,nxtm:nxtm+b2-1,:,:) = v(:,:,nxm:nxm+b2-1,:,:)
      sr(:,nxtn-b2+1:nxtn,:,:) = s(:,nxn-b2+1:nxn,:,:)       !right
      vr(:,:,nxtn-b2+1:nxtn,:,:) = v(:,:,nxn-b2+1:nxn,:,:)
      call prof_exit(2)

!.....work on central parts in x-direction..............................
      call prof_enter(8,'work on center')
      call update(forward,s,v,nx,ny,nz,dt)
      call prof_exit(8)

!.....no wait needed for single CPU version.............................
      call prof_enter(2,'copy and periodic boundaries')
      ! call comm_mpiwaitall(nx%requests) - not needed for single CPU

!.....fetch boundaries of sl,vl,sr,vr from receive buffer...............
      sshape = (/ ns,b1,ny%l,nz%l /)
      vshape = (/ 3,nv,b1,ny%l,nz%l /)
      sl(:,nxtm-b1:nxtm-1,:,:) = reshape(recv0(1:srs),sshape)
      vl(:,:,nxtm-b1:nxtm-1,:,:)=reshape(recv0(srs+1:srs+srv),vshape)
      sr(:,nxtn+1:nxtn+b1,:,:) = reshape(recv1(1:srs),sshape)
      vr(:,:,nxtn+1:nxtn+b1,:,:)=reshape(recv1(srs+1:srs+srv),vshape)
      call prof_exit(2)

!.....work on boundaries in x-direction.................................
      call prof_enter(9,'work on boundary')
      call update(forward,sl,vl,nxl,ny,nz,dt)
      call update(forward,sr,vr,nxr,ny,nz,dt)

!.....copy result on boundaries of s and v..............................
      s(:,nxm-yzb:nxm+xb-1,:,:) = sl(:,nxtm-yzb:nxtm+xb-1,:,:)
      v(:,:,nxm-yzb:nxm+xb-1,:,:) = vl(:,:,nxtm-yzb:nxtm+xb-1,:,:)
      s(:,nxn-xb+1:nxn+yzb,:,:) = sr(:,nxtn-xb+1:nxtn+yzb,:,:)
      v(:,:,nxn-xb+1:nxn+yzb,:,:) = vr(:,:,nxtn-xb+1:nxtn+yzb,:,:)
      deallocate(sl,sr,vl,vr)
      call prof_exit(9)
      call prof_exit(11)

      end subroutine sweep

!=======================================================================
