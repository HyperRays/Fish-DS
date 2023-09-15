!=======================================================================
!
!     init
!
!=======================================================================

      subroutine init(s,v,nx,ny,nz)

      use mhd_parameter_module
      use comm_module

      implicit none

      type(comm_wld), intent(in) :: nx,ny,nz
      real(4), dimension(ns,nx%l,ny%l,nz%l), intent(out) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l), intent(out) :: v

!-----------------------------------------------------------------------
!
!     scalar: s(md,:,:,:) ... density
!             s(ms,:,:,:) ... internal energy density
!             s(mg,:,:,:) ... gravitational potential
!
!             v(:,mv,:,:,:) ... momentum density
!             v(:,mb,:,:,:) ... magnetic field
!
!-----------------------------------------------------------------------

      integer :: im,i,j,k,status
      real, parameter :: epsilon=0.1
      real :: p0
  
      if (mhd_problem.eq.1) then
        p0=3./5.
        s=0
        v=0
        s(md,:,:,:)=1
        s(ms,:,:,:)=1.5*p0

        v(1,mv,:,1,1)=0.0001 &
     &    * sin( (/ (2*3.14159*(nx%r+i)/nx%g, i=1,nx%l) /) )
        s(md,:,1,1)=s(md,:,1,1)+v(1,mv,:,1,1)
        s(ms,:,1,1)=s(ms,:,1,1)+p0*v(1,mv,:,1,1)*5./3./(2./3.)
      endif

! circularly polarized alven wave:
! background : rho=B_x=1 all others zero
      if (mhd_problem.eq.2) then
        s=0
        v=0
        v(1,mb,:,:,:)=1
        s(md,:,:,:)=1
        s(ms,:,:,:)=0.001  ! to keep things stable
        do i=1,nx%l
          v(2,mv,i,:,:)=epsilon*sin( 2*3.14159*(nx%r+i)/nx%g )
          v(3,mv,i,:,:)=epsilon*cos( 2*3.14159*(nx%r+i)/nx%g )
        enddo
        v(2,mb,:,:,:)=-v(2,mv,:,:,:)
        v(3,mb,:,:,:)=-v(3,mv,:,:,:)

        if (ny%m.lt.2) write(6,*) 'Warning: ny%m<2 in init.f90'
        do i=2,ny%l
          im=i-1
          v(2,mb,:,i,:)=(v(2,mb,:,i,:)+v(2,mb,:,im,:))/2
        enddo
        if (nz%m.lt.2) write(6,*) 'Warning: nz%m<2 in init.f90'
        do i=2,nz%l
          im=i-1
          v(3,mb,:,:,i)=(v(3,mb,:,:,i)+v(3,mb,:,:,im))/2
        enddo
      endif

! alven wave:
! background : rho=B_x=1 all others zero
! \dot{v_y}+\div_x(-B_y)=0
! \dot{B_y}+\div_x (       -   v_y)     = 0
! let v_y=\epsilon sin(2 \pi (x-t)/L)
! then B_y=-v_y
      if (mhd_problem.eq.3) then
        s=0
        v=0
        v(1,mb,:,:,:)=1
        s(md,:,:,:)=1
        s(ms,:,:,:)=.001  ! to keep things stable
        do i=1,nx%l
          v(2,mv,i,:,:)=0.1*sin( 2*3.14159*(nx%r+i)/nx%g )
        enddo
        v(2,mb,:,:,:)=-v(2,mv,:,:,:)
        v(2,mb,:,:,:)=(v(2,mb,:,:,:)+cshift(v(2,mb,:,:,:),-1,1))/2
      endif

! magnetosonic
      if (mhd_problem.eq.4) then
        s=0
        v=0
        v(2,mb,:,:,:)=1
        s(md,:,:,:)=1
        s(ms,:,:,:)=.001  ! to keep things stable
        do i=1,nx%l
          v(1,mv,i,:,:)=0.01*sin( 2*3.14159*(nx%r+i)/nx%g )
        enddo
        v(2,mb,:,:,:)=v(2,mb,:,:,:)+v(1,mv,:,:,:)
        s(md,:,:,:)=v(2,mb,:,:,:)
      endif

      return
  
      end subroutine init

!=======================================================================
