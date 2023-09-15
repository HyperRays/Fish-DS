!***********************************************************************
!
!     mhd_module, provides module-access to FISH
!
!***********************************************************************
!
!     The indexing conventions are defined in comm_module.f
!     The solution vectors s,v are described in init.f90
!     The B field is stored on the left side of each cell
!
!     Version V
!
!     A similar version of this code is documented in
!     "FISH: A 3D parallel MHD code for astrophysical applications"
!     R. Kaeppeli, S. C. Whitehouse, S. Scheidegger,
!     U.-L. Pen, M. Liebendoerfer, http://arxiv.org/abs/0910.2854
!
!-----------------------------------------------------------------------

      module mhd_module

      use comm_module
      use mhd_parameter_module
      use prof_module

      implicit none

!.....interface for empty sweeps for boundary update....................
      interface empty
        subroutine empty(forward,s,v,nx,ny,nz,dummy)
        use mhd_parameter_module
        use comm_module
        type(comm_wld), intent(in) :: nx,ny,nz
        logical, intent(in) :: forward
        real, intent(in) :: dummy
        real(4), dimension(ns,nx%l,ny%l,nz%l) :: s
        real(4), dimension(3,nv,nx%l,ny%l,nz%l) :: v
        end subroutine empty
      end interface empty

!.....interface for Euler step..........................................
      interface euler
        subroutine euler(forward,s,v,nx,ny,nz,dt)
        use mhd_parameter_module
        use comm_module
        type(comm_wld), intent(in) :: nx,ny,nz
        logical, intent(in) :: forward
        real, intent(in) :: dt
        real(4), dimension(ns,nx%l,ny%l,nz%l) :: s
        real(4), dimension(3,nv,nx%l,ny%l,nz%l) :: v
        end subroutine euler
      end interface euler

!.....global variables..................................................
      type(comm_wld) :: nx,ny,nz
      real(4), dimension (:,:,:,:), allocatable :: s
      real(4), dimension (:,:,:,:,:), allocatable :: v
 
      contains

!=======================================================================
!
!     mhd_memory
!
!=======================================================================

      subroutine mhd_memory

!.....initialize MPI and figure out array dimensions....................
      call comm_mpinit(nx,ny,nz)
      allocate(s(ns,nx%l,ny%l,nz%l))
      allocate(v(3,nv,nx%l,ny%l,nz%l))

      end subroutine mhd_memory

!=======================================================================
!
!     mhd_initialization
!
!=======================================================================

      subroutine mhd_initialization(iter,time)

      integer,intent(out) :: iter
      real, intent(out) :: time

!-----------------------------------------------------------------------

      logical :: forward=.true.
      real :: dummy=0.

!.....read file or setup initial conditions.............................
      call readdmp(iter,time,s,v,nx,ny,nz)
!.....update buffers with initial conditions............................
      call sweep(forward,empty,0,1,yzbuf,s,v,nx,ny,nz,dummy)        !x
      call transposef(s,v,nx%l,ny%l,nz%l)
      call sweep(forward,empty,0,1,yzbuf,s,v,ny,nz,nx,dummy)        !y
      call transposef(s,v,ny%l,nz%l,nx%l)
      call sweep(forward,empty,0,1,yzbuf,s,v,nz,nx,ny,dummy)        !z
      call transposef(s,v,nz%l,nx%l,ny%l)

      end subroutine mhd_initialization

!=======================================================================
!
!     mhd_write
!
!=======================================================================

      subroutine mhd_write(time,iter)

      integer, intent(in) :: iter
      real, intent(in) :: time

      if ((iter.ne.restart*skip).and.(mod(iter,skip).eq.0)) then
        if (mod(iter/skip,dmpskip).eq.0.) then
          call writdmp(iter/skip,time,s,v,nx,ny,nz)
        endif
        call output(iter/skip,time,s,v,nx,ny,nz)
      endif

      end subroutine mhd_write

!=======================================================================
!
!     mhd_step
!
!=======================================================================

      subroutine mhd_step(dt)

      real, intent(in) :: dt

!-----------------------------------------------------------------------

      logical, parameter :: forward=.true.,backward=.false.
      integer :: i,j,k

      call prof_enter(2,'mhd_step')

!.....sweep forward.....................................................
!     write(6,*) 'fx'
      call sweep(forward,euler,0,7,yzbuf,s,v,nx,ny,nz,dt)
      call transposef(s,v,nx%l,ny%l,nz%l)

!     write(6,*) 'fy'
      call sweep(forward,euler,0,7,yzbuf,s,v,ny,nz,nx,dt)
      call transposef(s,v,ny%l,nz%l,nx%l)

!     write(6,*) 'fz'
      call sweep(forward,euler,0,7,yzbuf,s,v,nz,nx,ny,dt)

!.....sweep backward....................................................
!     write(6,*) 'bz'
      call sweep(backward,euler,0,7,yzbuf,s,v,nz,nx,ny,dt)
      call transposeb(s,v,nz%l,nx%l,ny%l)

!     write(6,*) 'by'
      call sweep(backward,euler,0,7,yzbuf,s,v,ny,nz,nx,dt)
      call transposeb(s,v,ny%l,nz%l,nx%l)

!     write(6,*) 'bx'
      call sweep(backward,euler,0,7,yzbuf,s,v,nx,ny,nz,dt)
      call prof_exit(2)

      end subroutine mhd_step

!=======================================================================
!
!     mhd_cfl
!
!=======================================================================

      subroutine mhd_cfl(dt)

      real, intent(out) :: dt

!-----------------------------------------------------------------------

      real :: cfl

!.....determine time step based on sound speed..........................
      call prof_enter(3,'calcfl')
      call calcfl(s,v,nx,ny,nz,cfl)
      call prof_exit(3)
      if (mr.eq.0) then
        if (cfl.lt.dt) then
          write(6,*) 'previous step: ',dt
          write(6,*) 'new cfl cond.: ',cfl
          write(6,*) 'Warning: cfl exceeded'
        endif
      endif
      dt = cflsav*cfl

      end subroutine mhd_cfl

!=======================================================================
!
!     further internal subroutines
!
!=======================================================================

      include 'calcfl.f90'
      include 'transpose.f90'
      include 'sweep.f90'
      include 'readdmp.f90'
      include 'writdmp.f90'
      include 'output.f90'

!=======================================================================

      end module mhd_module

!***********************************************************************
