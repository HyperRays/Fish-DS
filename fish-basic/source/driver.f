!=======================================================================
!
!     FISH driver
!
!     A similar version of this code is documented in
!     "FISH: A 3D parallel MHD code for astrophysical applications"
!     R. Kaeppeli, S. C. Whitehouse, S. Scheidegger,
!     U.-L. Pen, M. Liebendoerfer, http://arxiv.org/abs/0910.2854
!
!=======================================================================

      program fish

!.....driver............................................................
      use prof_module
      use mhd_module

      implicit none

!.....variables.........................................................
      character(len=4) :: fn
      integer :: iter
      real :: time,dt,vstep
      character(len=50) :: dir

      ! initialize dt to an value less than cflsav so that during 
      ! mhd_cfl subprocess it does not cmp to NaN value 
      ! Note: mhd_cfl anyways resets the dt value
      dt = 0

!.....initialize numerics...............................................
      call prof_initial
      call prof_directory(logdir)
      call mhd_memory
      call read_restart_loc

!.....redirect output...................................................
      if (mr.ne.0) then
        write(fn,22) mr
        do iter=1,4
          if (fn(iter:iter).eq.' ') fn(iter:iter) = '0'
        enddo
22      format(i4)
        write(dir, 77) trim(logdir), fn
77      format(A,'scratch_',A,'.log')
        open(6,file=trim(dir),status='unknown')
      endif

!-----build initial state-----------------------------------------------
      call mhd_initialization(iter,time)

!-----do time steps-----------------------------------------------------
      do

!.....determine time step...............................................
        if (time.ge.tf) goto 100
        iter=iter+1
        call mhd_cfl(dt)
        write (6,33) 'step/time/dt [dx/c]:',iter,time,dt
33      format(a24,i12,3g12.4)
        time = time + 2.*dt
        if (time.gt.tf) then
          dt = dt-time+tf
          time = tf
        endif

!.....MHD double time step..............................................
        call mhd_step(dt)

!.....writeout..........................................................
        call mhd_write(time,iter)
        call prof_write(mr)
        write(6,11)
11      format(72('-'))
      enddo

!.....finalize MPI......................................................
100   continue
      if (mr.eq.0) write(6,*) 'Info: done'
      call mpi_finalize(ierr)

      end program fish

!=======================================================================
