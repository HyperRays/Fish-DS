!=======================================================================
!
!     calculate the CFL condition
!
!=======================================================================

      subroutine calcfl(s,v,nx,ny,nz,cfl)

      use comm_module
      use prof_module

      implicit none

      type(comm_wld), intent(in) :: nx,ny,nz
      real(4) :: s(ns,nx%l,ny%l,nz%l)
      real(4), intent(in) :: v(3,nv,nx%l,ny%l,nz%l)
      real, intent(out) :: cfl

!-----------------------------------------------------------------------

      integer :: i,j,k,status
      real :: csig,bsqr,p,cs,rho,c,cm,cmax
      real, dimension(ns) :: ss
      real, dimension(3) :: vv,ww,bb
 
!.....calculate time step limit.........................................
      cm = 0.
      do k=1,nz%l-1
        do j=1,ny%l-1
          do i=1,nx%l-1

!.....limit on Euler step...............................................
            ss = s(:,i,j,k)
            call eos(ss(md),ss(ms),p,cs,status)
            if ((ss(ms).le.0.).or.(status.ne.0)) then
              if ((nz%m.le.k).and.(k.le.nz%n) &
     &          .and.(ny%m.le.j).and.(j.le.ny%n) &
     &          .and.(nx%m.le.i).and.(i.le.nx%n)) then
                write(6,*) 'Error: Problem with eos at zone'
                write(6,*) nx%r+i,ny%r+j,nz%r+k
                write(6,*) 'd,e: ',ss(md),ss(ms)/ss(md)
                call comm_mpistop('calcfl')
              endif
            endif
            bb(1) = 0.5*(v(1,mb,i,j,k)+v(1,mb,i+1,j,k))
            bb(2) = 0.5*(v(2,mb,i,j,k)+v(2,mb,i,j+1,k))
            bb(3) = 0.5*(v(3,mb,i,j,k)+v(3,mb,i,j,k+1))
            bsqr  = bb(1)*bb(1) + bb(2)*bb(2) + bb(3)*bb(3)
            csig = sqrt(cs*cs + 2.*bsqr/ss(md)) !empirical factor 2
            vv = v(:,mv,i,j,k)/ss(md)
            c = sqrt(sum(vv*vv)) + csig

!.....cfl time step.....................................................
            cm = max(cm,c)
          end do
        end do
      end do
      call prof_enter(31,'single CPU max in calcfl')
      cmax = cm  ! Single CPU version - no MPI reduction needed
      call prof_exit(31)
      write(6,*) 'cmax:',cmax
      cfl = 1./cmax

      end subroutine calcfl

!=======================================================================
