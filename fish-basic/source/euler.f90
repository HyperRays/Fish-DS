!=======================================================================
!
!     euler
!
!=======================================================================

      subroutine euler(forward,s,v,nx,ny,nz,dt)

      use mhd_parameter_module
      use prof_module
      use comm_module

      implicit none

      type(comm_wld), intent(in) :: nx,ny,nz
      logical, intent(in) :: forward
      real, intent(in) :: dt
      real(4), dimension(ns,nx%l,ny%l,nz%l) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l) :: v

!-----------------------------------------------------------------------
!
!     This update produces invalid results in 8 boundary zones at each
!     end of the x-axis
!
!     1 zone lost by gradient of gravitational potential (not used)
!     1 zone lost for first order hydrodynamics update
!     2 zones lost for second order hydrodynamics update
!     1 zone lost for velocity smoothing
!     1 zone lost for first order magnetic field update
!     2 zones lost for second order magnetic field update
!
!-----------------------------------------------------------------------
 
      integer i,j,k,ip,jm,km
      real :: tmp1,tmp2,tmp3
      real(4), dimension(3,nx%l) :: b
      real, dimension(nx%l) :: fluxbx,b1x,vx
      real(4), dimension(nx%l,ny%l,nz%l) :: bsqr

!.....fluid update......................................................
      call prof_enter(21,'euler fluid')
      call prof_enter(45,'OMP_euler forward')
!$OMP PARALLEL DO SHARED(nx,ny,nz,s,v,dt,bsqr) PRIVATE(i,j,k,b)
      do k=1,nz%l-1
        do j=1,ny%l-1

!.....put magnetic field on zone centers................................
          do i=1,nx%l-1
            b(:,i) = v(:,mb,i,j,k)
            b(1,i) = 0.5*(b(1,i) + v(1,mb,i+1,j,k))
            b(2,i) = 0.5*(b(2,i) + v(2,mb,i,j+1,k))
            b(3,i) = 0.5*(b(3,i) + v(3,mb,i,j,k+1))
          enddo

!.....calculate magnetic energy.........................................
          do i=1,nx%l-1
            bsqr(i,j,k) = b(1,i)*b(1,i)+b(2,i)*b(2,i)+b(3,i)*b(3,i)
          enddo

!.....tvd...............................................................
          if (forward) then
            call tvdeul(s(1,1,j,k),v(1,1,1,j,k),nx%l,dt,b,bsqr(1,j,k))
          endif
        end do
      end do
!$OMP END PARALLEL DO
      call prof_exit(45)
      call prof_exit(21)

!.....bx,by update......................................................
      call prof_enter(22,'euler b')
      do k=1,nz%l
        do j=2,ny%l
          jm=j-1
          tmp2 = 0.5*(v(1,mv,1,jm,k)/s(md,1,jm,k) &
     &              + v(1,mv,1,j ,k)/s(md,1,j ,k))
          tmp3 = 0.5*(v(1,mv,2,jm,k)/s(md,2,jm,k) &
     &              + v(1,mv,2,j ,k)/s(md,2,j ,k))
          do i=3,nx%l
            tmp1 = tmp2
            tmp2 = tmp3
            tmp3 = 0.5*(v(1,mv,i,jm,k)/s(md,i,jm,k) &
     &                + v(1,mv,i,j ,k)/s(md,i,j ,k))
            vx(i-1) = 0.25*(tmp1+2.*tmp2+tmp3)
          enddo
          b1x=v(2,mb,:,j,k)

!.....tvd...............................................................
          call tvdb(fluxbx,b1x,vx,nx%l,dt)
          do i=5,nx%l-4
            ip = i+1
            v(2,mb,i ,j ,k) = b1x(i)
            v(1,mb,ip,j ,k) = v(1,mb,ip,j ,k) - fluxbx(i)
            v(1,mb,ip,jm,k) = v(1,mb,ip,jm,k) + fluxbx(i)
          enddo
        enddo
      enddo

!.....bx,bz update......................................................
      do k=2,nz%l
        km=k-1
        do j=1,ny%l
          tmp2 = 0.5*(v(1,mv,1,j,km)/s(md,1,j,km) &
     &              + v(1,mv,1,j,k )/s(md,1,j,k ))
          tmp3 = 0.5*(v(1,mv,2,j,km)/s(md,2,j,km) &
     &              + v(1,mv,2,j,k )/s(md,2,j,k ))
          do i=3,nx%l
            tmp1 = tmp2
            tmp2 = tmp3
            tmp3 = 0.5*(v(1,mv,i,j,km)/s(md,i,j,km) &
     &                + v(1,mv,i,j,k )/s(md,i,j,k ))
            vx(i-1) = 0.25*(tmp1+2.*tmp2+tmp3)
          enddo
          b1x=v(3,mb,:,j,k)

!.....tvd...............................................................
          call tvdb(fluxbx,b1x,vx,nx%l,dt)
          do i=5,nx%l-4
            ip = i+1
            v(3,mb,i,j,k ) = b1x(i)
            v(1,mb,ip,j,k ) = v(1,mb,ip,j,k ) - fluxbx(i)
            v(1,mb,ip,j,km) = v(1,mb,ip,j,km) + fluxbx(i)
          enddo
        enddo
      enddo
      call prof_exit(22)

!.....fluid update......................................................
      call prof_enter(21,'euler fluid')
      call prof_enter(48,'OMP_euler backward')
!$OMP PARALLEL DO SHARED(nx,ny,nz,s,v,dt,bsqr) PRIVATE(i,j,k,b)
      do k=1,nz%l-1
        do j=1,ny%l-1

!.....put magnetic field on zone centers................................
          do i=1,nx%l-1
            b(:,i) = v(:,mb,i,j,k)
            b(1,i) = 0.5*(b(1,i) + v(1,mb,i+1,j,k))
            b(2,i) = 0.5*(b(2,i) + v(2,mb,i,j+1,k))
            b(3,i) = 0.5*(b(3,i) + v(3,mb,i,j,k+1))
          enddo

!.....tvd...............................................................
          if (.not.forward) then
            call tvdeul(s(1,1,j,k),v(1,1,1,j,k),nx%l,dt,b,bsqr(1,j,k))
          endif

!.....calculate magnetic energy.........................................
          do i=1,nx%l-1
            bsqr(i,j,k) = b(1,i)*b(1,i)+b(2,i)*b(2,i)+b(3,i)*b(3,i)
          enddo
        end do
      end do
!$OMP END PARALLEL DO
      call prof_exit(48)

!.....subtract magnetic energy from s(ms,:,:,:).........................
      s(ms,:,:,:) = s(ms,:,:,:) - 0.5*bsqr

      call prof_exit(21)
 
      end subroutine euler

!=======================================================================
