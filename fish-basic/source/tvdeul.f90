!=======================================================================
!
!     tvdeul
!
!=======================================================================

      subroutine tvdeul(s,v,n,dt,b,bsqrin)

      use mhd_parameter_module

      implicit none

      integer, intent(in) :: n
      real, intent(in) :: dt
      real(4), dimension(ns,n) :: s
      real(4), dimension(3,nv,n) :: v
      real(4), dimension(3,n), intent(in) :: b
      real(4), dimension(n), intent(in) :: bsqrin

!-----------------------------------------------------------------------
!
!     Input:
!     s(:,1:n)    ... array of scalars along current x-direction
!     v(3,:,1:n)  ... array of vectors along current x-direction
!     n           ... number of points in x-direction
!     dt          ... time step
!     b(3,1:n)    ... magnetic field vector
!     bsqrin(1:n) ... 2*(magnetic energy density)
!     
!     s(md,:)   ... density
!     s(ms,:)   ... internal energy density
!     v(1,mv,:) ... x-component of momentum density
!     v(2,mv,:) ... y-component of momentum density
!     v(3,mv,:) ... z-component of momentum density
!
!     Output:
!     s(:,4:n-3)   ... updated array of scalars
!     v(3,:,4:n-3) ... updated array of vectors
!     
!     s(md,:)   ... density
! !!  s(ms,:)   ... internal energy density + magnetic energy density !!
!     v(1,mv,:) ... x-component of momentum density
!     v(2,mv,:) ... y-component of momentum density
!     v(3,mv,:) ... z-component of momentum density
!
!     The output is by 3 points short at each boundary because the
!     first order advection half timestep consumes 1 point,
!     and the second order full timestep consumes 2 points.
!
!     In order to keep memory access efficient, the routine is organised
!     as shown in the following example calculating
!
!     s_i = s_i + 0.5*(grad_i+1/2 + grad_i-1/2),
!
!     where grad_i+1/2 = u_i+1 - u_i:
!
!     do i=1,n
!       u2 = u3
!       u3 = s(mx,i)                         !fetch data point from array
!       grad1p = grad2p
!       grad2p = u3 - u2
!       s(mx,i-1) = u2 + 0.5*(grad2p+grad1p) !put result into array
!     enddo
!
!-----------------------------------------------------------------------

      integer, parameter :: Lvx=1, Lvy=2, Lvz=3
      integer, parameter :: Ld=4,Le=5
      integer, parameter :: nu=5
      integer :: i,j,status
      real :: e,p,ps,cs,cf,c4,c5
      real :: vsqr,bsqr5,bsqr6,bdotv
      real, dimension(3) :: b5,b6,vadv
      real, dimension(nu) :: u4,u5,uu3,uu4,uu5,uu6
      real, dimension(nu) :: ff5,ff6,f3,f4,f5,cfdu3p,cfdu4p
      real, dimension(nu) :: fflx4p,fflx5p,flx2p,flx3p
      real, dimension(nu,2) :: ddn,dup,d3p,tmp !1: right-going waves
                                               !2: left-going waves

!.....initialisations....................................................
      b6 = 0.
      bsqr6 = 0.
      c5 = 1.
      cfdu4p = 0.
      dup = 0.
      f4 = 0.
      f5 = 0.
      ff6 = 0.
      fflx5p = 0.
      flx3p = 0.
      u5 = 0.
      uu4 = 0.
      uu5 = 0.
      uu6 = 0.

!.....loop along x-direction............................................
      do i=1,n

!.....preparation: fetch data from array................................
        uu3 = uu4
        uu4 = uu5
        uu5 = uu6
        uu6((/Lvx,Lvy,Lvz,Ld,Le/)) = (/v(:,mv,i),s((/md,ms/),i)/)
        b5 = b6
        b6 = b(:,i)

!.....half step: call equation of state.................................
        call eos(uu6(Ld),uu6(Le),p,cs,status)
        vadv = uu6(Lvx:Lvz)/uu6(Ld)
        vsqr = uu6(Lvx)*vadv(1) + uu6(Lvy)*vadv(2) + uu6(Lvz)*vadv(3)
        bsqr5 = bsqr6
        bsqr6  = b6(1)*b6(1) + b6(2)*b6(2) + b6(3)*b6(3)
        uu6(Le) = uu6(Le) + 0.5*(vsqr + bsqrin(i)) !uu6(Le) is now total energy
        ps = p + 0.5*bsqr6

!.....half step: calculate fluxes at zone centers.......................
        ff5 = ff6
        ff6 = uu6*vadv(1)
        ff6(Lvx) = ff6(Lvx) + ps
        ff6(Lvx:Lvz) = ff6(Lvx:Lvz) - b6*b6(1)
        bdotv = b6(1)*vadv(1) + b6(2)*vadv(2) + b6(3)*vadv(3)
        ff6(Le) = ff6(Le) + ps*vadv(1) - b6(1)*bdotv

!.....half step: fluxes on zone boundaries..............................
        fflx4p = fflx5p
        fflx5p = 0.5*(ff5+ff6)

!.....half step: update zone center values..............................
        if (i.le.2) cycle
        u4 = u5
        u5 = uu5 - (fflx5p-fflx4p)*0.5*dt

!.....full step: call equation of state.................................
        vadv = u5(Lvx:Lvz)/u5(Ld)
        vsqr = u5(Lvx)*vadv(1) + u5(Lvy)*vadv(2) + u5(Lvz)*vadv(3)
        e = u5(Le) - 0.5*(vsqr + bsqr5)
        call eos(u5(Ld),e,p,cs,status)
        ps = p + 0.5*bsqr5
        c4 = c5
        c5 = abs(vadv(1)) + sqrt(cs*cs + bsqr5/u5(Ld))
        cf = max(c4,c5)                    !freezing speed
        cfdu3p = cfdu4p
        cfdu4p = cf*(u5-u4)

!.....full step: calculate fluxes at zone centers.......................
        f3 = f4
        f4 = f5
        f5 = u5*vadv(1)
        f5(Lvx) = f5(Lvx) + ps
        f5(Lvx:Lvz) = f5(Lvx:Lvz) - b5*b5(1)
        bdotv = b5(1)*vadv(1) + b5(2)*vadv(2) + b5(3)*vadv(3)
        f5(Le) = f5(Le) + ps*vadv(1) - b5(1)*bdotv

!.....full step: calculate flux differences for flux limiter............
        ddn = dup
        dup(:,1) = 0.5*(f4-f3+cfdu3p)
        dup(:,2) = 0.5*(f5-f4-cfdu4p)

!.....full step: Van Leer flux limiter..................................
        tmp = ddn*dup
        where (tmp.gt.0.)
          d3p = 2.*tmp/(ddn+dup)
        else where
          d3p = 0.
        end where

!.....full step: update zone center values according to limited fluxes..
        flx2p = flx3p
        flx3p = 0.5*(f3+f4-cfdu3p + d3p(:,1)-d3p(:,2))
        j = i-3
        if (j.le.3) cycle
        uu3 = uu3 - (flx3p-flx2p)*dt

!.....finalising: subtracting kinetic energy from total energy..........
        vadv = uu3(Lvx:Lvz)/uu3(Ld)
        vsqr = uu3(Lvx)*vadv(1) + uu3(Lvy)*vadv(2) + uu3(Lvz)*vadv(3)
        uu3(Le) = uu3(Le) - 0.5*vsqr !uu3(Le) is now E(internal+magnetic)

!.....finalising: writing result to output array........................
        v(:,mv,j) = uu3(Lvx:Lvz)
        s((/md,ms/),j) = uu3((/Ld,Le/))
      enddo

      end subroutine tvdeul

!=======================================================================
