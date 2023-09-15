!=======================================================================
!
!     tvdb
!
!=======================================================================

      subroutine tvdb(flux,b,vg,n,dt)

      implicit none

      integer, intent(in) :: n
      real, intent(in) :: dt
      real, dimension(n) :: flux,b,vg

!-----------------------------------------------------------------------
!
!     unlike the B field, the flux lives on the right cell boundary
!
!     This routine requires input for zones 2 ... n-1 and provides
!     results for zones 5 ... n-4.
!
!-----------------------------------------------------------------------

      integer :: i,i4,i1
      real :: b4,b3,b2,b1,bb3,tmp
      real :: f3,f2,f1,f0,w3,w2,w1,dw2,dw1,d2,d1
      real :: vg4,vg3,vh3,vh2,vh1

!.....initialize........................................................
      f1 = 0.
      w2 = 0.
      dw2 = 0.
      d2 = 0.
      vh2 = 0.
      b2 = 0.
      vh3 = 0.
      f3 = 0.
      w3 = 0.
      b3 = 0.
      vg4 = 0.
      b4 = 0.

!.....down shift variables belonging to second update...................
      do i=2,n-1
        f0 = f1
        vh1 = vh2
        w1 = w2
        dw1 = dw2
        d1 = d2
        b1 = b2
        vh2 = vh3
        w2 = w3
        b2 = b3

!.....down shift variables belonging to first update....................
        f2 = f3
        vg3 = vg4
        b3 = b4
        vg4 = vg(i)
        b4 = b(i)

!.....first update......................................................
        vh3 = 0.5*(vg3+vg4)
        if (vh3.gt.0.) then
          f3 = b3*vg3
        else
          f3 = b4*vg4
        endif
        bb3 = b3 - 0.5*(f3-f2)*dt

!.....second update.....................................................
        w3 = vg3*bb3
        dw2 = 0.5*(w3-w2)
        tmp = dw1*dw2
        if (tmp.gt.0.) then
          d2 = 2.*tmp/(dw1+dw2)
        else
          d2 = 0.
        endif
        if (vh1.gt.0) then
          f1=(w1+d1)*dt
        else
          f1=(w2-d2)*dt
        end if
        if (i.gt.7) then
          flux(i-3) = f1
          b(i-3) = b1 - (f1-f0)
        endif
      end do

      end subroutine tvdb

!=======================================================================
