!=======================================================================
!
!     eos
!
!=======================================================================

      subroutine eos(d,e,p,cs,status)

      implicit none

      real, intent(in) :: d,e
      real, intent(out) :: p,cs
      integer, intent(out) :: status

!-----------------------------------------------------------------------
!
!     Input:
!     e  ... internal energy density [(c/dx)^2*c^2/G]
!     d  ... density [(c/dx)^2/G]
!
!     Output:
!     p  ... pressure [(c/dx)^2*c^2/G]
!     cs ... sound speed [c^2]
!     status ... 0 for success
!
!-----------------------------------------------------------------------

      real, parameter :: gamma=5./3.

      if ((d.gt.0.).and.(e.gt.0.)) then
        p = (gamma-1.)*e
        cs = sqrt(gamma*p/d)
        status = 0
      else
        status = 1
      endif

      end subroutine eos

!=======================================================================
