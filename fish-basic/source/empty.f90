!=======================================================================
!
!     empty
!
!=======================================================================

      subroutine empty(forward,s,v,nx,ny,nz,dummy)

      use mhd_parameter_module
      use comm_module

      implicit none

      type(comm_wld), intent(in) :: nx,ny,nz
      logical, intent(in) :: forward
      real, intent(in) :: dummy
      real(4), dimension(ns,nx%l,ny%l,nz%l) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l) :: v

      end subroutine empty

!=======================================================================
