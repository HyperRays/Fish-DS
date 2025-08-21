!=======================================================================
!
!     transpose forward
!
!=======================================================================

      subroutine transposef(s,v,nx,ny,nz)
      
        use mhd_parameter_module, only: ns, nv
        use prof_module

      implicit none

      integer, intent(in) :: nx,ny,nz
      real(4), dimension(ns,nx,ny,nz) :: s
      real(4), dimension(3,nv,nx,ny,nz) :: v

!-----------------------------------------------------------------------
!
!     This routine would be faster if tmp would have the size
!     dimension(3,nv,nx,ny,nz). However, this would add memory for
!     a copy of the largest array in the code. I hope the solution
!     below is a valuable compromise between memory and speed.
!
!-----------------------------------------------------------------------

      integer :: i,j,k,m,n,nt
      real(4), dimension(3,ny,nz,nx) :: tmp

      call prof_enter(6,'transpose')
!.....scalar............................................................
      do m=1,ns,3
        n = min(m+2,ns)
        nt = 1+n-m
        do k=1,nz
          do j=1,ny
            do i=1,nx
              tmp(1:nt,j,k,i) = s(m:n,i,j,k)
            enddo
          enddo
        enddo
        s(m:n,:,:,:) = reshape(tmp(1:nt,:,:,:),(/nt,nx,ny,nz/))
      enddo

!.....vector............................................................
      do m=1,nv
        do k=1,nz
          do j=1,ny
            do i=1,nx
              tmp(:,j,k,i) = v((/2,3,1/),m,i,j,k)
            enddo
          enddo
        enddo
        v(:,m,:,:,:) = reshape(tmp,(/3,nx,ny,nz/))
      enddo

      call prof_exit(6)

      end subroutine transposef

!=======================================================================
!
!     transpose backward
!
!=======================================================================

      subroutine transposeb(s,v,nx,ny,nz)

        use mhd_parameter_module, only: ns, nv
        use prof_module

      implicit none

      integer, intent(in) :: nx,ny,nz
      real(4), dimension(ns,nx,ny,nz) :: s
      real(4), dimension(3,nv,nx,ny,nz) :: v

!-----------------------------------------------------------------------

      integer :: i,j,k,m,n,nt
      real(4), dimension(3,nz,nx,ny) :: tmp

      call prof_enter(6,'transpose')

!.....scalar............................................................
      do m=1,ns,3
        n = min(m+2,ns)
        nt = 1+n-m
        do k=1,nz
          do j=1,ny
            do i=1,nx
              tmp(1:nt,k,i,j) = s(m:n,i,j,k)
            enddo
          enddo
        enddo
        s(m:n,:,:,:) = reshape(tmp(1:nt,:,:,:),(/nt,nx,ny,nz/))
      enddo

!.....vector............................................................
      do m=1,nv
        do k=1,nz
          do j=1,ny
            do i=1,nx
              tmp(:,k,i,j) = v((/3,1,2/),m,i,j,k)
            enddo
          enddo
        enddo
        v(:,m,:,:,:) = reshape(tmp,(/3,nx,ny,nz/))
      enddo

      call prof_exit(6)

      end subroutine transposeb

!=======================================================================
