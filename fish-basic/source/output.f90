!=======================================================================
!
!     output
!
!=======================================================================

      subroutine output(num,time,s,v,nx,ny,nz)

      implicit none

      integer, intent(in) :: num
      type(comm_wld), intent(in) :: nx,ny,nz
      real, intent(in) :: time
      real(4), dimension(ns,nx%l,ny%l,nz%l), intent(in) :: s
      real(4), dimension(3,nv,nx%l,ny%l,nz%l), intent(in) :: v

!-----------------------------------------------------------------------
!
!     We write out binary data on planes in each direction.
!
!     dx ... length of one zone side
!     c  ... speed of light
!     G  ... gravitational constant
!
!-----------------------------------------------------------------------

      logical :: found
      integer, parameter :: np=3 !number of segments to write
      integer, parameter :: nred=5 !data reduction factor
      character(8) :: date
      character(11) :: suffix
      integer :: i0,j0,k0,ired,jred,kred
      integer :: status,i,j,k,m,is,imin,jmin,kmin,imax,jmax,kmax
      integer, dimension(np) :: xmin,ymin,zmin,xmax,ymax,zmax
      real(4) :: time4,dx4
      real(4), dimension(:,:,:,:), allocatable :: dov

      call prof_enter(14,'output')

!-----3D data-----------------------------------------------------------

!.....define range of segments in terms of global indices...............
      xmin(1) = 1               !xy-plane
      xmax(1) = nx%g
      ymin(1) = 1
      ymax(1) = ny%g
      zmin(1) = nz%g/2+1
      zmax(1) = nz%g/2+1

      xmin(2) = 1               !xz-plane
      xmax(2) = nx%g
      ymin(2) = ny%g/2+1
      ymax(2) = ny%g/2+1
      zmin(2) = 1
      zmax(2) = nz%g

      xmin(3) = nx%g/2+1        !yz-plane
      xmax(3) = nx%g/2+1
      ymin(3) = 1
      ymax(3) = ny%g
      zmin(3) = 1
      zmax(3) = nz%g

!.....open output file and write general information....................
      write(suffix,11) num,mr
11    format('_',i5,'_',i4)
      do m=1,11
        if (suffix(m:m).eq.' ') suffix(m:m) = '0'
      enddo
      open(1,file=trim(binpath)//'proc'//suffix(8:11)//'/'// &
     &  trim(name)//suffix//'.dat', &
     &  status='unknown',form='unformatted')
      call date_and_time(date)
      write(1) date
      time4 = time
      dx4 = dx
      write(1) time4,dx4
      write(1) np,num,mr

!.....loop over segments and determine overlap with local cube..........
      do is=1,np
        imin = max(xmin(is),nx%r+nx%m)  !overlap in local indices
        imax = min(xmax(is),nx%r+nx%n)
        jmin = max(ymin(is),ny%r+ny%m)
        jmax = min(ymax(is),ny%r+ny%n)
        kmin = max(zmin(is),nz%r+nz%m)
        kmax = min(zmax(is),nz%r+nz%n)
        write(1) 3*nv+ns-ndbg           !number of variables per zone
        write(1) imin,imax,jmin,jmax,kmin,kmax!overlap in global indices
        write(1) ((((/reshape(v(:,:,i,j,k),(/3*nv/)), &
     &    s(1:ns-ndbg,i,j,k)/), &
     &    i=imin-nx%r,imax-nx%r), &
     &    j=jmin-ny%r,jmax-ny%r), &
     &    k=kmin-nz%r,kmax-nz%r)
      enddo

!.....offset so that symmetry center belongs to domain overview.........
      i0 = mod(nx%g/2,nred)          !offset in standard global index
      j0 = mod(ny%g/2,nred)
      k0 = mod(nz%g/2,nred)
      imin = 1+(nx%r+nx%m-1-i0)/nred !reduced global index
      imax = 1+(nx%r+nx%n-1-i0)/nred
      jmin = 1+(ny%r+ny%m-1-j0)/nred
      jmax = 1+(ny%r+ny%n-1-j0)/nred
      kmin = 1+(nz%r+nz%m-1-k0)/nred
      kmax = 1+(nz%r+nz%n-1-k0)/nred
      allocate(dov(3*nv+ns-ndbg,1+imax-imin,1+jmax-jmin,1+kmax-kmin))
      do kred=1,1+kmax-kmin               !reduced local index
        k = 1+i0+(kmin+kred-2)*nred-nz%r    !standard local index
        do jred=1,1+jmax-jmin
          j = 1+j0+(jmin+jred-2)*nred-ny%r
          do ired=1,1+imax-imin
            i = 1+i0+(imin+ired-2)*nred-nx%r
            dov(:,ired,jred,kred) = (/reshape(v(:,:,i,j,k),(/3*nv/)), &
     &        s(1:ns-ndbg,i,j,k)/)
          enddo
        enddo
      enddo

!.....write domain overview.............................................
      write(1) nred,i0,j0,k0          !offset to fetch symmetry center
      write(1) 3*nv+ns-ndbg           !number of variables per zone
      write(1) imin,imax,jmin,jmax,kmin,kmax!overlap in global indices
      write(1) ((( dov(:,ired,jred,kred), &
     &  ired=1,1+imax-imin), &
     &  jred=1,1+jmax-jmin), &
     &  kred=1,1+kmax-kmin)
      deallocate(dov)

      close(1)
      call prof_exit(14)

!.....stop here if requested from outside...............................
      inquire(file=trim(txtpath)//'stop',exist=found)
      if (found) then
        write(6,*) 'Info: Stop file detected in ',trim(txtpath), &
     &    ', process ',mr
        call prof_write(mr)
        call mpi_finalize(ierr)
        stop
      endif

      end subroutine output

!=======================================================================
