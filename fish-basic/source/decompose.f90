!***********************************************************************
!
!     subroutine decompose calculates optimum domain decomposition
!     (translated from subroutine written by Markus Langer(?), Basel 2008)
!
!***********************************************************************

      subroutine decompose(ncpu,nx,ny,nz,nxmn,nymn,nzmn)

      implicit none

      integer, intent(in) :: ncpu            !number of CPU's
      integer, intent(in) :: nx,ny,nz        !size of full domain
      integer, intent(out) :: nxmn,nymn,nzmn !size of decomposed domain

!-----------------------------------------------------------------------
!
!     Input:
!     ncpu ... number of cpu's
!     nx   ... number of zones in x-direction
!     ny   ... number of zones in y-direction
!     nz   ... number of zones in z-direction
!
!     Output:
!     nxmn ... number of zones per domain in x-direction
!     nymn ... number of zones per domain in y-direction
!     nzmn ... number of zones per domain in z-direction
!
!-----------------------------------------------------------------------

      
      logical :: valid                  !property of solution
      integer :: lx,ly,lz               !subdomain lx*ly*lz
      integer, allocatable, dimension(:,:) :: solution
      integer, allocatable, dimension(:) :: value
      integer, allocatable, dimension(:) :: xfactors,yfactors,zfactors
      integer :: fx,fy,fz               !number of factors
      integer :: i,n,n1,n2,n3           !loop indices
      integer, dimension(1) :: pos      !position of minimum
      integer, parameter :: maxsol=100  !maximum number of solutions

!.....determine factors in each direction...............................
      allocate(xfactors(nx/2),yfactors(ny/2),zfactors(nz/2))
      call factor(nx,xfactors,fx)
      call factor(ny,yfactors,fy)
      call factor(nz,zfactors,fz)

!.....evaluating all possible decompositions............................
      allocate(solution(maxsol,3))
      i=0
      do n1=1,fx
        do n2=1,fy
          do n3=1,fz
            call check(valid, &
     &        nx,xfactors(n1),ny,yfactors(n2),nz,zfactors(n3),ncpu)

!.....saving valid decompositions......................................
            if (valid) then
              i=i+1
              solution(i,1) = xfactors(n1)
              solution(i,2) = yfactors(n2)
              solution(i,3) = zfactors(n3)
            endif
          enddo
        enddo
      enddo

!.....determine best decomposition.....................................
      allocate(value(i))
      do n=1,i
        call evaluate(value(n), &
     &    solution(n,1),solution(n,2),solution(n,3))
      enddo
      if (i.ge.1) then
        pos(1) = minloc(value,dim=1)
        nxmn = solution(pos(1),1)
        nymn = solution(pos(1),2)
        nzmn = solution(pos(1),3)
      else
        nxmn = 0
        nymn = 0
        nzmn = 0
      endif

      deallocate(xfactors,yfactors,zfactors,solution,value) 

      contains

!=======================================================================
!
!     subroutine evaluate determines the value of a decomposition
!
!=======================================================================

      subroutine evaluate(value,lx,ly,lz)

      integer, intent(out) :: value
      integer :: lx,ly,lz

!.....calculate the surface of the domain...............................
      value = lx*ly+ly*lz+lz*lx                 

      end subroutine evaluate

!=======================================================================
!
!     subroutine factor returns possible factors of n
!
!=======================================================================

      subroutine factor(n,factors,i)

      implicit none

      integer, intent(in) :: n
      integer, intent(out), dimension((n/2))::factors
      integer, intent(inout) :: i
      integer :: f

      i = 0
      do f=1,n
        if (mod(n,f) == 0) then
          i = i+1
          factors(i) = f
        endif
      enddo

      end subroutine factor

!=======================================================================
!
!     subroutine check determines the validity of a decomposition
!
!=======================================================================

      subroutine check(valid,nx,lx,ny,ly,nz,lz,ncpu)

      implicit none

      integer :: nx,lx,ny,ly,nz,lz,ncpu
      logical ,intent(out):: valid

      valid = ((mod(nx,lx).eq.0) .and. &
     &    (mod(ny,ly).eq.0) .and. &
     &    (mod(nz,lz).eq.0) .and. &
     &    ((nx/lx)*(ny/ly)*(nz/lz)-ncpu.eq.0))

      end subroutine check

!=======================================================================

      end subroutine decompose

!***********************************************************************
