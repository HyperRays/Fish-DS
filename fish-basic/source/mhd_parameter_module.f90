!***********************************************************************
!
!     parameter_module
!
!***********************************************************************

      module mhd_parameter_module

      implicit none

!.....choice of initial values for init.f90.............................
      integer, parameter :: mhd_problem=2

!.....solution vector...................................................
      integer, parameter :: nv=2 !number of vectors
      integer, parameter :: ns=2 !number of scalars
      integer, parameter :: ndbg=0 !number of debug scalars

!.....indizes of individual variables in solution vector................
      integer, parameter :: mv=1 !momentum vector
      integer, parameter :: mb=2 !magnetic field vector/sqrt(4pi)
      integer, parameter :: md=1 !density scalar
      integer, parameter :: ms=2 !entropy scalar

!-----------------------------------------------------------------------
!
!     Hint: if you want to track a single zone for debugging, this can
!     be difficult with all the rotations and buffers involved.
!     Set ns=5, ndbg=3 and write at the beginning the global zone
!     indices into s(3:5). The simple statement
!     if (trouble) write(6,*) 'I am zone ',nint(s(3:5))
!     then gives the correct information wherever it stands.
!
!-----------------------------------------------------------------------

!.....numerics..........................................................
      integer, parameter :: nxg=240 !number of zones along x direction
      integer, parameter :: nyg=240 !number of zones along y direction
      integer, parameter :: nzg=240 !number of zones along z direction
      integer, parameter :: nmax=240 !maximum n for a local cube
      integer, parameter :: yzbuf=4 !permanent xyz overlap between cubes
      integer, parameter :: buf=11 !permanent+volatile buffer length

!.....physics...........................................................
      real, parameter :: cflsav=0.75 !safety for cfl condition
      real, parameter :: dx=1.       !physical width of one zone
      real, parameter :: tf=20. !stop time

!.....file handling.....................................................
      character(12), parameter :: name='data' !name of run
      character(44), parameter :: &
     &  binpath='./data/' !path for binaries
      character(44), parameter :: &
     &  txtpath='./data/'  !path for ascii
      character(44), parameter :: &
     &  glbpath='./data/'   !path for restart files
     character(20), parameter :: &
     &  logdir='./logs/'   !path for logfiles
      integer, parameter :: skip=5    !skip for bin and txt output 5
      integer, parameter :: dmpskip=5 !skip for restart files 5
      integer :: restart=0 !where to restart

!-----------------------------------------------------------------------
!
!     Note: if after data output a file named "stop" is detected
!     in txtpath, the code dumps restart files and exits nicely.
!
!     Note: restart=0 starts the run from scratch, restart=x restarts
!     the run from restartfiles name_xxxxx_process.dmp. To automatically
!     delete all restart files of previous runs with the same name and
!     txtpath, create a file named "delete" in txtpath and start a new
!     run from scratch. Thus, don't edit/remove the logfiles
!     txtpath/name_process.dmp by hand unless you know what you do!
!
!-----------------------------------------------------------------------

      contains

!=======================================================================
!
!     Reading number of last restart file from last_restart.file
!
!=======================================================================
 
        subroutine read_restart_loc
      
        implicit none
        logical :: found
        integer :: res

        inquire(file=trim(txtpath)//'last_restart.file',exist=found)
        if (found) then
          open(unit=1,file=trim(txtpath)//'last_restart.file', &
     &      status='old')
          read (1,*) res
          write(6,*) 'Reading restart location ',res
          close(1)
          restart = res
        endif
      
        end subroutine read_restart_loc

!=======================================================================
!
!     Writing number of last restart file to last_restart.file
!
!=======================================================================
      
         subroutine write_restart_loc(number)
      
         implicit none
         integer, intent(in) :: number
 
         open(unit=1,file=trim(txtpath)//'last_restart.file', &
     &     status='unknown')
         write(6,*) 'Writing restart location ',number
         write (1,*) number
         close(1)

         end subroutine write_restart_loc

!=======================================================================

      end module mhd_parameter_module

!***********************************************************************
