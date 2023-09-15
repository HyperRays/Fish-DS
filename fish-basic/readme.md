------------------------------------------------------------------------
      A similar version of this code is documented in
      "FISH: A 3D parallel MHD code for astrophysical applications"
      R. Kaeppeli, S. C. Whitehouse, S. Scheidegger,
      U.-L. Pen, M. Liebendoerfer, http://arxiv.org/abs/0910.2854
------------------------------------------------------------------------

Some more informal introductionary comments for new students
************************************************************

Important variables
===================

All scalar variables (density, energy density, debug indices)
are contained in the array s(ns,i,j,k) and all vector variables
momentum, magnetic field) are contained in the array
v(3,nv,i,j,k). The file mhd_parameter_module.f90 specifies parameters
to identify the physical quantities within the array. For example, the
first array element s(1,i,j,k) is the density and s(2,i,j,k) the internal
energy density. Let's assume I now code everything in terms of s(1,i,j,k)
and s(2,i,j,k). So far ok. But let's assume I would like to introduce a
two-component fluid with density1 and density2 of two components, so I
would have to change all places that refer to the array s and am likely
to forget some and introduce errors.  To avoid this, I define once in
mhd_parameter_module.f90 md=1 and ms=2 and then refer to s(md,i,j,k)
for density and s(ms,i,j,k) for energy density. Now it is very simple
to extend the code for two components: md1=1, md2=2, ms=3. With this
approach I don't have to change anything with respect to the references
to energy density and if I forget to update some reference to s(md,i,j,k)
the compiler will remind me that md is not defined.

For large simulations it is not possible to keep the full array in the
memory of one node. Therefore it is necessary to divide (decompose)
the array into different sections where each node can hold one section
and work on it. Of course it has to exchange data occasionally with
neighboring nodes. Our MHD code works with cubic domain decomposition,
i.e. the computational domain can be cut in each direction so that the
computational domain on one node becomes again more or less cubic (not
to be taken too strictly, more to be understood in contrast to slabs
or bars).  So each node contains a block of data and all blocks together
form the original array s or v. One node can only update the variables on
its full block if it also knows the content of the neighboring zones on
a neighboring block (e.g. the pressure gradient at the block edge depends
on the pressure on the neighbor block).  So each block has to be enveloped
by some buffer zones that mirror the data from the neigboring blocks. The
decomposition into blocks with buffer zones padded on each side makes
the indexing a bit intricate because one has to know the local index
in the decomposed array as well as the global array index to identify
the physical position of the zone.  I explain the index conventions in
x-direction, they are exactly the same in the y-direction and z-direction:

The indexing is handled by a data structure named comm_wld that is defined
in comm_module.f90:
```
      type comm_wld
        integer :: d !global label of direction {x=1,y=2,z=3}
        integer :: g !global number of zones in that direction
        integer :: r !global index of index 0 in local array
        integer :: m,n !start and end of local section without buffers
        integer :: l !dimension of array inclusive buffers
        integer, dimension(4) :: requests !communication handles
      end type comm_wld

For example in x-direction (*** indicates zones in a block of the
computaitonal domain while ... indicates buffer zones):

         node 1:
         d=1, g=48, r=-8, l=32, m=9, n=24
         ........****************........
         |       |              |       |
local:   1       m              n       l
global: r+1     r+m            r+n     r+l

                         node 2:
                         d=1, g=48, r= 8, l=32, m=9, n=24
                         ........****************........
                         |       |              |       |
                local:   1       m              n       l
                global: r+1     r+m            r+n     r+l

                                         node 3:
                                         d=1, g=48, r=24, l=32, m=9, n=24
                                         ........****************........
                                         |       |              |       |
                                local:   1       m              n       l
                                global: r+1     r+m            r+n     r+l

```
These conventions apply to all directions and are stored in variables
nx,ny,nz, each being of above type comm_wld. For example nx%g is the
total number of zones in x-direction, ny%l is the local array length
of a decomposed block including buffer zones in y-direction, or
nz%r+nz%m is the global index of the first local zone that is not a buffer.
If for example the x-direction and y-direction are swapped in transpose.f90
(see below), it is sufficient to also swap nx with ny and all dimensions
and indexing is again consistent.

By swapping directions in the data it is possible to write all routines
only in the x-direction. The basic approach is to rotate the data and
apply updates always across the current x-direction instead of keeping the
data fixed and coding updates for the x-, y-, and z-direction separately.
This makes the coding of subroutines very simple, but complicates the
indexing and debugging (brain gymnastics required to identify the global
zone index if a zone causes trouble in a rotated state on a local block).

Source code
===========

main directory
--------------
This directory contains the mhd_parameter_module.f90. It contains
all the settings and initial values for the run. You have to
recompile the whole code to see an effect of changes in
mhd_parameter_module.f90!

mhd_parameter_module.f90:

- solution vector
The code stores all variables on a Cartesian grid. Scalars are
stored in the array s(ns,:,:,:), where ns specifies how many
scalar quantities there are. Vectors are stored in the array
v(3,nv,:,:,:), where nv specifies how many vector quantities
there are. Each vector has three components. Parameters in the
mhd_parameter_module.f90 define which component in the array
belongs to which physical quantity. For example:
s(md,i,j,k) is the density at grid point x(i),y(j),z(k)
v(1,mv,i,j,k) is the x-component of the momentum at above grid point
v(3,mb,i,j,k) is the z-component of the magnetic field ibidem etc.
Sometimes the index 'ms' is used for entropy, sometimes for
energy density. In the example given with fish, it is only
used for energy density (even if the labels mostly say 'entropy').

- numerics
The computational domain is decomposed into blocks to be treated
by different processors in parallel. 'nmax' is relevant for the
memory consumption of each process, it must be set to the maximum
dimension of the decomposed blocks. Here it is set to the overall
maximum lengths nxg,nyg,nzg to allow also to run the code on one
processor without any decomposition.

- physics
The code is written in natural units. All quantities are expressed
in units of natural constants so that the constants do not appear
anymore in the code. However, one length scale of the code must be
set, it is set in 'dx' which specifies the physical width of the
computational grid (distance between neighboring grid points).

- file handling
Each run has a name that is specified here. Then there is a
path for the binary output files with array segments (some large
disk) and a different path for the ascii files (require less space
and preferentially somewhere where they are backed up). Then there
is a local path for the restart files. Because the restart files
are very large they are best written onto local scratch disks. In
this example, however, the files are not yet large and they can all
be written into the same directory. 'skip' defines how many time
steps should be done between file dumps and 'dmpskip' defines how
many dumps should be done between saved restart files. If the
parameter 'restart' is zero, the run will start from the beginning,
if 'restart' is set to some number, the run will try to restart
from the restart file with that number.

prof_module.f:
This routine measures the execution time of code-blocks and is used
to otpimize the code.

source
------

calcfl.f90:
Calculates the Courant-Friedrichs-Levy condition for the time step.

comm_module.f90:
Setup of all MPI communications between the decomposed blocks. The
connection between the buffer zones on different processes is installed
permanently and named by handles. Whenever one or more of these handles
are passed to MPI_STARTALL, the corresponding data is exchanged between
the processes.

driver.f:
Main program.

empty.f90:
Does nothing. See sweep.f90 to learn about its purpous...

eos.f90:
Polytropic equation of state, calculates the pressure and sound speed
as a function of density and internal energy.

euler.f90: 
Manages the update of the fluid variables and the magnetic field due
to one time step update in x-direction. It calls tvdeul.f90 and tvdb.f90
for the details of the one-dimensional update. This is the main (this
version 'only') place where OpenMP is activated for a the hybrid
execution of the code.

init.f90:
Provides the initial data if the run is not restarted from a previously
generated restart file.

mhd_module.f90:
Module that encapsulates the interior working of the MHD code and
provides simple interfaces for external calls. This avoids doubly-
declared variables and other inconveniences if the code is to be
merged with other code parts.

output.f90:
Creates the output files. This is also the location where it is
determined which information goes into the files
name_stepnumber_processnumber.dat (currently the three main
plains).

readdmp.f90:
Manages the restart from restartfiles.

sweep.f90:
One of the most important routines. Physically it manages the update
of variables in x-direction. Technically it receives the name of a function
as argument. Then it triggers the exchange of buffer zones and, at the
same time, calls the function to work on the inner array-section that
is not affected by buffer zones. When this work is done and the buffer
zones have been exchanged, it completes the update by calling the
function again for an update of the border zones that are affected
by the buffer zones. If for testing one wants to just exchange the
buffer zones without doing anything, one can call sweep with the
argument 'empty'.

transpose.f90:
Swaps directions in the 3D computational domain. Instead of programming
all routines for the x-, y-, and z-direction, the routines are only
written for the x-direction. Then the directions in the data are swapped
between sweeps, so that the updates proceed in the following order:
x-update,y-update,z-update,z-update,y-update,x-update. The ordering is
symmetric to maintain second order accuracy in time.

tvdb.f90:
Does the details of one-dimensional update of the magnetic field
in x-direction.

tvdeul.f90:
Does the details of one-dimensional hydrodynamics step in
x-direction. This routine is critical for speed and should be
optimised...

writdmp.f90
Writes restart files.

Compile & Run
=============

go to the directory with the makefile and type 'make' at the prompt. This
should compile the code and result in an executable fish. There are
also many object files *.o and *.mod, they are not used for running
the program (but can be used to compile the next time faster). However,
if the parameters in mhd_parameter_module.f90 have been changed, it is
a good idea to rebuild the whole code with 'make -B'.

Then run the code by typing 'mpirun ./fish'. It will run in one process.
If you want to run it in parallel then assign more CPU's with the flags
listed in the batch script qfish.

Results
=======

Output files
------------
The most important output files contain the name of the run,
the time step number, and the process number from the cubic
domain decomposition. Output files about performance and
logfiles of the progress are written to the main directory.
Datafiles are written to the path given in mhd_parameter_module.f90,
which is currently set to ./data. The following outputfiles are generated:

- name_stepnumber_processnumber.dat
This file contains the most useful data, i.e. array-segments
out of the full volumetric data of all variables in the s and
v arrays. The segments are defined in the output.f90 routine.
Currently these are the three main planes at the center of the
computational domain. The data is in binary format to save
space on the disk.

- scratch_processnumber.log
This is the logfile for process processnumber. Info's, Warnings,
Errors that occur during the run are written there. The logfile
of process zero is written to standard output, not to file
scratch_0000.dat. If something goes wrong, it sometimes help
to look at the last messages in scratch_*.log

- name_stepnumber_processnumber.dmp
This file contains the full arrays s and v with some additional
information so that runs can be restarted from these files. Due
to the large size of this file, it is not kept for all time steps.
Parameters in mhd_parameter_module.f90 specify how often the
restart files are kept.

- name_processnumber.dmp
This is an important file if the run runs on many processors and the
data is written to local discs. It saves a record about which data has
been saved in restart files on which node (all restart files are kept
locally on the nodes).  If one wants to restart a run, each new process
has to look up in this file, on which node its precessor has run and
then read the corresponding restart file from that node.

- profile_processnumber.dat
This file contains performance information. It is helpful to
identify the slowest parts in the program (respectively the
parts where it is most beneficial to think about optimization).

Visualization
-------------

Already available are some Matlab scripts (functions) to visualize the
results. If you make a directory 'matlab' in your home-directory and put
the scripts directly into it, matlab will automatically find them. There
should be following functions:

- drawslice.m

draws slices through the computational domain. With the current settings
in output.f90, segment 1 is the xy-plane, segment 2 the xz-plane and
segment 3 the yz-plane. Four windows are opened with the density,
specific energy, velocity, and magnetic field. drawslice is based on
the data contained in the binary files name_stepnumber_processnumber.dat

- readslice.m
reads the binary files name_stepnumber_processnumber.dat for
all processnumbers of the run and merges the domain-decomposed
data to one big array 'u' that can be used in matlab to plot
or analyse the data.

Usage:
Start matlab e.g. with 'matlab -nojvm' from the directory where the output
data are (or from anywhere and then use 'cd' within matlab to change
the directory to where the data are). At the matlab prompt, use e.g.

drawslice for slices through the main planes:
```
>> help drawslice
>> drawslice('data',0,1,1e6*[-1 1 -1 1 -1 1]);
>> drawslice('data',0,1,5e8*[-1 1 -1 1 -1 1])
```