LDLIBS= -I/soft/intel/mkl721/include -L/soft/intel/mkl721/lib/em64t \
        -I/soft/mpich-1.2.7p1/intel_v9/include \
        -qmkl -lpthread

FC = mpiifort -c -heap-arrays                  #compile
LK = mpiifort -o fish              #link
SF = -mcmodel=medium -shared-intel -no-ip -double-size 64 -fpp -traceback    #compile flags
CF = $(SF) -O2               #optimization flags
# CF = $(SF) -g -debug -fp-stack-check -C -ftrapuv -warn interfaces #debug flags

mhd = source/

OBJS =\
  comm_module.o\
  decompose.o\
  driver.o\
  empty.o\
  eos.o\
  euler.o\
  init.o\
  mhd_module.o\
  mhd_parameter_module.o\
  prof_module.o\
  tvdb.o\
  tvdeul.o\

fish:\
  $(OBJS); $(LK) $(CF) $(OBJS) $(LDLIBS);\
  mv *.o dest;\
  mv *.mod dest;\
  mv *genmod* dest;\


comm_module.o:\
  mhd_parameter_module.o\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

decompose.o:\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

driver.o:\
  mhd_module.o\
  prof_module.o\
  $(mhd)$(@:.o=.f); $(FC) $(CF) $(mhd)$(@:.o=.f)

empty.o:\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

eos.o:\
  mhd_parameter_module.o\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

euler.o:\
  mhd_parameter_module.o\
  prof_module.o\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

init.o:\
  mhd_parameter_module.o\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

mhd_module.o:\
  comm_module.o\
  mhd_parameter_module.o\
  prof_module.o\
  $(mhd)calcfl.f90\
  $(mhd)transpose.f90\
  $(mhd)sweep.f90\
  $(mhd)readdmp.f90\
  $(mhd)writdmp.f90\
  $(mhd)output.f90\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

mhd_parameter_module.o:\
  $(@:.o=.f90); $(FC) $(CF) $(@:.o=.f90)

prof_module.o:\
  $(@:.o=.f); $(FC) $(CF) $(@:.o=.f)

tvdb.o:\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)

tvdeul.o:\
  mhd_parameter_module.o\
  mhd_module.o\
  $(mhd)$(@:.o=.f90); $(FC) $(CF) $(mhd)$(@:.o=.f90)
