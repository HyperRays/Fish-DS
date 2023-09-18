LDLIBS= -I/soft/intel/mkl721/include -L/soft/intel/mkl721/lib/em64t \
        -I/soft/mpich-1.2.7p1/intel_v9/include \
        -qmkl -lpthread

mhd = source/


# ---------- COMPILE FLAGS ----------

FC = mpiifort -c                  #compile
LK = mpiifort -o fish             #link

HA = -heap-arrays                 # stack overflow/limit workaround
OF = -O2                          # Optimization Flags

SF = -mcmodel=medium -shared-intel -no-ip -double-size 64 -fpp #Standard Flags

DF = -g -debug -fp-stack-check -C -ftrapuv -warn interfaces -traceback #Debug Flags

CF = $(HA) $(OF) $(SF) $(SF) $(DF) # Combined Flags

# ------- END OF COMPILE FLAGS ------


# --------  MAKEFILE ORDERS ---------

DEFAULT_RUN = \
  fish\
  EOMF\

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

# ------ ENF OF MAKEFILE ORDERS -----



# ---------------- BUILD ORDER RULES

default:\
  $(DEFAULT_RUN)

clean-all:\                              # quality-of-life clean rule
  clean_logs\
  clean_data\

# ---------------- BUILD RULES

fish:\
  $(OBJS); $(LK) $(CF) $(OBJS) $(LDLIBS);\

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


# ---- FILE CLEANUP & BASH SETTING RULES

FORCE:                                  # force the execution of a rule

EOMF:\                                  # rule to run at end of make file (default flow)
  set_stack_limit\
  cp_out_files\
  create_aux_folders\

set_stack_limit:\
  FORCE;\
  ulimit -s unlimited;\

cp_out_files:\                          # copies all obj, mod, etc. files to dest folder
  FORCE;\
  mv *.o dest;\
  mv *.s dest;\
  mv *.mod dest;\
  mv *genmod* dest;\
  true                                  # continue regardless of success of previous commands


clean_data: FORCE; rm -rf data/*; true
clean_logs: FORCE; rm -rf logs/*; true

create_aux_folders:\
  FORCE;\
  mkdir dest;\
  mkdir logs;\
  true                                  # continue regardless of success of previous commands
  



