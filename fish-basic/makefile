LDLIBS= -I/soft/intel/mkl721/include -L/soft/intel/mkl721/lib/em64t \
				-I/soft/mpich-1.2.7p1/intel_v9/include \
				-qmkl -lpthread

mhd = /app/fish-basic/source
dest = /app/fish-basic/dest

cngdir = cd $(dest)


# ---------- COMPILE FLAGS ----------
FC = mpiifort -c                  # compile
LK = mpiifort -o fish             # link

HA = -heap-arrays                 # stack overflow/limit workaround 
																	# Due to memory limitations even with unlimited
																	# stack size it is necessary to enable this flag
																	# to prevent segfault due to stack-overflow

OMPF = -qopenmp                   # OpenMP Flags
OF = -O2                          # Optimization Flags
SF = -mcmodel=medium -shared-intel -no-ip -double-size 64 -fpp # Standard Flags
# DF = -g -debug -fp-stack-check -C -ftrapuv -warn interfaces -traceback # Debug Flags
CUF = $(HA) $(OF) $(SF) $(DF) $(OMPF) # Combined Uncleaned Flags
CF = $(shell echo $(CUF) | tr -d '\t')  # Combined Flags


# ------- END OF COMPILE FLAGS ------

# --------  MAKEFILE ORDERS ---------
DEFAULT_RUN = flags fish EOMF

OBJS = \
	$(dest)/comm_module.o \
	$(dest)/decompose.o \
	$(dest)/driver.o \
	$(dest)/empty.o \
	$(dest)/eos.o \
	$(dest)/euler.o \
	$(dest)/init.o \
	$(dest)/mhd_module.o \
	$(dest)/mhd_parameter_module.o \
	$(dest)/prof_module.o \
	$(dest)/tvdb.o \
	$(dest)/tvdeul.o

# ----- END OF MAKEFILE ORDERS ------

# ---------------- BUILD ORDER RULES
.PHONY: default flags clean-all
default: $(DEFAULT_RUN)

flags: FORCE
	@echo "Using flags: $(CF) \n"


clean-all clean_all ca: clean_data clean_logs clean_out
	@echo "deleting executable"
	@rm -f fish

# ---------------- BUILD RULES
fish: $(OBJS)
	@echo "Building: "$@ 
	@$(LK) $(CF) $(OBJS) $(LDLIBS)

$(dest)/comm_module.o: $(mhd)/comm_module.f90 $(dest)/mhd_parameter_module.o 
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/decompose.o: $(mhd)/decompose.f90
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/driver.o: $(mhd)/driver.f $(dest)/mhd_module.o $(dest)/prof_module.o 
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/empty.o: $(mhd)/empty.f90 $(dest)/comm_module.o
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/eos.o: $(mhd)/eos.f90 $(dest)/mhd_parameter_module.o 
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/euler.o: $(mhd)/euler.f90 $(dest)/mhd_parameter_module.o $(dest)/prof_module.o $(dest)/comm_module.o
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/init.o: $(mhd)/init.f90 $(dest)/mhd_parameter_module.o $(dest)/comm_module.o
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/mhd_module.o: $(mhd)/mhd_module.f90 $(dest)/comm_module.o $(dest)/mhd_parameter_module.o $(dest)/prof_module.o \
							$(mhd)/calcfl.f90 $(mhd)/transpose.f90 $(mhd)/sweep.f90 \
							$(mhd)/readdmp.f90 $(mhd)/writdmp.f90 $(mhd)/output.f90 
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/mhd_parameter_module.o: $(mhd)/mhd_parameter_module.f90
	@echo $(@:.o=.f90)
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/prof_module.o: $(mhd)/prof_module.f
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/tvdb.o: $(mhd)/tvdb.f90
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

$(dest)/tvdeul.o: $(mhd)/tvdeul.f90 $(dest)/mhd_parameter_module.o $(dest)/mhd_module.o 
	@echo "Building: "$@ \($<\)
	@$(cngdir);$(FC) $(CF) $<

# ---- FILE CLEANUP & BASH SETTING RULES
.PHONY: FORCE EOMF clean_data clean_logs clean_out create_aux_folders
FORCE:                                  # force the execution of a rule

EOMF: create_aux_folders                               # rule to run at the end of make file (default flow)

# continue regardless of the success of command: || true
# supress output: 2>/dev/null

clean_data: FORCE
	@echo "cleaning data files"
	@rm -rf data/*


clean_logs: FORCE
	@echo "cleaning log files"
	@rm -rf logs/*

clean_out: FORCE
	@echo "cleaning out files"
	@rm -rf dest/*

create_aux_folders: FORCE
	@mkdir -p dest
	@mkdir -p logs
