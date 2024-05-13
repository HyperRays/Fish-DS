LDLIBS= -I/soft/intel/mkl721/include -L/soft/intel/mkl721/lib/em64t \
				-I/soft/mpich-1.2.7p1/intel_v9/include \
				-qmkl -lpthread

# ---------- COMPILE FLAGS ----------
FC = mpiifort -fc=ifx -c                  # compile
LK = mpiifort -fc=ifx -o octree             # link

HA = -heap-arrays                 # stack overflow/limit workaround 
																	# Due to memory limitations even with unlimited
																	# stack size it is necessary to enable this flag
																	# to prevent segfault due to stack-overflow

OMPF = -qopenmp                   # OpenMP Flags
OF = -O2                          # Optimization Flags
SF = -mcmodel=medium -shared-intel -double-size 64 -fpp # Standard Flags
# DF = -g -debug -fp-stack-check -C -ftrapuv -warn interfaces -traceback # Debug Flags
CF = $(HA) $(OF) $(SF) $(DF) $(OMPF) # Combined Flags


# ------- END OF COMPILE FLAGS ------

OBJS = main.o octree.o 

octree: $(OBJS)
	@echo "Building: "$@
	@$(LK) $(CF) $(OBJS) $(LDLIBS)


main.o: octree.o main.f90
	@echo "Building: "$@
	@$(FC) $(CF) $(mhd)$(@:.o=.f90)

octree.o: octree.f90
	@echo "Building: "$@
	@$(FC) $(CF) $(mhd)$(@:.o=.f90)


clean: 
	@rm -f *.o
	@rm -f *.mod
	@rm -f octree