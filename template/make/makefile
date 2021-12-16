# work_standard_makefile

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

#################################################################

ifndef STAR
STAR = star
endif

# STAR_WORK_OBJS defined in makefile_header
OBJS = $(STAR_WORK_OBJS) run.o diffusivities.o utilities.o identification.o magnetism.o mdot.o penetration.o structure.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job
WARN = -Wall -Wdo-subscript -Wzerotrip -Wtarget-lifetime -Wcompare-reals -Wrealloc-lhs-all -Wrealloc-lhs -Wfunction-elimination -Wintrinsic-shadow -Wundefined-do-loop -Wsurprising -Wintrinsics-std -Winteger-division

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header



ifeq ($(SKIP_NAN_TRAPS),YES)
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<
else
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $(FCtrapNANs) $<
endif

diffusivities.o: $(WORK_SRC_DIR)/diffusivities.f90
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/diffusivities.f90

utilities.o: $(WORK_SRC_DIR)/utilities.f90
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/utilities.f90

identification.o: $(WORK_SRC_DIR)/identification.f90
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/identification.f90

magnetism.o: $(WORK_SRC_DIR)/magnetism.f90 diffusivities.o
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/magnetism.f90

penetration.o: $(WORK_SRC_DIR)/penetration.f90
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/penetration.f90

structure.o: $(WORK_SRC_DIR)/structure.f90
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/structure.f90


mdot.o: $(WORK_SRC_DIR)/mdot.f90
	$(WORK_COMPILE) $(WARN) $(WORK_SRC_DIR)/mdot.f90

run_star_extras.o: diffusivities.o utilities.o identification.o magnetism.o mdot.o penetration.o structure.o $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) $(WORK_SRC_DIR)/run_star_extras.f90

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)









