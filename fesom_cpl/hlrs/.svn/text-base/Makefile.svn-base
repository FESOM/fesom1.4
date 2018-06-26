######################################################
# Makefile                
######################################################

######################################################
# Include architecture-specific definitions

include Makefile.in

###### Objects for Mesh Partitioning ################################################
# modules
MOD_INI =	oce_modules.o                 	 	\
		gen_modules_config.o			


# objects
OBJ_INI =   	gen_partitioning.o          		\
		fesom_partition_init.o \
		fort_part.o


# modules
MODULES =       oce_modules.o                           \
                ice_modules.o                           \
		gen_modules_config.o			\
                gen_modules_clock.o                     \
                gen_modules_rotate_grid.o               \
                gen_modules_read_NetCDF.o               \
                gen_modules_forcing.o                   \
                gen_modules_diagnose.o                  \
                oce_mixing_kpp.o                        \
                oce_mixing_pp.o                         \
                oce_mixing_my2p5.o                      \
                oce_mixing_tidal.o                      \
                oce_age_tracer.o                        \
                oce_passive_tracer.o  			\
		cpl_modules.o				\
                cpl_driver.o

# objects
OBJECTS=        oce_step.o                              \
                oce_init.o                              \
                oce_init_special_case.o                 \
                oce_dens_press.o                        \
                oce_neutral_physics.o                   \
                oce_mesh_setup.o                        \
                oce_read_mesh.o                         \
		tide_potential.o			\
                oce_matrices.o                          \
                oce_rhs_dyn.o                           \
                oce_rhs_tra.o                           \
                oce_fct.o                               \
                oce_convection.o                        \
                oce_open_boundary.o                     \
                oce_cavity.o                            \
                oce_tidal_boundary.o                    \
                oce_shortwave_pene.o                    \
                oce_salt_balance.o                      \
                oce_landice_water.o                     \
                oce_solver.o                            \
		oce_brine_rejection.o			\
		oce_mesh_aux.o				\
                ice_coupling.o                          \
                ice_step.o                              \
                ice_init.o                              \
                ice_thermo_oce.o                        \
                ice_thermo_cpl.o                        \
                ice_advection.o                         \
                ice_fct.o                               \
                ice_dyn.o                               \
                ice_solver.o                            \
	        psolve.o                                \
                gen_partitioning.o                      \
                gen_communication.o                     \
                gen_forcing_init.o                      \
                gen_forcing_ocean.o                     \
                gen_forcing_couple.o                    \
                gen_bulk_formulae.o                     \
                gen_input.o                             \
                gen_output.o                            \
                gen_calcmeans.o                         \
                gen_interpolation.o                     \
                gen_setup_model.o                       \
                gen_diag.o                              \
                gen_warning.o                           \
                gen_support.o                           \
		gen_output_mesh.o			\
                fesom_main.o          			\
		fort_part.o
######################################################



# Name of executables

EXE     = fesom.x
EXE_INI = fesom_ini.x


######################################################"
# TARGETS

default: fesom

fesom:  $(MODULES) $(OBJECTS) 
	@echo "======= Building FESOM =========="
	$(LD) $(OPT) -o $(EXE) $(FOS_INC) $(MODULES) $(OBJECTS) \
	 $(MPI_LIB) $(LIB_LAP) $(LIB_PARMS) $(LIB_METIS) $(OASIS3_LIB) $(NCLIB)


fesom_ini: cleanomod $(MOD_INI) $(OBJ_INI)
	@echo "======= Building FESOM paritioning program =========="
	$(LD)  $(OPT) -o $(EXE_INI) $(FOS_INC) $(MOD_INI) $(OBJ_INI) \
	$(MPI_LIB) $(LIB_LAP) $(LIB_PARMS) $(LIB_METIS) $(OASIS3_LIB) $(NCLIB)

.SUFFIXES: .c .F90 .o 

.c.o :
	$(CC) $(COPT) $(PARMS_DEFS) $(PARMS_INC) $(MPI_INCLUDE) -c $*.c	

.F90.o :
	@echo $@
	$(FC) $(OPT) $(CPP_SOL) $(CPP_DEFS) $(OASIS3_INC) $(FOS_INC)\
	$(PARMS_INC) $(MPI_INCLUDE) $(NCINC)  -c $*.F90

clean :
	rm -f *.o *.mod *~ *.f90 fesom fesom_ini

cleanomod:
	rm -f *.o *.mod
	CPP_DEFS=; export CPP_DEFS




