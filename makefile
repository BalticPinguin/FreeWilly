######################################################################
#
# Template libMesh application Makefile
LIBMESH_DIR ?= /home/tm162/bin/libmesh/my_libmesh
#LIBMESH_DIR ?= /home/tm162/bin/libmesh/ifem_map_test


# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

target	   := ./FrWll-$(METHOD)

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

this_files	:= FreeWilly.C assembles.C Mesh.C radial_interpolation.C NN_interpolation.C read_DO.C normalisation.C Cube_IO.C SlepcConfig.C plotSH.C
fsu_files       := fsu_soft/r8lib.C fsu_soft/rbf_interp_nd.C fsu_soft/legendre_polynomial.C fsu_soft/sphere_lebedev_rule.C fsu_soft/sphere_design_rule.C grids/geodesic.C fsu_soft/fn_prb.C fsu_soft/besselj.C
grd_files       := grids/Wom_gen.C grids/Wom_ev.C grids/Wom_md.C grids/Wom_me.C grids/Wom_mn.C 
srcfiles        := $(this_files) $(fsu_files) $(grd_files) bas_pars_library/build/libbas_pars.a

objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

.PHONY: dust clean distclean

###############################################################################
# Target:
#

all:: $(notdir $(target))

# Production rules:  how to make the target - depends on library configuration
$(notdir $(target)): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS)


# Useful rules.
dust:
	@echo "Deleting old output and runtime files"
	@rm -f out*.m job_output.txt output.txt* *.gmv.* *.plt.* out*.xdr* out*.xda* PI* complete

clean: dust
	@rm -f $(objects) *.$(obj-suffix) *.lo

clobber: clean 
	@rm -f $(target)

distclean: clean
	@rm -rf *.o .libs .depend

echo:
	@echo srcfiles = $(srcfiles)
	@echo objects = $(objects)
	@echo target = $(target)

# include the dependency list
-include .depend

#
# Dependencies
#
.depend: $(srcfiles) $(LIBMESH_DIR)/include/libmesh/*.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(LIBMESH_DIR)/include $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) -I bas_pars_library/build -lgfortran > .depend

###############################################################################
