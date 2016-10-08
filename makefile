######################################################################
#
# Template libMesh application Makefile
LIBMESH_DIR ?= /home/tm162/bin/libmesh/my_libmesh
#LIBMESH_DIR ?= /home/tm162/bin/libmesh/ifem_map_test


# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

#target	   := ./InSpect-$(METHOD)
target	   := ./FrWll-$(METHOD)
#target	   := ./test-$(METHOD)

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#srcfiles	:= FreeWilly.C assembles.C Mesh.C radial_interpolation.cpp
srcfiles	:= FreeWilly.C assembles.C Mesh.C radial_interpolation.C NN_interpolation.C fsu_soft/r8lib.C fsu_soft/rbf_interp_nd.C read_DO.C bas_pars_library/build/libbas_pars.a normalisation.C fsu_soft/legendre_polynomial.cpp Cube_IO.C fsu_soft/sphere_lebedev_rule.cpp
#   read_DO.C bas_pars_library/build/libbas_pars.a
#srcfiles	:= test.C radial_interpolation.C r8lib.C rbf_interp_nd.C

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

run: complete

complete: $(wildcard *.in)
#	@$(MAKE) dust
	@$(MAKE) -C $(dir $(target)) $(notdir $(target))
	@echo "***************************************************************"
	@echo "* Running App " $(notdir $(target))
	@echo "***************************************************************"
	@echo " "
	${LIBMESH_RUN} $(target) ${LIBMESH_OPTIONS} 2>&1 | tee output.txt
	@bzip2 -f output.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running App " $(notdir $(target))
	@echo "***************************************************************"

# include the dependency list
-include .depend

#
# Dependencies
#
.depend: $(srcfiles) $(LIBMESH_DIR)/include/libmesh/*.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(LIBMESH_DIR)/include $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) -I bas_pars_library/build -lgfortran > .depend

###############################################################################
