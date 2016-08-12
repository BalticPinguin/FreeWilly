##############
# compiler:
# currently working: gcc 6.1.0 and mkl
# also working: gcc 4.9.3 and acml 5.3.1
FF = gfortran
##############
# development flags as recommended by http://www.fortran90.org/src/faq.html
# This warns about undefined symbols, stops at the first error, 
# turns on all debugging checks (bounds checks, array temporaries, ...) 
# and turns on backtrace printing when something fails at runtime (typically accessing an array out of bounds). 
# You can use -Werror to turn warnings into errors (so that the compilation stops when undefined symbol is used).
# With gfortran 4.4 and older, replace fcheck=all with -fbounds-check -fcheck-array-temporaries
FFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace 
##############
# build directory:
builddir = build
##############
# install directory
installdir = bin
##############
# library directory:
libdir = build
# object files to be compiled
# ordered by dependency, the latter depends on at least one of the former modules
constants = constants.o
global    = global_data.o global_utils.o
#intlib    = libCint_interface.o
parseInp  = parse_utils.o basis_utils.o parse_inp.o
utils     = utils.o
main      = Auger.o
lib       = libbas_pars.a
#
obj = $(constants) $(utils) $(global) $(intlib) $(parseInp)
buildObj = $(patsubst %.o, $(builddir)/%.o, $(obj))
##############
# Linker flags
#LDFLAGS = -Wl,-rpath,/home/gg114/lib/GTO_int_lib/LibCInt/libcint-master-23.03.2016/build_int64  -L/home/gg114/lib/GTO_int_lib/LibCInt/libcint-master-23.03.2016/build_int64 -lcint 
#LDFLAGS = -Wl,-rpath,/home/gg114/lib/GTO_int_lib/LibCInt/libcint_int64/build_int64  -L/home/gg114/lib/GTO_int_lib/LibCInt/libcint_int64/build_int64 -lcint 
#LDFLAGS = -Wl,-rpath,/scratch/tm162/bin/FreeWilly/libcint/bin/lib/x86_64-linux-gnu -L/scratch/tm162/bin/FreeWilly/libcint/bin/lib/x86_64-linux-gnu -lcint -Wl,-rpath,/scratch/tm162/bin/FreeWilly/bas_pars_library/$(libdir) -L/scratch/tm162/bin/FreeWilly/bas_pars_library/$(libdir) -lbas_pars
LDFLAGS = -Wl,-rpath,/scratch/tm162/bin/FreeWilly/bas_pars_library/$(libdir) -L/scratch/tm162/bin/FreeWilly/bas_pars_library/$(libdir) -lbas_pars
##############
# compilation block
# link the main program:
$(builddir)/auger.exe :  $(libdir)/$(lib) $(builddir)/$(main)
	$(FF) -fopenmp $(LDFLAGS) -I $(builddir) -o $@ $(builddir)/$(main) $(libdir)/$(lib)
#
# create the library:
$(libdir)/$(lib) : $(buildObj)
	ar rcs $@ $^
# compile files in the /constants dir
constantdir = src/constants
$(builddir)/%.o: $(constantdir)/%.f95
	$(FF) $(FFLAGS) -c $< -o $@ -J $(builddir)
#
# compile files in the /utils dir
utildir = src/utils
$(builddir)/%.o: $(utildir)/%.f95
	$(FF) $(FFLAGS) -c $< -o $@ -J $(builddir)
#
# compile files in the /global_data dir
globdir = src/global_data
$(builddir)/%.o: $(globdir)/%.f95
	$(FF) $(FFLAGS) -c $< -o $@ -J $(builddir)
#
# compile files in the /libCint_interface dir
#libCdir = src/libCint_interface
$(builddir)/%.o: $(libCdir)/%.f95
	$(FF) $(FFLAGS) -c $< -o $@ -J $(builddir)
#
# compile files in the /parse_inp dir
parsedir = src/parse_inp
$(builddir)/%.o: $(parsedir)/%.f95
	$(FF) $(FFLAGS) -c $< -o $@ -J $(builddir)
#
# compile files in the /main dir
maindir = src/main
$(builddir)/%.o: $(maindir)/%.f95 $(buildobj)
	$(FF) $(FFLAGS) -c $< -o $@ -J $(builddir)
#
.PHONY : install
install:
	mv $(builddir)/*.exe $(installdir)

.PHONY : clean
clean:
	rm $(builddir)/*.o $(builddir)/*.mod ; rm $(installdir)/*.exe ; rm $(libdir)/*.a

