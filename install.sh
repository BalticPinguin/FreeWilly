export SLEPC_DIR='~/bin/libmesh/slepc'
export PETSC_DIR='~/bin/libmesh/petsc'
export PETSC_ARCH="arch-linux2-c-debug"

./configure --prefix=~/bin/libmesh \
   --enable-default-comm-world \
   --disable-strict-lgpl   \
   --disable-maintainer-mode \
   --enable-complex        \
   --enable-perflog        \
   --enable-qhull          \
   --enable-slepc          \
   --enable-petsc          \
   --enable-trilinos       \
   --enable-triangle       \
   --enable-tetgen         \
   --enable-tecplot        \
   --enable-ifem           \
   --disable-second        \
   --enable-everything     >& configs.conf

make         >& make.log
make check   >& check.log
make -C tests check    >& unitcheck.log
make install >& install.log
