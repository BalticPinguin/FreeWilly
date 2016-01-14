#!/bin/bash

export SLEPC_DIR=/home/tm162/lib/libmesh/slepc-3.6.2
export PETSC_DIR=/home/tm162/lib/libmesh/petsc/
export PETSC_ARCH="arch-linux2-c-debug"

./configure --prefix=/home/tm162/bin/libmesh \
   --enable-default-comm-world \
   --disable-strict-lgpl   \  # enables tetgen and others.
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
   --enable-everything     >& configs.conf

make         >& make.log
make check   >& check.log
make -C tests check    >& unitcheck.log
make install >& install.log
