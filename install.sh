#!/bin/bash

# built without nodeconstraint
#  with infinite elements, without mpi
#  use tripleprecision scalars
#  Provide global libMesh::CommWorld
#  build without periodic boundary condition suppport
#  build to support complex-number solutions
#  build with performance logging turned on
#  build with slepc library...
#  build with Qhull API support

   #--disable-nodeconstraint \

CXXFLAGS=-std=c++11 # try to get rid of cast-errors
export PETSC_DIR=/home/hubert/bin/libmesh/petsc-3.6.2 
export PETSC_ARCH=arch-linux2-c-debug
export SLEPC_DIR=/home/hubert/bin/libmesh/slepc-3.6.2
./configure --prefix=/usr/local \
   --disable-mpi            \
   --enable-tripleprecision \
   --enable-default-comm-world \
   --disable-periodic      \
   --enable-complex        \
   --enable-perflog        \
   --enable-qhull          \
   --enable-slepc          \
   --enable-petsc          \
   --enable-tetgen         \
   --enable-triangle       \
   --enable-vtk            \
   --enable-trilinos       \
   --enable-ifem            >& configs.conf

make         >& make.log
make check   >& check.log
make -C tests check    >& unitcheck.log
#make install >& install.log
