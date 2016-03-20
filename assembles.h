#ifndef ASSEMBLES_H
#define ASSEMBLES_H

// libMesh include files.
#include <math.h> // needed for sqrt function in Coul()
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe_interface.h" // for dirichlet boundary conditions
#include "libmesh/error_vector.h" // for dirichlet boundary conditions
#include "libmesh/elem.h"
// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
#include <complex.h> // the infinite element version requires complex numbers explicitly.

//prototypes of functions in assemble.C (that are called from out of it)
double Harm(libMesh::Real , libMesh::Real, libMesh::Real);
double Coul(libMesh::Real , libMesh::Real, libMesh::Real);
double Morse(libMesh::Real, libMesh::Real, libMesh::Real);
void get_dirichlet_dofs(libMesh::EquationSystems& , const std::string& , std::set<unsigned int>&);
void assemble_EigenSE(libMesh::EquationSystems& , const std::string&);
void assemble_InfSE(libMesh::EquationSystems & es, const std::string & system_name);
// this one maybe will never work!?
void assemble_InfFullSE(libMesh::EquationSystems & es, const std::string & system_name);

#endif
