#ifndef FREEWILLY_H
#define FREEWILLY_H

#include <assert.h> // 
#include <stdlib.h>    /* atof */
#include <vector>
#include <math.h> // needed for sqrt function
#include <iostream> 
#include <complex.h> // the infinite element version requires complex numbers explicitly.
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
//for mesh generation
#include "libmesh/face_tri3.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_tetgen_interface.h"
#include "libmesh/mesh_triangle_interface.h"
#include "libmesh/node.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/fe_interface.h" // for dirichlet boundary conditions
#include "libmesh/error_vector.h" // for dirichlet boundary conditions

// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
// for refinement:
#include "libmesh/mesh_refinement.h"
#include "libmesh/kelly_error_estimator.h"
//for interpolation
#include "libmesh/mesh_function.h"
#include "libmesh/meshfree_interpolation.h"
#include "radial_interpolation.h"
#include "NN_interpolation.h"
// for the ESP system
#include "libmesh/explicit_system.h"
// for the lebedev/womersley-grids
# include "fsu_soft/sphere_lebedev_rule.hpp"
# include "fsu_soft/sphere_design_rule.hpp"
# include "grids/Wom.h"
# include "grids/geodesic.h"

#include "read_DO.h"


// own enumeration type for integral types.
enum IntegralType: int{
   LENGTH=0,
   VELOCITY,
   OVERLAP
};

//Function that reads the geometry and charge of the molecule
std::vector<libMesh::Node> getGeometry(std::string filename);

//prototypes of functions in assemble.C 
void get_dirichlet_dofs(libMesh::EquationSystems& , const std::string& , std::set<unsigned int>&);
void assemble_InfSE(libMesh::EquationSystems & es, const std::string & system_name);
void assemble_ESP(libMesh::EquationSystems & es, const std::string & system_name);
void assemble_DO(libMesh::EquationSystems & es, const std::string & system_name);

// self-written output formats:
void cube_io(libMesh::EquationSystems& es, std::vector<libMesh::Node> geom, 
             std::string output, std::string SysName);
void line_out(libMesh::EquationSystems& es, std::string output, std::string SysName);

//This in the tetrahedralisation of a sphere
void tetrahedralise_sphere(libMesh::UnstructuredMesh& mesh, 
                           std::vector<libMesh::Node> geometry, 
                           std::string creator, 
                           libMesh::Real r, 
                           std::string scheme, 
                           libMesh::Real p, 
                           libMesh::Real VolConst, 
                           libMesh::Real L,
                           unsigned int N);

// Helper routine for tetrahedralize_domain().  Adds the points and elements
// of a convex hull generated by TetGen to the input mesh
void add_sphere_convex_hull_to_mesh(libMesh::MeshBase& mesh, 
                                    libMesh::Real r_max, 
                                    std::string scheme, 
                                    libMesh::Real p, 
                                    std::vector<libMesh::Node> geometry, 
                                    std::string creator, 
                                    const libMesh::Real L, 
                                    const unsigned int N);

unsigned int point_size(unsigned int iterations, int n);

// Functions that read and evaluate the Dyson orbital
void getDyson(const char *filename, 
              int namelength, 
              std::vector<std::vector<double> >& do_j, 
              std::vector<unsigned int>& l,
              std::vector<double>& alpha, 
              double&  energy, 
              double& normDO);

double evalDO(const std::vector<std::vector<double> >& , 
              const std::vector<unsigned int>&, 
              const std::vector<double>& , 
              const std::vector<libMesh::Node>& ,
              const libMesh::Point);

// functions to calculate norms and overlaps:
libMesh::Real overlap_DO(libMesh::EquationSystems& eq_sys,
                           const std::string sys1, int var1, IntegralType int_type, bool infinite=false);
libMesh::Number norm_DO(libMesh::EquationSystems& , bool );
libMesh::Real normalise(libMesh::EquationSystems& , bool infel=false);
libMesh::Number calculate_overlap(libMesh::EquationSystems& eq_sys, 
                                  const std::string sys1, int var1, 
                                  const std::string sys2, int var2,
                                  IntegralType int_type);

void ProjectSphericals (libMesh::EquationSystems&, int ,int);

#endif //define FREEWILLY_H
