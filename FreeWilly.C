#include <assert.h> // 
#include <stdlib.h>    /* atof */
#include <vector>
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
//#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
//#include "libmesh/quadrature_gauss.h"
//#include "libmesh/dense_matrix.h"
//#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
//#include "libmesh/dof_map.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/fe_interface.h" // for dirichlet boundary conditions
#include "libmesh/error_vector.h" // for dirichlet boundary conditions
//#include "libmesh/explicit_system.h"
// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
#include <complex.h> // the infinite element version requires complex numbers explicitly.
// for refinement:
#include "libmesh/mesh_refinement.h"
#include "libmesh/kelly_error_estimator.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;
   
std::vector<Node> getGeometry(std::string);

void tetrahedralise_sphere(UnstructuredMesh& mesh, std::vector<Node> geometry, std::string creator, Real r, int NrBall, Real VolConst, Real L, unsigned int N);

int main (int argc, char** argv){
   // Initialize libMesh and the dependent libraries.
   LibMeshInit init (argc, argv);
   GetPot cl(argv[1]);

   #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
    // SLEPc currently gives us a nasty crash with Real==float
    libmesh_example_requires(false, "--disable-singleprecision");
   #endif
   // Check for proper usage.
   if (argc < 2)
      libmesh_error_msg("\nUsage: " << argv[0] << " <input-filename>");
   // Tell the user what we are doing.
   else {
      out<<"WELCOME TO FREE WILLY."<<std::endl;
      out << "Running " << argv[0];
      for (int i=1; i<argc; i++)
          out << " " << argv[i];
      out << std::endl << std::endl;
      out<<"Input file: "<<std::endl;
      out<<"--------------------------------------------------------------..."<<std::endl;
      std::string getcontent;
      std::ifstream openfile (argv[1]);
      if(openfile.is_open()){
        while(getline(openfile, getcontent)){
            out<<"| ";
            out << getcontent << std::endl;
        }
      }
      out<<"--------------------------------------------------------------..."<<std::endl;
      out<<"End of input file "<<std::endl<<std::endl;
   }
   
   int dim = 3;
   // Create a mesh, with dimension to be overridden later, on the
   // default MPI communicator.
   Mesh mesh(init.comm(), dim);
   
   std::string molec_file=cl("mol_file", "invalid_file"); // this file contains all informations on the molecule
   std::string angular_creator=cl("angular", "invalid"); 
   Real r=cl("radius", 20.);
   int NrBall=cl("points", 50);
   Real VolConst= cl("maxVol", 1./32 );
   Real L=cl("bending", 2.);
   int N=cl("circles", 5);
   std::vector<Node> geometry=getGeometry(molec_file);

   // the function below creates a mesh using the molecular structure.
   tetrahedralise_sphere(mesh, geometry, angular_creator, r, NrBall, VolConst, L, N);

   // All done.
   return 0;
}
