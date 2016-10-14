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

// Function prototype.  This is the function that will assemble
// the eigen system. Here, we will simply assemble a mass matrix.
double normalise(EquationSystems& equation_systems, bool infel);
std::vector<libMesh::Node> getGeometry(std::string filename);

//prototypes of functions in assemble.C (that are called from out of it)
void get_dirichlet_dofs(libMesh::EquationSystems& , const std::string& , std::set<unsigned int>&);
void assemble_InfSE(libMesh::EquationSystems & es, const std::string & system_name);
void assemble_ESP(EquationSystems & es, const std::string & system_name);
void assemble_DO(EquationSystems & es, const std::string & system_name);
void cube_io(EquationSystems& es, std::vector<Node> geom, std::string output, std::string SysName);

//This in the tetrahedralisation of a sphere
//void tetrahedralise_sphere(libMesh::MeshBase& mesh, const libMesh::Parallel::Communicator& comm);
// void tetrahedralise_sphere(mesh, const libMesh::Parallel::Communicator& comm);
//void tetrahedralise_sphere(libMesh::UnstructuredMesh& mesh, std::vector<Point> geometry, std::string creator);
void tetrahedralise_sphere(UnstructuredMesh& mesh, std::vector<Node> geometry, std::string creator, Real r, int NrBall, Real VolConst, Real L, unsigned int N);

enum IntegralType: int{
   MU=0,
   OVERLAP
};

Number overlap_DO(EquationSystems& eq_sys, const std::string sys1, int var1, IntegralType int_type);
Real normalise(EquationSystems& equation_systems, bool infel);
Number calculate_overlap(EquationSystems& eq_sys, const std::string sys1, int var1, const std::string sys2, int var2 , IntegralType int_type);

int main (int argc, char** argv){
   // Initialize libMesh and the dependent libraries.
   LibMeshInit init (argc, argv);
   GetPot cl(argv[1]);

   // Skip SLEPc examples on a non-SLEPc libMesh build
   #ifndef LIBMESH_HAVE_SLEPC
     libmesh_example_requires(false, "--enable-slepc");
     }
   #else
  
   #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
    // SLEPc currently gives us a nasty crash with Real==float
    libmesh_example_requires(false, "--disable-singleprecision");
   #endif
   bool infel = cl("infinite", false);
   if (infel){
      #ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
      libmesh_example_requires(false, "--enable-ifem");
      #endif
    }
   
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

   // Get the number of eigen values to be computed from argv[2]
   //const unsigned int nev = std::atoi(argv[2]);
   const unsigned int nev = cl("nev",10);

   // Skip this 2D example if libMesh was compiled as 1D-only.
   libmesh_example_requires(3 <= LIBMESH_DIM, "2D support");
   
   int dim = 3;
   // Create a mesh, with dimension to be overridden later, on the
   // default MPI communicator.
   Mesh mesh(init.comm(), dim);
   
   Real E = cl("Energy", 0.0);
   std::string molec_file=cl("mol_file", "invalid_file"); // this file contains all informations on the molecule
   std::string angular_creator=cl("angular", "invalid"); 
   Real r=cl("radius", 20.);
   int NrBall=cl("points", 50);
   Real VolConst= cl("maxVol", 1./(32.*sqrt(E*E*E)) );
   Real L=cl("bending", 2.);
   int N=cl("circles", 5);
   int maxiter=cl("maxiter", 700);
   bool cap = cl("cap", false);
   std::vector<Node> geometry=getGeometry(molec_file);
   // Use the internal mesh generator to create a uniform
   // 2D grid on a square.
   // be aware: it is pot file, not pot whale!
   std::string pot_file=cl("mesh_file", "none");
   assert(pot_file!="none");

   // the function below creates a mesh using the molecular structure.
   tetrahedralise_sphere(mesh, geometry, angular_creator, r, NrBall, VolConst, L, N);
   
   int order=cl("order", 1);
   FEType fe_type(FIRST);
   if (order==2){
      //convert element to second-order mesh.
      // In case of tetrahedra: from Tet4 to Tet10
      mesh.all_second_order();
      // Create an FEType describing the approximation
      // characteristics of the InfFE object. Note that
      // the constructor automatically defaults to some
      // sensible values.  But use FIRST order
      // approximation.
      FEType fe_type(SECOND);
      //Order radial_order=THIRD; // default value.
      //FEType fe_type(FIRST, LAGRANGE, radial_order);
   }

   // in case of infinite elements, they are added now. This is done in the following by an automatized interface
   // that finds the center of finite elemnts and so on.
   //ExodusII_IO (mesh).write("Molec_Mesh2.e");
   if (infel){
      // determine geometric center of molecule:
      InfElemBuilder::InfElemOriginValue com_x;
      InfElemBuilder::InfElemOriginValue com_y;
      InfElemBuilder::InfElemOriginValue com_z;
      com_x.first=true;
      com_y.first=true;
      com_z.first=true;
      for (unsigned int atom=0;atom<geometry.size(); atom++){
         com_x.second+=geometry[atom](0);
         com_y.second+=geometry[atom](1);
         com_z.second+=geometry[atom](2);
      }
      com_x.second/=geometry.size();
      com_y.second/=geometry.size();
      com_z.second/=geometry.size();
      InfElemBuilder builder(mesh);
      builder.build_inf_elem(com_x, com_y, com_z, 
                             false,  false,  false,
                             true, libmesh_nullptr);
   
      // Reassign subdomain_id() of all infinite elements.
      // Otherwise, the exodus-api will fail on the mesh.
      MeshBase::element_iterator       elem_it  = mesh.elements_begin();
      const MeshBase::element_iterator elem_end = mesh.elements_end();
      for (; elem_it != elem_end; ++elem_it){
         Elem* elem = *elem_it;
         if(elem->infinite()){
            elem->subdomain_id() = 1;
         }
      }

      // find the neighbours; for correct linking the two areas
      mesh.find_neighbors();
   }
   // Print information about the mesh to the screen.
   mesh.print_info();

   // Create an equation systems object.
   EquationSystems equation_systems (mesh);
   
   // Create a EigenSystem named "Eigensystem" and (for convenience)
   // use a reference to the system we create.
   CondensedEigenSystem & eigen_system = equation_systems.add_system<CondensedEigenSystem> ("EigenSE");
   LinearImplicitSystem & ESP = equation_systems.add_system<LinearImplicitSystem> ("ESP");
   LinearImplicitSystem & DO = equation_systems.add_system<LinearImplicitSystem> ("DO");

   equation_systems.parameters.set<std::string >("origin_mesh")=cl("mesh_geom", "sphere");
   equation_systems.parameters.set<bool >("cap")=cap;

   equation_systems.parameters.set<std::string>("potential")=pot_file;
   equation_systems.parameters.set<std::string>("DO_file")=molec_file;
   // Declare the system variables.
   // Adds the variables to the different equation systems.
   eigen_system.add_variable("phi", fe_type);
   ESP.add_variable("esp", fe_type);
   DO.add_variable("do", fe_type);
 
   // Give the system a pointer to the matrix assembly
   // function defined below.
   eigen_system.attach_assemble_function (assemble_InfSE);
   ESP.attach_assemble_function (assemble_ESP);
   DO.attach_assemble_function (assemble_DO);

   eigen_system.set_eigenproblem_type(GHEP);
   //eigen_system.set_eigenproblem_type(GNHEP);
   
   // Set necessary parametrs used in EigenSystem::solve(),
   // i.e. the number of requested eigenpairs \p nev and the number
   // of basis vectors \p ncv used in the solution algorithm. Note that
   // ncv >= nev must hold and ncv >= 2*nev is recommended.
   equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
   equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3+4;

   bool refinement=cl("refine", false);
   
   eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); // this is default
   //eigen_system.eigen_solver->set_eigensolver_type(ARNOLDI);
   //eigen_system.eigen_solver->set_eigensolver_type(LANCZOS);

   const std::string spect = cl("spect","sm");
   if (spect=="sm"){
      eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);
   }
   else if (spect=="lm"){
      eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
   }
   else if (spect=="sr"){
      eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_REAL);
   }
   else if (spect=="lr"){
      eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_REAL);
   }
   else if (spect=="si"){
      eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_IMAGINARY);
   }
   else if (spect=="li"){
      eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_IMAGINARY);
   }
   else{
      libmesh_error_msg("\nUnknown position in spectrum given. \n Only sm, lm, sr, lr, si, li are valid .\n");
   }
   // an alternative would be here: TARGET_MAGNITUDE --> get values closest to ...
   
   // Set the solver tolerance and the maximum number of iterations.
   equation_systems.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
   equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = maxiter;
   equation_systems.parameters.set<Real> ("radius") = r;
   Real power=cl("power",12.);
   equation_systems.parameters.set<Real> ("power") = power;
   Real gamma=cl("gamma",0.0);
   equation_systems.parameters.set<Real> ("gamma") = gamma;
   equation_systems.parameters.set<std::vector<Node>> ("mol_geom") = geometry;
   equation_systems.parameters.set<bool> ("cap") = cap;

   // Prints information about the system to the screen.
   //equation_systems.print_info();
         
   // Initialize the data structures for the equation system.
   eigen_system.init();
   ESP.init();
   DO.init();

    // add boundary conditions if not infinite elements used. In the latter case ...
   if (not infel){
      std::set<unsigned int> dirichlet_dof_ids;
      get_dirichlet_dofs(equation_systems, "EigenSE" ,dirichlet_dof_ids);
      eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
   }

   // Solve the system "Eigensystem".
   
   // In Do.solve(): set energy-offset and dyson norm
   DO.solve(); 
   ESP.solve();
   // set the photone energy:
   equation_systems.parameters.set<Real>("energy")=E-equation_systems.parameters.get<Real>("E_do");
   out<<"E_ph: "<<E<<"  ";
   out<<"E_do: "<<equation_systems.parameters.get<Real>("E_do")<<std::endl;
   out<<"E_kin:"<<equation_systems.parameters.get<Real>("energy")<<std::endl;
   // set the ESP as initial guess for solution vector.
   * eigen_system.solution = * ESP.solution; 

   //now, do refinement loop, if refinement is allowd:
   if (refinement){
      MeshRefinement mesh_refinement(mesh);
      //refine and coarsen elements with errors in 30-ths percentile:
      // in general: more coarsening than refinement!
      mesh_refinement.refine_fraction()=0.75;
      mesh_refinement.coarsen_fraction()=0.3;
      // maximum number of refinements for single element/:
      mesh_refinement.max_h_level()=7;
      unsigned int r_max=7;
      for(unsigned int r=0; r<=r_max; r++){
         eigen_system.solve();
         if (r<r_max){
            ErrorVector error;
            KellyErrorEstimator error_estimator;
            error_estimator.estimate_error(eigen_system, error);
            out<<"error estimate \n l2 norm="
               <<error.l2_norm()
               <<"\n maximum norm = "
               <<error.maximum()
               <<std::endl;

            mesh_refinement.flag_elements_by_error_fraction(error);
            // do the actual work:
            mesh_refinement.refine_and_coarsen_elements();

            equation_systems.reinit();
         }
      }
      // reinitialise and estimate the esp-system
      //ESP.reinit();
      ESP.solve();
      //DO.reinit();
      DO.solve();
   }
   else{
      // else: simply solve the system
      eigen_system.solve();

      // the other two are solved already!
      //ESP.solve();
      //DO.solve();
   }

   // Get the number of converged eigen pairs.
   unsigned int nconv = eigen_system.get_n_converged();

   std::cout << "Number of converged eigenpairs: " << nconv << "\n" << std::endl;
   
   std::ostringstream eigenvector_output_name;
   eigenvector_output_name<< "esp.cube";
   cube_io(equation_systems, geometry, eigenvector_output_name.str(), "ESP");
   eigenvector_output_name.str(std::string());
   eigenvector_output_name<< "do.cube";
   cube_io(equation_systems, geometry, eigenvector_output_name.str(), "DO");

   // Write the eigen vector to file.
   for(unsigned int i=0; i<nconv; i++){
      std::pair<Real,Real> eigpair = eigen_system.get_eigenpair(i);
      std::cout<<"kinetic energy: "<<i<<" = "<<eigpair.first+equation_systems.parameters.get<Real>("energy")<<std::endl;
      #ifdef LIBMESH_HAVE_EXODUS_API
         eigenvector_output_name.str(std::string());
         if (infel)
            eigenvector_output_name<< i <<"-"<<cl("pot","unknwn")<<"_inf.e" ;
         else
            eigenvector_output_name<< i <<"-"<<cl("pot","unknwn")<<".e" ;
         ExodusII_IO (mesh).write_equation_systems(eigenvector_output_name.str(), equation_systems);
      #endif // #ifdef LIBMESH_HAVE_EXODUS_API
      //eigenvector_output_name<< i <<"_err.e";
      //ErrorVector::plot_error(eigenvector_output_name.str(), equation_systems.get_mesh() );
      eigenvector_output_name.str(std::string());
      eigenvector_output_name<< "phi-"<<i <<".cube" ;
      cube_io(equation_systems, geometry, eigenvector_output_name.str(), "EigenSE");
   }
   if (nconv==0){
      // that one can look at the mesh and some properties...
      #ifdef LIBMESH_HAVE_EXODUS_API
         eigenvector_output_name.str(std::string());
         if (infel)
            eigenvector_output_name<<"U"<<"-"<<cl("pot","unknwn")<<"_inf.e" ;
         else
            eigenvector_output_name<<"U"<<"-"<<cl("pot","unknwn")<<".e" ;
         ExodusII_IO (mesh).write_equation_systems(eigenvector_output_name.str(), equation_systems);
      #endif // #ifdef LIBMESH_HAVE_EXODUS_API
      Real normDO = 0;
      if (!infel) 
         normDO= DO.calculate_norm(*DO.solution, 0, L2);
      out<<"norm of DO:   "<< normDO <<"  ";
      out<< sqrt(calculate_overlap(equation_systems, "DO", 0, "DO", 0, OVERLAP))<<"  ";
      out<< sqrt(overlap_DO(equation_systems, "DO", 0, OVERLAP))<<std::endl;
   }
   if(nconv>0){
      eigen_system.get_eigenpair(0);
      Real intensity=normalise(equation_systems, infel);
      out<<"intensity:  "<<intensity<<std::endl;
      for(unsigned int i=1; i<nconv; i++){
         out<<" for the solution nr "<<i<<":"<<std::endl;
         eigen_system.get_eigenpair(i);
         intensity=normalise(equation_systems, infel);
         out<<"intensity:  "<<intensity<<std::endl;
      }
   }
   Number normDO=equation_systems.parameters.get<Real>("DOnorm");
   out<<"Norm DO:    "<<normDO<<std::endl;

   // All done.
   return 0;
}
#endif // LIBMESH_HAVE_SLEPC
