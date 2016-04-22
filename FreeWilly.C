#include <assert.h> // 
#include <stdlib.h>    /* atof */
#include <vector>
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
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
#include "libmesh/explicit_system.h"
// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
#include <complex.h> // the infinite element version requires complex numbers explicitly.

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the eigen system. Here, we will simply assemble a mass matrix.
std::vector<Point> getGeometry(GetPot cl);

//prototypes of functions in assemble.C (that are called from out of it)
void get_dirichlet_dofs(libMesh::EquationSystems& , const std::string& , std::set<unsigned int>&);
void assemble_EigenSE(libMesh::EquationSystems& , const std::string&);
void assemble_InfSE(libMesh::EquationSystems & es, const std::string & system_name);
// this one maybe will never work!?
void assemble_InfFullSE(libMesh::EquationSystems & es, const std::string & system_name);

//This in the tetrahedralisation of a sphere
//void tetrahedralise_sphere(libMesh::MeshBase& mesh, const libMesh::Parallel::Communicator& comm);
// void tetrahedralise_sphere(mesh, const libMesh::Parallel::Communicator& comm);
void tetrahedralise_sphere(libMesh::UnstructuredMesh& mesh, std::vector<Point> geometry);

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
      std::cout << "Running " << argv[0];

      for (int i=1; i<argc; i++)
          std::cout << " " << argv[i];
      std::cout << std::endl << std::endl;
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

   std::string mesh_geom = cl("mesh_geom", "sphere");
   // Use the internal mesh generator to create a uniform
   // 2D grid on a square.
   //MeshTools::Generation::build_square (mesh, 40, 40, 0., 1., 0, 1.1, QUAD4);
   //MeshTools::Generation::build_cube (mesh, 20, 20, 20, 0., 1., 0, 1.1, 0, 1.1, PRISM15);
   //MeshTools::Generation::build_sphere(mesh, 1., 10, HEX8, 20, false);
   //MeshTools::Generation::build_cube (mesh, 50, 50, 50, -20., 20., -20., 20., -20., 20., PRISM6);
   //MeshTools::Generation::build_cube (mesh, 1, 1, 1, -2., 2., -2., 2., -2., 2., HEX8);
   //MeshTools::Generation::build_cube (mesh, 20, 20, 20, -20., 20., -20., 20., -20., 20., PRISM6);
   //MeshTools::Generation::build_cube (mesh, 1, 1, 1, -2., 2., -2., 2., -2., 2., PRISM6);
   if (mesh_geom=="sphere"){
      MeshTools::Generation::build_sphere(mesh, 7., 3, HEX8, 2, true);
      out<<"sphere"<<std::endl;}
   else if (mesh_geom=="box"){
      MeshTools::Generation::build_cube (mesh, 3, 3, 3, -2., 2., -2., 2., -2., 2., PRISM6);
      out<<"box"<<std::endl;}
   else if (mesh_geom=="own"){
      std::vector<Point> geometry;
      geometry=getGeometry(cl);
      // the function below creates a mesh using the molecular structure.
      tetrahedralise_sphere(mesh, geometry);
      std::string pot_file=cl("mesh_file", "none");
      assert(pot_file!="none");
   }
   else
      libmesh_error_msg("\nUnknown geometry of the finite element region.\n");

   // Print information about the mesh to the screen.
   mesh.print_info();

   // in case of infinite elements, they are added now. This is done in the following by an automatized interface
   // that finds the center of finite elemnts and so on.
   ExodusII_IO (mesh).write("Molec_Mesh.e");
   if (infel){
      InfElemBuilder builder(mesh);
      builder.build_inf_elem(true);
   
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

      // print info on new mesh
      mesh.print_info();

      // find the neighbours; for correct linking the two areas
      mesh.find_neighbors();
   }

   // Create an equation systems object.
   EquationSystems equation_systems (mesh);
   
   // Create a EigenSystem named "Eigensystem" and (for convenience)
   // use a reference to the system we create.
   CondensedEigenSystem & eigen_system = equation_systems.add_system<CondensedEigenSystem> ("EigenSE");

   equation_systems.parameters.set<std::string >("origin_mesh")=cl("mesh_geom", "sphere");
   if (mesh_geom=="own"){
      // can I set a parameter to a MeshFunction object? seemingly not!
      std::string pot_file=cl("mesh_file", "none");
      equation_systems.parameters.set<std::string>("potential")=pot_file;
      //out<<pot_file<<std::endl;
   }
   else
      equation_systems.parameters.set<std::string>("potential")=cl("pot","coul");
   // Declare the system variables.
   // Adds the variable "p" to "Eigensystem".   "p"
   // will be approximated using second-order approximation.
   //eigen_system.add_variable("phi", SECOND);
   eigen_system.add_variable("phi", FIRST);
   
   // Give the system a pointer to the matrix assembly
   // function defined below.
   if( infel ){
     eigen_system.attach_assemble_function (assemble_InfSE);
   }
   else {
     eigen_system.attach_assemble_function (assemble_InfSE);
     //eigen_system.attach_assemble_function (assemble_EigenSE);
   }
   eigen_system.set_eigenproblem_type(GHEP);
   //eigen_system.set_eigenproblem_type(GNHEP);
   
   // Set necessary parametrs used in EigenSystem::solve(),
   // i.e. the number of requested eigenpairs \p nev and the number
   // of basis vectors \p ncv used in the solution algorithm. Note that
   // ncv >= nev must hold and ncv >= 2*nev is recommended.
   equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
   equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3+4;

   // set energy-offset
   equation_systems.parameters.set<Number>("offset")= cl("Energy", 0.0);
   
   eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); // this is default
   //eigen_system.eigen_solver->set_eigensolver_type(ARNOLDI); // this is default
   //eigen_system.eigen_solver->set_eigensolver_type(LANCZOS); // this is default
   // other options: Power, Lapack, subscape, Arnoldi, Lanczoc, Krylovschur 

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
   equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 10000;
   
   // Initialize the data structures for the equation system.
   equation_systems.init();

   // Prints information about the system to the screen.
   equation_systems.print_info();

    // add boundary conditions if not infinite elements used. In the latter case ...
   if (not infel){
      std::set<unsigned int> dirichlet_dof_ids;
      get_dirichlet_dofs(equation_systems, "EigenSE" ,dirichlet_dof_ids);
      eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
   }
   // Solve the system "Eigensystem".
   eigen_system.solve();

   // Get the number of converged eigen pairs.
   unsigned int nconv = eigen_system.get_n_converged();

   std::cout << "Number of converged eigenpairs: " << nconv << "\n" << std::endl;

   #ifdef LIBMESH_HAVE_EXODUS_API
       // Write the eigen vector to file.
       for(unsigned int i=0; i<nconv; i++){
          std::pair<Real,Real> eigpair = eigen_system.get_eigenpair(i);
          std::cout<<"energy of state "<<i<<" = "<<eigpair.first+equation_systems.parameters.set<Number>("offset")<<std::endl;
          std::ostringstream eigenvector_output_name;
          if (infel){
             eigenvector_output_name<< i <<"-"<<cl("pot","unknwn")<<"_inf.e" ;
          }
          else{
             //eigenvector_output_name<< i <<"-"<<cl("pot","unknwn")<<"_inf.e" ;
             eigenvector_output_name<< i <<"-"<<cl("pot","unknwn")<<".e" ;
             //eigenvector_output_name<< i <<"-"<<cl("pot","unknwn")<<".e" ;
          }
          ExodusII_IO (mesh).write_equation_systems ( eigenvector_output_name.str(), equation_systems);
          //eigenvector_output_name<< i <<"_err.e";
          //ErrorVector::plot_error(eigenvector_output_name.str(), equation_systems.get_mesh() );
       }
   #endif // #ifdef LIBMESH_HAVE_EXODUS_API

   // All done.
   return 0;
}
#endif // LIBMESH_HAVE_SLEPC
      
std::vector<Point> getGeometry(GetPot cl){
   // extract the given geometry:
   std::string x=cl("x", ".");
   std::string y=cl("y", ".");
   std::string z=cl("z", ".");
   // get number of values. Therefore, count spaces in the string.
   int length_x=0;
   int length_y=0;
   int length_z=0;
   int strlength=x.length();
   for(int i=0; i<strlength; i++){
      if (x[i] == ',')
            length_x++; 
      if (y[i] == ',')
            length_y++; 
      if (z[i] == ',')
            length_z++; 
   }
   assert (length_x== length_y);
   assert (length_x== length_z);
   length_x++;
   //put them into a vector:
   std::string x_i, y_i, z_i;
   int x_j=0, y_j=0, z_j=0;
   std::vector<double> x_v(length_x), y_v(length_x), z_v(length_x);
   for(int i=0; i<strlength; i++){
      if (x[i] == ','){
         x_v[x_j]=atof(x_i.c_str());
         x_i="";
         x_j++; }
      else
         x_i.append(& x[i]);
      if (y[i] == ','){
         y_v[y_j]=atof(y_i.c_str());
         y_i="";
         y_j++;}
      else
         y_i.append(& y[i]);
      if (z[i] == ','){
         z_v[z_j]=atof(z_i.c_str());
         z_i="";
         z_j++;}
      else
         z_i.append(& z[i]);
   }
   x_v[x_j]=atof(x_i.c_str());
   y_v[y_j]=atof(y_i.c_str());
   z_v[z_j]=atof(z_i.c_str());
   
   std::vector<Point> geometry(length_x);
   for(int i=0; i<length_x; i++){
      geometry[i]=Point(x_v[i], y_v[i], z_v[i]);
   }
   return geometry;
}
