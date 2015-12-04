// libMesh include files.
#include "libmesh/libmesh.h"
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

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the eigen system. Here, we will simply assemble a mass matrix.
void assemble_EigenSE(EquationSystems& es, const std::string& system_name);
void get_dirichlet_dofs(EquationSystems& es, const std::string& system_name, std::set<unsigned int>& global_dirichlet_dofs_set);

int main (int argc, char** argv){
   // Initialize libMesh and the dependent libraries.
   LibMeshInit init (argc, argv);

   // Skip SLEPc examples on a non-SLEPc libMesh build
   #ifndef LIBMESH_HAVE_SLEPC
     libmesh_example_requires(false, "--enable-slepc");
     }
   #else
  
   #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
    // SLEPc currently gives us a nasty crash with Real==float
    libmesh_example_requires(false, "--disable-singleprecision");
   #endif

   // Check for proper usage.
   if (argc < 3)
      libmesh_error_msg("\nUsage: " << argv[0] << " -n <number of eigen values>");
   // Tell the user what we are doing.
   else {
      std::cout << "Running " << argv[0];

      for (int i=1; i<argc; i++)
          std::cout << " " << argv[i];
      std::cout << std::endl << std::endl;
   }

   // Get the number of eigen values to be computed from argv[2]
   const unsigned int nev = std::atoi(argv[2]);

   // Skip this 2D example if libMesh was compiled as 1D-only.
   libmesh_example_requires(3 <= LIBMESH_DIM, "2D support");
   
   int dim = 3;
   // Create a mesh, with dimension to be overridden later, on the
   // default MPI communicator.
   Mesh mesh(init.comm(), dim);

   // Use the internal mesh generator to create a uniform
   // 2D grid on a square.
   //MeshTools::Generation::build_square (mesh, 40, 40, 0., 1., 0, 1.1, QUAD4);
   MeshTools::Generation::build_cube (mesh, 20, 20, 20, 0., 1., 0, 1.1, 0, 1.1, PRISM15);
   //MeshTools::Generation::build_sphere(mesh, 1., 10, QUAD4, 20, false);

   // Print information about the mesh to the screen.
   mesh.print_info();

   // Create an equation systems object.
   EquationSystems equation_systems (mesh);
   
   // Create a EigenSystem named "Eigensystem" and (for convenience)
   // use a reference to the system we create.
   CondensedEigenSystem & eigen_system = equation_systems.add_system<CondensedEigenSystem> ("EigenSE");

   // Declare the system variables.
   // Adds the variable "p" to "Eigensystem".   "p"
   // will be approximated using second-order approximation.
   eigen_system.add_variable("phi", SECOND);
   
   // Give the system a pointer to the matrix assembly
   // function defined below.
   eigen_system.attach_assemble_function (assemble_EigenSE);
   
   // Set necessary parametrs used in EigenSystem::solve(),
   // i.e. the number of requested eigenpairs \p nev and the number
   // of basis vectors \p ncv used in the solution algorithm. Note that
   // ncv >= nev must hold and ncv >= 2*nev is recommended.
   equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
   equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3;
   
   //eigen_system.set_eigenproblem_type(GNHEP);
   eigen_system.set_eigenproblem_type(GHEP);
   eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); // this is default
   eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);
   
   // Set the solver tolerance and the maximum number of iterations.
   equation_systems.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
   equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
   
   // Initialize the data structures for the equation system.
   equation_systems.init();

   // Prints information about the system to the screen.
   equation_systems.print_info();
    // add boundary conditions
    std::set<unsigned int> dirichlet_dof_ids;
    get_dirichlet_dofs(equation_systems, "EigenSE" ,dirichlet_dof_ids);
    eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
   // Solve the system "Eigensystem".
   eigen_system.solve();

   // Get the number of converged eigen pairs.
   unsigned int nconv = eigen_system.get_n_converged();

   std::cout << "Number of converged eigenpairs: " << nconv << "\n" << std::endl;

   #ifdef LIBMESH_HAVE_EXODUS_API
       // Write the eigen vector to file.
       for(unsigned int i=0; i<nconv; i++){
          eigen_system.get_eigenpair(i);
          std::ostringstream eigenvector_output_name;
          eigenvector_output_name<< i <<"_ep.e";
          ExodusII_IO (mesh).write_equation_systems ( eigenvector_output_name.str(), equation_systems);
       }
   #endif // #ifdef LIBMESH_HAVE_EXODUS_API

   // All done.
   return 0;
}
#endif // LIBMESH_HAVE_SLEPC

void assemble_EigenSE(EquationSystems& es, const std::string& system_name){
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "EigenSE");
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();
   std::cout<<dim<<std::endl;

  // Get a reference to our system.
  EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  // A reference to the system matrix
  SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
  SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  libMesh::Number E=0.0;
  libMesh::Number V=0.00;
  libMesh::Number co0_5= 0.5;

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap& dof_map = eigen_system.get_dof_map();

  // The element mass matrix.
  DenseMatrix<Number> Se;
  DenseMatrix<Number> H;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el){
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrices and rhs before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Se.resize (dof_indices.size(), dof_indices.size());
      H.resize (dof_indices.size(), dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //
      // We will build the element matrix.  This involves
      // a double loop to integrate the test funcions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<qrule.n_points(); qp++){
        for (unsigned int i=0; i<phi.size(); i++){
          for (unsigned int j=0; j<phi.size(); j++){
            Se(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
            H(i,j)  += JxW[qp]*(co0_5*(dphi[i][qp]*dphi[j][qp]) +(V- E)*phi[i][qp]*phi[j][qp]);
          }
        }
      }

      // On an unrefined mesh, constrain_element_matrix does
      // nothing.  If this assembly function is ever repurposed to
      // run on a refined mesh, getting the hanging node constraints
      // right will be important.  Note that, even with
      // asymmetric_constraint_rows = false, the constrained dof
      // diagonals still exist in the matrix, with diagonal entries
      // that are there to ensure non-singular matrices for linear
      // solves but which would generate positive non-physical
      // eigenvalues for eigensolves.
      dof_map.constrain_element_matrix(Se, dof_indices, false);
      dof_map.constrain_element_matrix(H, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrix.
      matrix_A.add_matrix (H, dof_indices);
      matrix_B.add_matrix (Se, dof_indices);

  } // end of element loop

  // Avoid compiler warnings
  libmesh_ignore(es);

  /**
   * All done!
   */
  return;
}

void get_dirichlet_dofs(EquationSystems& es, const std::string& system_name, std::set<unsigned int>& dirichlet_dof_ids){
   dirichlet_dof_ids.clear();
   // It is a good idea to make sure we are assembling
   // the proper system.
   libmesh_assert_equal_to (system_name, "EigenSE");
   
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
   std::cout<<dim<<std::endl;
   
   // Get a reference to our system.
   EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
   
   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = eigen_system.get_dof_map().variable_type(0);
   const DofMap& dof_map = eigen_system.get_dof_map();
   
   // This vector will hold the degree of freedom indices for
   // the element.  These define where in the global system
   // the element degrees of freedom get mapped.
   std::vector<dof_id_type> dof_indices;
   
   // Now we will loop over all the elements in the mesh that
   // live on the local processor. We will compute the element
   // matrix and right-hand-side contribution.  In case users
   // later modify this program to include refinement, we will
   // be safe and will only consider the active elements;
   // hence we use a variant of the \p active_elem_iterator.
   MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   
   for ( ; el != end_el; ++el){
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      {
        // All boundary dofs are Dirichlet dofs in this case
        for (unsigned int s=0; s<elem->n_sides(); s++){
            if (elem->neighbor(s) == NULL){
               std::vector<unsigned int> side_dofs;
               FEInterface::dofs_on_side(elem, dim, fe_type, s, side_dofs);

               for(unsigned int ii=0; ii<side_dofs.size(); ii++){
                  dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
               }
            }
         }
      }
   } // end of element loop
   
   /**
   * All done!
   */
   return;
}
