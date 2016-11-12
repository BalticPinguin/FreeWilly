# include "FreeWilly.h"
# include "SlepcConfig.h"

// for finding element for point
#include "libmesh/point_locator_tree.h"
// for dirichlet boundary conditions
#include "libmesh/fe_interface.h" 
#include "libmesh/fe_compute_data.h"
#include "libmesh/error_vector.h" 

// Bring in everything from the libMesh namespace
using namespace libMesh;

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
   
   // Check for proper usage and print the input file to output.
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
   

   //Read in the options provided in the input-file and store them in
   // local variables.
   // \p molec_file contains all informations on the molecule
   std::string molec_file=cl("mol_file", "invalid_file"); 
   std::string angular_creator=cl("angular", "invalid"); 
   Real r=cl("radius", 20.);
   std::string scheme=cl("scheme", "tm");
   Real p=cl("p", 1.0);
   Real E = cl("Energy",0.0);
   Real VolConst= cl("maxVol", 1./(2*E*E*E));
   Real L=cl("bending", 2.);
   int N=cl("circles", 5);
   int maxiter=cl("maxiter", 700);
   bool cap = cl("cap", false);
   bool refinement=cl("refine", false);
   bool quadrature_print = cl("print_quadrature", false);
   bool pictorious = cl("pictorious", false);
   int spherical_analysis= cl("spherical_analysis", -1);
   Real r_0=cl("r_0",12.);
   Real gamma=cl("gamma",0.0);

   // it is pot file, not pot whale!
   std::string pot_file=cl("mesh_file", "none");
   assert(pot_file!="none");

   DOrbit dyson(molec_file);
   Real Energy= E-dyson.get_energy();

   // make sure that the distance between two spheres is 
   // at least ~ 1/(4*lambda)
   if (N<=(int)(sqrt(Energy)*r/17.8))
      N=(int)(r*sqrt(Energy)/17.8);
 

   if (scheme=="tm" || scheme=="tm_300" ||
       scheme=="const_tm" || scheme=="sqrt_tm"){
      // assure that L is not larger than ~lambda/3 (i.e. 2*pi/(3*k)).
      if (L<sqrt(2./Energy))
         L=sqrt(2./Energy);

      // this is necessary to be in numerically stable regime.
      if( N<r/L)
         N=r/L;
   }
  
   // the function below creates a mesh using the molecular structure.
   tetrahedralise_sphere(mesh, dyson.geometry, angular_creator, r, scheme, p, VolConst, L, N);
   
   int order=cl("order", 1);
 
   // define the fe_type including infinite eleemnt parameters:
   //FEType fe_type(FIRST, LAGRANGE, FIRST, JACOBI_20_00, CARTESIAN);
   FEType fe_type(FIRST, LAGRANGE, FIFTH, JACOBI_20_00, CARTESIAN);
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
   if (infel){
      // determine geometric center of molecule:
      InfElemBuilder::InfElemOriginValue com_x;
      InfElemBuilder::InfElemOriginValue com_y;
      InfElemBuilder::InfElemOriginValue com_z;
      com_x.first=true;
      com_y.first=true;
      com_z.first=true;
      for (unsigned int atom=0;atom<dyson.geometry.size(); atom++){
         com_x.second+=dyson.geometry[atom](0);
         com_y.second+=dyson.geometry[atom](1);
         com_z.second+=dyson.geometry[atom](2);
      }
      com_x.second/=dyson.geometry.size();
      com_y.second/=dyson.geometry.size();
      com_z.second/=dyson.geometry.size();
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
   
   // Create the systems needed here: 
   // \p eigen_system solves the 1-p SE for the free electron
   CondensedEigenSystem & eigen_system = equation_systems.add_system<CondensedEigenSystem> ("EigenSE");

   // set the parameters of the calculation now as (globally available) paramers:
   equation_systems.parameters.set<std::string >("origin_mesh")=cl("mesh_geom", "sphere");
   equation_systems.parameters.set<bool >("cap")=cap;
   equation_systems.parameters.set<std::string>("potential")=pot_file;
   equation_systems.parameters.set<std::string>("DO_file")=molec_file;
   equation_systems.parameters.set<bool>("quadrat_print") = quadrature_print;

   // Declare the system variables.
   // Adds the variables to the different equation systems.
   eigen_system.add_variable("phi", fe_type);
 
   eigen_system.set_eigenproblem_type(GHEP);
   //eigen_system.set_eigenproblem_type(GNHEP);
   
   // Set necessary parametrs used in EigenSystem::solve(),
   // i.e. the number of requested eigenpairs \p nev and the number
   // of basis vectors \p ncv used in the solution algorithm. Note that
   // ncv >= nev must hold and ncv >= 2*nev is recommended.
   equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
   equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3+4;
   
   // chose among the solver options.  
   eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); // this is default
   //eigen_system.eigen_solver->set_eigensolver_type(LAPACK);  // this seems to be quite good.
   //eigen_system.eigen_solver->set_eigensolver_type(ARNOLDI);
   //eigen_system.eigen_solver->set_eigensolver_type(LANCZOS);
   
   // Set the solver tolerance and the maximum number of iterations.
   equation_systems.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
   equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = maxiter;
   equation_systems.parameters.set<Real> ("radius") = r;
   equation_systems.parameters.set<Real> ("r_0") = r_0;
   equation_systems.parameters.set<Real> ("gamma") = gamma;
   equation_systems.parameters.set<std::vector<Node>> ("mol_geom") = dyson.geometry;
   equation_systems.parameters.set<bool> ("cap") = cap;

   equation_systems.parameters.set<Real>("energy")=Energy;
   equation_systems.parameters.set<Number>("momentum")=sqrt((Energy)*2.);
   equation_systems.parameters.set<Real>("E_do")=dyson.get_energy();

   // Give the system a pointer to the matrix assembly
   // function defined below.
   eigen_system.attach_assemble_function (assemble_InfSE);
   if (pictorious){
      // \p ESP solves a simple matrix equation for the electrostatic potential;
      //   this system is for illustrative purposes only.
      LinearImplicitSystem & ESP = equation_systems.add_system<LinearImplicitSystem> ("ESP");
      // \p DO solves a simple matrix equation for the dyson orbital, also for illustrative purposes
      //   but here also some important options are set: the energy of the free electron and othels.
      LinearImplicitSystem & DO = equation_systems.add_system<LinearImplicitSystem> ("DO");

      ESP.add_variable("esp", fe_type);
      DO.add_variable("do", fe_type);

      ESP.attach_assemble_function (assemble_ESP);
      DO.attach_assemble_function (assemble_DO);

      ESP.init();
      DO.init();
   
      // In Do.solve(): set energy-offset and dyson norm
      DO.solve(); 
      ESP.solve();

      if (infel) 
         // set the ESP as initial guess for solution vector.
         // does not work for finite element due to different boundary conditions.
         eigen_system.eigen_solver->set_initial_space(*ESP.solution);
   }

   // Prints information about the system to the screen.
   //equation_systems.print_info();
         
   // Initialize the data structures for the equation system.
   eigen_system.init();
   //Real energy = E -equation_systems.parameters.get<Real>("E_do");
   libmesh_assert(Energy ==E -equation_systems.parameters.get<Real>("E_do"));

   eigen_system.eigen_solver->set_position_of_spectrum(
                          equation_systems.parameters.get<Real>("energy"));
 
   // add boundary conditions if not infinite elements used. 
   if (not infel){
      std::set<unsigned int> dirichlet_dof_ids;
      get_dirichlet_dofs(equation_systems, "EigenSE" ,dirichlet_dof_ids);
      eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
   }

   // print the energy values: photon energy, dyson orbital energy and kinetic energy
   out<<"E_ph: "<<cl("Energy",0.0)<<"  ";
   out<<"E_do: "<<equation_systems.parameters.get<Real>("E_do")<<"  ";
   out<<"E_kin:"<<Energy<<std::endl;

   // set in addition set a spectral transformation to stabilise the numerical scheme.
   SlepcEigenSolver<Number>* solver = 
                 libmesh_cast_ptr<SlepcEigenSolver<Number>* >( &(*eigen_system.eigen_solver) );

   SlepcSolverConfiguration ConfigSolver( *solver);
  
   // set the spectral transformation:
   ConfigSolver.SetST(SINVERT);
   //ConfigSolver.SetST(CAYLEY);
   //ConfigSolver.SetST(SHIFT); // this is default
   solver ->set_solver_configuration(ConfigSolver);

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
     // if (pictorious){
     //    //ESP.reinit();
     //    ESP.solve();
     //    //DO.reinit();
     //    DO.solve();
     // }
   }
   else
      // else: simply solve the system
      eigen_system.solve();

   // Get the number of converged eigen pairs.
   unsigned int nconv = eigen_system.get_n_converged();

   std::cout << "Number of converged eigenpairs: " << nconv << "\n" << std::endl;
   
   std::ostringstream eigenvector_output_name;
   if (pictorious){
      eigenvector_output_name<< "esp.cube";
      cube_io(equation_systems, dyson.geometry, eigenvector_output_name.str(), "ESP");
      eigenvector_output_name.str(std::string());
      eigenvector_output_name<< "do.cube";
      cube_io(equation_systems, dyson.geometry, eigenvector_output_name.str(), "DO");
   }

   // Write the eigen vector to file.
   nconv=std::min(nconv, nev);
   std::pair<Real,Real> eigpair;
   for(unsigned int i=0; i<nconv; i++){
      eigpair = eigen_system.get_eigenpair(i);
      
      std::cout<<"kinetic energy: "<<i<<" = "<<eigpair.first<<std::endl;
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
  
      // frequency=k/2*pi.
      equation_systems.parameters.set<Real>("current frequency")=sqrt(eigpair.first/pi);
      eigenvector_output_name.str(std::string());
      //eigenvector_output_name<< "phi-"<<i <<".cube";
      eigenvector_output_name<< "phi-"<<i <<".line";
      //cube_io(equation_systems, dyson.geometry, eigenvector_output_name.str(), "EigenSE");
      line_out(equation_systems, eigenvector_output_name.str(), "EigenSE");
   }
   if (nconv==0){
      // that one can look at the mesh and some properties...
      #ifdef LIBMESH_HAVE_EXODUS_API
         equation_systems.parameters.set<Real>("current frequency")=sqrt(eigpair.first/pi);
         if (infel)
            eigenvector_output_name<<"U"<<"-"<<cl("pot","unknwn")<<"_inf.e" ;
         else
            eigenvector_output_name<<"U"<<"-"<<cl("pot","unknwn")<<".e" ;
         ExodusII_IO (mesh).write_equation_systems(eigenvector_output_name.str(), equation_systems);
      #endif // #ifdef LIBMESH_HAVE_EXODUS_API
      //out<<"norm of DO: ";
      //out<< sqrt(norm_DO(equation_systems))<<std::endl;
   }
   else{
      Real intensity;
      for(unsigned int i=0; i<nconv; i++){
         out<<" for the solution nr "<<i<<":"<<std::endl;
         eigpair = eigen_system.get_eigenpair(i);
         equation_systems.parameters.set<Real>("current frequency")=sqrt(eigpair.first/pi);
         intensity=normalise(equation_systems, infel);
         out<<"intensity:  "<<intensity<<std::endl;
      }
   }
   out<<"norm of DO: ";
   out<< sqrt(norm_DO(equation_systems, false))<<std::endl;
   out<<" exact: "<<dyson.get_norm()<<std::endl;

   // this will become an option for the input later.
   //bool spherical_analysis=true;
   if (spherical_analysis>=0){
      for(unsigned int i=0; i<nconv; i++){
         eigpair = eigen_system.get_eigenpair(i);
         equation_systems.parameters.set<Real>("current frequency")=sqrt(eigpair.first/pi);
         ProjectSphericals (equation_systems, 5, i);
      }
   }

   // All done.
   return 0;
}
#endif // LIBMESH_HAVE_SLEPC

