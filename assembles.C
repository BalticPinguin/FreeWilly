// libMesh include files.
// #include "assembles.h"
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


// Bring in everything from the libMesh namespace
using namespace libMesh;

double Harm(Real x, Real y, Real z){
   const Real x0=0.0;
   const Real y0=0.0;
   const Real z0=0.0;
   const Real a=10;
   const Real b=10;
   const Real c=10;
   const Real V0=(-a-b-c)*10;
   return V0+a*(x-x0)*(x-x0)+b*(y-y0)*(y-y0)+c*(z-z0)*(z-z0);
}

double Coul(Real x, Real y, Real z){
   const Real x0=0.0;
   const Real y0=0.0;
   const Real z0=0.0;
   const Real Z=1.;
   return -Z/sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
}

double Morse(Real x, Real y, Real z){
   const Real x0=0.0;
   const Real y0=0.0;
   const Real z0=0.0;
   const Real D=5.3;
   const Real a=0.5;
   return -D*( 1-exp(-a*sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))) );
}

void assemble_EigenSE(EquationSystems & es, const std::string & system_name){
   // It is a good idea to make sure we are assembling
   // the proper system.
   libmesh_assert_equal_to (system_name, "EigenSE");
      
   #ifdef LIBMESH_HAVE_SLEPC

   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();

   // Get a reference to our system.
   //EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);
   CondensedEigenSystem & eigen_system = es.get_system<CondensedEigenSystem> (system_name);
      
   const std::string & Pot = es.parameters.get<std::string>("potential");
   Number E = es.parameters.get<Number>("offset");
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
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
      
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
   const std::vector<Point>& q_point = fe->get_xyz();
      
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
      
   Real pot=0;
         
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
        if (Pot == "harm"){
		pot=Harm(q_point[qp](0), q_point[qp](1), q_point[qp](2));
	}
        else if (Pot == "coul"){
		pot=Coul(q_point[qp](0), q_point[qp](1), q_point[qp](2));
	}
        else if (Pot == "morse"){
		pot=Morse(q_point[qp](0), q_point[qp](1), q_point[qp](2));
	}
	else{
		pot=0;
        }
        //const Real pot=V(q_point[qp](0), q_point[qp](1), q_point[qp](2));
        for (unsigned int i=0; i<phi.size(); i++){
          for (unsigned int j=0; j<phi.size(); j++){
            Se(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
            H(i,j)  += JxW[qp]*(co0_5*(dphi[i][qp]*dphi[j][qp]) +(pot- E)*phi[i][qp]*phi[j][qp]);
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
   #endif //ifdef LIBMESH_HAVE_SLEPC
   return;
}

void get_dirichlet_dofs(EquationSystems & es, const std::string & system_name, std::set<unsigned int>& dirichlet_dof_ids){
   dirichlet_dof_ids.clear();
   // It is a good idea to make sure we are assembling
   // the proper system.
   libmesh_assert_equal_to (system_name, "EigenSE");
   
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
   
   // Get a reference to our system.
   CondensedEigenSystem & eigen_system = es.get_system<CondensedEigenSystem> (system_name);
   
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
        for (unsigned int s=0; s< elem->n_sides(); s++){
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

void assemble_InfSE(EquationSystems & es, const std::string & system_name){
   // It is a good idea to make sure we are assembling
   // the proper system.
   libmesh_assert_equal_to (system_name, "EigenSE");
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
      
   // Get a reference to our system.
   CondensedEigenSystem & eigen_system = es.get_system<CondensedEigenSystem> (system_name);
      
   const std::string & Pot = es.parameters.get<std::string>("potential");
   Number E = es.parameters.get<Number>("offset");
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
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
   //AutoPtr<FEBase> inf_fe (FEBase::build_InfFE(dim,fe_type));
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
      
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   libMesh::Number co0_5= 0.5;
   libMesh::Number co2= 2.;
   //libMesh::Number k=omega; //divided by c which is 0 in atomic units.
   // -->ik = -i*k => for neg. energy: exp(-i*sqrt(2E)*mu(x))= exp(-sqrt(2|E|)*mu(x)) ==> expon. decay in function.
   libMesh::Number ik=sqrt(co2*E)*(std::complex<double>)_Complex_I; // -->try this for now...
   if (real(E)<0){ // E<0:
     ik=sqrt(co2*abs(E)); // this gives exponential decay .
   }
   libMesh::Number temp; // -->try this for now...
      
   // A reference to the \p DofMap object for this system.  The \p DofMap
   // object handles the index translation from node and element numbers
   // to degree of freedom numbers.
   const DofMap& dof_map = eigen_system.get_dof_map();
      
   // The element mass matrix and Hamiltonian
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
   MeshBase::const_element_iterator       el  = mesh.active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      
   Real pot=0;
      
   for ( ; el != end_el; ++el){
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // unifyging finite and infinite elements
      FEBase * cfe = libmesh_nullptr;

      if (elem->infinite()){
         // We have an infinite element.  Let \p cfe point
         // to our \p InfFE object.  This is handled through
         // an UniquePtr.  Through the \p UniquePtr::get() we "borrow"
         // the pointer, while the \p  UniquePtr \p inf_fe is
         // still in charge of memory management.
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }
     
      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real>& JxW = cfe->get_JxW();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >& phi = cfe->get_phi();
      const std::vector<std::vector<RealGradient> >& dphi = cfe->get_dphi();
      const std::vector<Point>& q_point = cfe->get_xyz();
      // get extra data needed for infinite elements
      const std::vector<RealGradient>& dphase = cfe->get_dphase();
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D
      const std::vector<RealGradient>& dweight = cfe->get_Sobolev_dweight();
   
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
   
      // Zero the element matrices and rhs before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Se.resize (dof_indices.size(), dof_indices.size());
      H.resize (dof_indices.size(), dof_indices.size());

      //out<<"index size: "<<dof_indices.size()<<std::endl;
   
      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         //if (elem->infinite()){
         //   std::cout<<q_point[qp](0)<<"  ";
         //   std::cout<<q_point[qp](1)<<"  ";
         //   std::cout<<q_point[qp](2)<<"  ";
         //   std::cout<<dphase[qp]<<"  "<<weight[qp]<<std::endl;
         //}
         if (Pot == "harm")
            pot=Harm(q_point[qp](0), q_point[qp](1), q_point[qp](2));
         else if (Pot == "coul")
            pot=Coul(q_point[qp](0), q_point[qp](1), q_point[qp](2));
         else if (Pot == "morse")
            pot=Morse(q_point[qp](0), q_point[qp](1), q_point[qp](2));
         else
            pot=0;
         // Now, get number of shape functions:
         unsigned int n_sf = cfe->n_shape_functions();
         // loop over it:
         for (unsigned int i=0; i<n_sf; i++){
            for (unsigned int j=0; j<n_sf; j++){
            // this is changed here due the Petrov-Galerkin scheme. and works with finite and infinite elements.
            Se(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
            temp= dweight[qp]*phi[i][qp]*(dphi[j][qp]-ik*dphase[qp]*phi[j][qp])+
                  weight[qp]*(dphi[j][qp]*dphi[i][qp]-ik*ik*dphase[qp]*dphase[qp]*phi[i][qp]*phi[j][qp]+
                  ik*dphase[qp]*(phi[i][qp]*dphi[j][qp]-phi[j][qp]*dphi[i][qp]));
            H(i,j) += JxW[qp]*( co0_5*temp + (pot- E)*weight[qp]*phi[i][qp]*phi[j][qp]);
         }  }  
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
   //matrix_A.close();
   //matrix_B.close();
   //matrix_A.print();
   //matrix_B.print();
         
   /**
   * All done!
   */
   return;
}

void assemble_InfFullSE(EquationSystems & es, const std::string & system_name){
   // It is a good idea to make sure we are assembling
   // the proper system.
   libmesh_assert_equal_to (system_name, "EigenSE");
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
      
   // Get a reference to our system.
   CondensedEigenSystem & eigen_system = es.get_system<CondensedEigenSystem> (system_name);
      
   const std::string & Pot = es.parameters.get<std::string>("potential");
   //Number E = es.parameters.get<Number>("offset");
   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = eigen_system.get_dof_map().variable_type(0);
      
   // A reference to the system matrix
   //SparseMatrix<Number>&  matrix_A = *eigen_system.matrix_A;
   // SparseMatrix<Number>&  matrix_B = *eigen_system.matrix_B;
   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
   AutoPtr<FEBase> inf_fe (FEBase::build_InfFE(dim,fe_type));
      
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   libMesh::Number co0_5= 0.5;
      
   //libMesh::Number k=omega; //divided by c which is 0 in atomic units.
   //  // -->ik = -i*k => for neg. energy: exp(-i*sqrt(2E)*mu(x))= exp(-sqrt(2|E|)*mu(x)) ==> expon. decay in function.
   //  if (real(E)<0){ // E<0:
   //     libMesh::Number ik=sqrt((std::complex<double>)2.*abs(E)); // this gives exponential decay .
   //  }
   //  else{
   //    libMesh::Number ik=sqrt((std::complex<double>)2.*E)*(std::complex<double>)_Complex_I; // -->try this for now...
   //  }
   libMesh::Number temp; // -->try this for now...
      
   // A reference to the \p DofMap object for this system.  The \p DofMap
   // object handles the index translation from node and element numbers
   // to degree of freedom numbers.
   const DofMap& dof_map = eigen_system.get_dof_map();
      
   // The three matrices of interest for quadratic eigensystem
   DenseMatrix<Number> Se;
   DenseMatrix<Number> Ce;
   DenseMatrix<Number> Me;
      
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
   MeshBase::const_element_iterator       el  = mesh.active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      
   Real pot=0;
      
   for ( ; el != end_el; ++el){
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

     // unifyging finite and infinite elements
     FEBase* cfe=NULL; 

     if (elem->infinite()){
        cfe = inf_fe.get();
     }
     else{
        cfe = fe.get();
     }
     
     // The element Jacobian * quadrature weight at each integration point.
     const std::vector<Real>& JxW = cfe->get_JxW();

     // The element shape functions evaluated at the quadrature points.
     const std::vector<std::vector<Real> >& phi = cfe->get_phi();
     const std::vector<std::vector<RealGradient> >& dphi = cfe->get_dphi();
     const std::vector<Point>& q_point = cfe->get_xyz();
     // get extra data needed for infinite elements
     const std::vector<RealGradient>& dphase = cfe->get_dphase();
     const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D
     const std::vector<RealGradient>& dweight = cfe->get_Sobolev_dweight();

     // Compute the element-specific data for the current
     // element.  This involves computing the location of the
     // quadrature points (q_point) and the shape functions
     // (phi, dphi) for the current element.
     cfe->reinit (elem);

     // Zero the element matrices and rhs before
     // summing them.  We use the resize member here because
     // the number of degrees of freedom might have changed from
     // the last element.  Note that this will be the case if the
     // element type is different (i.e. the last element was a
     // triangle, now we are on a quadrilateral).
     Se.resize (dof_indices.size(), dof_indices.size());
     Ce.resize (dof_indices.size(), dof_indices.size());
     Me.resize (dof_indices.size(), dof_indices.size());

     // Now loop over the quadrature points.  This handles
     // the numeric integration.
     //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
     unsigned int max_qp = cfe->n_quadrature_points();
     for (unsigned int qp=0; qp<max_qp; qp++){
        //const Real pot=V(q_point[qp](0), q_point[qp](1), q_point[qp](2));
        if (Pot == "harm")
            pot=Harm(q_point[qp](0), q_point[qp](1), q_point[qp](2));
        else if (Pot == "coul")
            pot=Coul(q_point[qp](0), q_point[qp](1), q_point[qp](2));
        else if (Pot == "morse")
            pot=Morse(q_point[qp](0), q_point[qp](1), q_point[qp](2));
	else
            pot=0;
        // Now, get number of shape functions:
        unsigned int n_sf = cfe->n_shape_functions();
        // loop over it:
        for (unsigned int i=0; i<n_sf; i++){
          for (unsigned int j=0; j<n_sf; j++){
            // this is changed here due the Petrov-Galerkin scheme...
            Se(i,j) +=-JxW[qp]*( (1-pot)*weight[qp]*phi[i][qp]*phi[j][qp]
                      +co0_5*dweight[qp]*phi[i][qp]*dphi[j][qp]
                      +co0_5*weight[qp]*dphi[j][qp]*dphi[i][qp]);
            Ce(i,j) +=JxW[qp]*(-co0_5*dphase[qp]*phi[j][qp]*phi[i][qp]*dweight[qp]
                      +co0_5*dphase[qp]*weight[qp]*phi[i][qp]*dphi[j][qp]
                      -co0_5*dphase[qp]*weight[qp]*phi[j][qp]*dphi[i][qp]);
            Me(i,j) +=-JxW[qp]*co0_5*weight[qp]*phi[i][qp]*phi[j][qp]*(dphase[qp]*dphase[qp]-1);
        }  }  
      }
      // for mesh-refinement:
      dof_map.constrain_element_matrix(Se, dof_indices, false);
      dof_map.constrain_element_matrix(Ce, dof_indices, false);
      dof_map.constrain_element_matrix(Me, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrix.
      //
      //(Se+ik*Ce+ik*ik*Me)c_vec=0 -> search for  ik, c_vec
   
   } // end of element loop
      
   // Avoid compiler warnings
   libmesh_ignore(es);
         
   /**
   * All done!
   */
   return;
}
