#include "FreeWilly.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

struct ESP{
   // these vectors store the points and potential (at the points) as given.
   std::vector<Point> node;
   std::vector<Number> potential;
   unsigned int size;
};

void Read(ESP& esp, std::string input_file){
   double noonecares;
   double x,y,z,v;
   std::ifstream esp_file;

   esp_file.open(input_file.c_str() );
   if (!esp_file.is_open()){
      //give an error message.
      libmesh_error_msg("Unable to open the file containing the electronstatic potential");
   }

   esp_file>>esp.size>>noonecares;
   // this is some error... maybe this way it works as well....
   esp.node.resize(esp.size+1);
   esp.potential.resize(esp.size+1);
   //for(std::string line; getline(esp_file, line); i++){
   for(unsigned int i=0; !esp_file.eof() ; i++){
     esp_file >> x>>y>>z>>v;
      esp.node[i]=Point(x,y,z);
      // the given potential seems to be the negative of the
      // normal potential; at least the numbers of He+ are all positive...
      esp.potential[i]=-v;
   }
   esp_file.close();
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
   EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);

   //const std::string & mesh_origin = es.parameters.get<std::string >("origin_mesh");
   const std::string & potfile = es.parameters.get<std::string>("potential");
   const std::string & formulation = es.parameters.get<std::string>("formulation");
   bool cap=es.parameters.get<bool >("cap");
   Real radius=es.parameters.get<Real>("radius");
   Real gamma =es.parameters.get<Real>("gamma");
   Real r_0=es.parameters.get<Real>("r_0");
   std::vector<Node> mol_geom=es.parameters.get<std::vector<Node>> ("mol_geom");
   bool quadrature = es.parameters.get<bool>("quadrat_print");
   Real power=es.parameters.get<Real> ("power");
   assert(power>0.);
   assert(power<1.);

   struct ESP esp;
   Read(esp, potfile);

   //InverseDistanceInterpolation<3> potential(mesh.comm(), 8, r_0);
   RBFInterpolation<3> potential(mesh.comm(), 9, r_0, mol_geom);
   //NeNeInterpolation<3> potential(mesh.comm(), 1, r_0, mol_geom);
   const std::vector<std::string> esp_data(1);
   potential.set_field_variables(esp_data);
   potential.add_field_data(esp_data, esp.node, esp.potential);

   potential.prepare_for_use();
     
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
   //QGauss qrule (dim, SIXTH);
   //QGauss qrule (dim, TWENTIETH);
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   //libMesh::Number k=omega; //divided by c which is 1 in atomic units.
   // -->ik = -i*k 
   Number ik=sqrt(-(Number)1.)*es.parameters.get<Number>("momentum");
   //out<<" i*k= "<<ik<<es.parameters.get<Number>("momentum")<<std::endl;
   //ik=sqrt(-es.parameters.get<Real>("energy")*(Number)2);
   //out<<" i*k= "<<ik<<sqrt(-(Number)1.)<<std::endl;

   // set parameters for infinite elements:
   es.parameters.set<Real>("speed")=1.;

   Number temp; 
      
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
      
   Number pot;
      
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
      
      std::vector<Number> potval;  //initialisation just for dummy reasons
   
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
   
      // Zero the element matrices before summing them.  
      // We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Se.resize (dof_indices.size(), dof_indices.size());
      H.resize (dof_indices.size(), dof_indices.size());
      potential.interpolate_field_data(esp_data, q_point, potval);

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      out.setf(std::ios::fixed, std::ios::floatfield);
      out<<std::setprecision(10);
      for (unsigned int qp=0; qp<max_qp; qp++){
	 bool away=true;
         if(quadrature){
	    for(int i=0; i<mol_geom.size(); i++){
                // 0.01 as threshold does work...
		if ((q_point[qp]-mol_geom[i]).norm()<0.0005){
		   away=false;
                   break;
                }
            }
            if(away){
            out<<"quadrature:";
            out<<std::setw(20)<<q_point[qp](0);
            out<<std::setw(20)<<q_point[qp](1);
            out<<std::setw(20)<<q_point[qp](2);
            out<<std::setw(20)<<potval[qp].real()<<std::endl;
            }
         }
         //out<<1/(q_point[qp].norm())<<std::endl;
         pot=potval[qp];
         if (cap){
            Real mindist=6000;
            for(unsigned int site=0; site<mol_geom.size(); site++){
               if ((q_point[qp]-mol_geom[site]).norm()<mindist)
                  mindist=(q_point[qp]-mol_geom[site]).norm();
            }
            if (mindist>=radius)
               pot=potval[qp]+Number(0,gamma*mindist*mindist);
         }

         // Now, get number of shape functions that are nonzero at this point::
         unsigned int n_sf = cfe->n_shape_functions();
         // loop over them:
         Number factor= 0.25*dweight[qp]*dweight[qp]/weight[qp]-ik*ik*dphase[qp]*dphase[qp]*weight[qp];

         for (unsigned int i=0; i<n_sf; i++){
            for (unsigned int j=0; j<n_sf; j++){
               if (formulation=="original"){
                  Se(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
                  temp= dweight[qp]*phi[i][qp]*(dphi[j][qp]-ik*dphase[qp]*phi[j][qp])+
                        weight[qp]*(dphi[j][qp]*dphi[i][qp]-ik*ik*dphase[qp]*dphase[qp]*phi[i][qp]*phi[j][qp]+
                        ik*dphase[qp]*(phi[i][qp]*dphi[j][qp]-dphi[i][qp]*phi[j][qp]));
                  H(i,j) += JxW[qp]*(0.5*temp + pot*weight[qp]*phi[i][qp]*phi[j][qp]);
               }
               else if (formulation=="squared"){
                  Se(i,j) += JxW[qp]*weight[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
                  temp= 2.*dweight[qp]*phi[i][qp]*(dphi[j][qp]-ik*dphase[qp]*phi[j][qp])+
                        weight[qp]*(dphi[j][qp]*dphi[i][qp]-ik*ik*dphase[qp]*dphase[qp]*phi[i][qp]*phi[j][qp]+
                        ik*dphase[qp]*(phi[i][qp]*dphi[j][qp]-dphi[i][qp]*phi[j][qp]));
                  H(i,j) += JxW[qp]*weight[qp]*(0.5*temp + pot*weight[qp]*phi[i][qp]*phi[j][qp]);
               }
               else if (formulation=="symmetric"){
                  Se(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
                  temp= phi[i][qp]*phi[j][qp]*factor+
                        0.5*dweight[qp]*(phi[i][qp]*dphi[j][qp]+dphi[i][qp]*phi[j][qp]) +
                        weight[qp]*(dphi[j][qp]*dphi[i][qp]
                                    +ik*dphase[qp]*(dphi[i][qp]*phi[j][qp]-phi[i][qp]*dphi[j][qp]));
                  H(i,j) += JxW[qp]*(0.5*temp + pot*weight[qp]*phi[i][qp]*phi[j][qp]);
               }
               else if (formulation=="root"){
                  Se(i,j) += JxW[qp]*sqrt(weight[qp])*phi[i][qp]*phi[j][qp];
                  temp= (dweight[qp]*dweight[qp]*phi[i][qp]*phi[j][qp]*0.25/weight[qp] +
                         dweight[qp]*(phi[i][qp]*dphi[j][qp]+dphi[i][qp]*phi[j][qp]))/(4*weight[qp])+
                        dphi[i][qp]*dphi[j][qp]+
                        ik*dphase[qp]*(dphi[i][qp]*phi[j][qp]-phi[i][qp]*dphi[j][qp]-
                                      ik*dphase[qp]*phi[i][qp]*phi[j][qp]);
                  H(i,j) += JxW[qp]*sqrt(weight[qp])*(0.5*temp + pot*phi[i][qp]*phi[j][qp]);
               }
               else if (formulation=="power"){
                  Se(i,j) += JxW[qp]*pow(weight[qp],2.*power)*phi[i][qp]*phi[j][qp];
                //  temp= power*power*dweight[qp]*dweight[qp]*phi[i][qp]*phi[j][qp]*pow(weight[qp],2*power-2)+
                //        power*dweight[qp]*pow(dweight[qp],2*power-1)*(phi[i][qp]*dphi[j][qp]+dphi[i][qp]*phi[j][qp])+
                //        pow(dweight[qp],2*power)*dphi[i][qp]*dphi[j][qp]+
                //        pow(dweight[qp],2*power)*ik*dphase[qp]*(dphi[i][qp]*phi[j][qp]-phi[i][qp]*dphi[j][qp])-
                //        pow(dweight[qp],2*power)*ik*dphase[qp]*ik*dphase[qp]*phi[i][qp]*phi[j][qp];
                //  H(i,j) += JxW[qp]*(0.5*temp + pow(dweight[qp],2*power)*pot*phi[i][qp]*phi[j][qp]);
                  temp= (power*dweight[qp]*phi[i][qp]*phi[j][qp]/weight[qp]+
                         phi[i][qp]*dphi[j][qp]+dphi[i][qp]*phi[j][qp])*dweight[qp]*power/weight[qp]+
                        dphi[i][qp]*dphi[j][qp]+
                        ik*dphase[qp]*(dphi[i][qp]*phi[j][qp]-phi[i][qp]*dphi[j][qp]-
                                       ik*dphase[qp]*phi[i][qp]*phi[j][qp]);
                  H(i,j) += JxW[qp]*pow(weight[qp],2.*power)*(0.5*temp + pot*phi[i][qp]*phi[j][qp]);
               }
               else{
                  std::cerr<<"Formulation not known.";
                  assert(false);
               }
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
      /* the copy is needed since dof_indices is changed
       * in constrain_element_matrix().
       */
      std::vector<dof_id_type> dof_indices2=dof_indices;
      dof_map.constrain_element_matrix(H, dof_indices2, false);
      dof_map.constrain_element_matrix(Se, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrix.
      matrix_A.add_matrix (H, dof_indices);
      matrix_B.add_matrix (Se, dof_indices2);

   } // end of element loop
   matrix_A.close();
   matrix_B.close();

   // Aout and Bout need to be initialised differently
   // std::ostream Aout;
   // Aout<<"sparsity_A";
   // std::ostream Bout;
   // Bout<<"sparsity_B";

   //matrix_A.print(out,true);
   //matrix_B.print(out,true);
      
   /**
   * All done!
   */
   return;
}

void assemble_ESP(EquationSystems & es, const std::string & system_name){
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
      
   // Get a reference to our system.
   LinearImplicitSystem & eigen_system = es.get_system<LinearImplicitSystem> (system_name);

   Real r_0=es.parameters.get<Real>("r_0");
   const std::string & potfile = es.parameters.get<std::string>("potential");
   std::vector<Node> mol_geom=es.parameters.get<std::vector<Node>> ("mol_geom");
   bool cap=es.parameters.get<bool >("cap");
   Real radius=es.parameters.get<Real>("radius")*0.9;
   Real gamma =es.parameters.get<Real>("gamma");
   
   struct ESP esp;
   Read(esp, potfile);

   //InverseDistanceInterpolation<3> potential(mesh.comm(), 8, r_0);
   RBFInterpolation<3> potential(mesh.comm(), 9, r_0, mol_geom);
   //NeNeInterpolation<3> potential(mesh.comm(), 1, r_0, mol_geom);
   const std::vector<std::string> esp_data(1);
   potential.set_field_variables(esp_data);
   potential.add_field_data(esp_data, esp.node, esp.potential);

   potential.prepare_for_use();
      
   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = eigen_system.get_dof_map().variable_type(0);
      
   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
   //QGauss qrule (dim, SIXTH);
   //QGauss qrule (dim, TWENTIETH);
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   // A reference to the \p DofMap object for this system.  The \p DofMap
   // object handles the index translation from node and element numbers
   // to degree of freedom numbers.
   const DofMap& dof_map = eigen_system.get_dof_map();
      
   // The element mass matrix and Hamiltonian
   DenseMatrix<Number> M;
   DenseVector<Number> f;
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
      
   Number pot=0;
      
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
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }
     
      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real>& JxW = cfe->get_JxW();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >& phi = cfe->get_phi();
      const std::vector<Point>& q_point = cfe->get_xyz();
      std::vector<Number> potval;
      // get extra data needed for infinite elements
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D
   
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
      M.resize (dof_indices.size(), dof_indices.size());
      f.resize (dof_indices.size());
      potential.interpolate_field_data(esp_data, q_point, potval);

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         //out<<1/(q_point[qp].norm())<<std::endl;
         pot=potval[qp];
         if (cap){
            Real mindist=6000;
            for(unsigned int site=0; site<mol_geom.size(); site++){
               if ((q_point[qp]-mol_geom[site]).norm()<mindist)
                  mindist=(q_point[qp]-mol_geom[site]).norm();
            }
            if (mindist>=radius)
               pot=potval[qp]+Number(0,gamma*mindist*mindist);
         }
         
         // Now, get number of shape functions:
         unsigned int n_sf = cfe->n_shape_functions();
         for (unsigned int i=0; i<n_sf; i++){
            f(i)+=JxW[qp]*weight[qp]*phi[i][qp]*pot;
            for (unsigned int j=0; j<n_sf; j++){
               M(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
            }
         }
      }
      //dof_map.heterogenously_constrain_element_matrix_and_vector (M, f, dof_indices);
      dof_map.constrain_element_matrix_and_vector (M, f, dof_indices);

      eigen_system.matrix->add_matrix (M, dof_indices);
      eigen_system.rhs->add_vector (f, dof_indices);

   } // end of element loop
         
   /**
   * All done!
   */
   return;
}

std::vector<Number> evalSphWave(int l_max, Point qp, Real k);

void assemble_DO(EquationSystems & es, const std::string & system_name){
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
   
   // Get a reference to our system.
   LinearImplicitSystem & eigen_system = es.get_system<LinearImplicitSystem> (system_name);

   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = eigen_system.get_dof_map().variable_type(0);

   DOrbit dyson(es.parameters.get<std::string>("DO_file"));
      
   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
   //QGauss qrule (dim, SIXTH);
   //QGauss qrule (dim, TWENTIETH);

   //Real k=es.parameters.get<Real>("current frequency")*2*3.1415926;
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   // A reference to the \p DofMap object for this system.  The \p DofMap
   // object handles the index translation from node and element numbers
   // to degree of freedom numbers.
   const DofMap& dof_map = eigen_system.get_dof_map();
      
   // The element mass matrix and Hamiltonian
   DenseMatrix<Number> M;
   DenseVector<Number> f;
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
      
   Number do_val;
      
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
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }
     
      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real>& JxW = cfe->get_JxW();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >& phi = cfe->get_phi();
      const std::vector<Point>& q_point = cfe->get_xyz();
      // get extra data needed for infinite elements
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D
   
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
      M.resize (dof_indices.size(), dof_indices.size());
      f.resize (dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         // get the value of DO at q_point[qp]
         do_val=dyson.evalDO(q_point[qp]);

         // 
         //do_val=evalSphWave(1, 0, q_point[qp], k);
         
         // Now, get number of shape functions:
         unsigned int n_sf = cfe->n_shape_functions();
         //out<<"init: "<<q_point[qp];
         //out<<"  "<<do_val<<std::endl;
         for (unsigned int i=0; i<n_sf; i++){
            f(i)+=JxW[qp]*weight[qp]*phi[i][qp]*do_val;
            for (unsigned int j=0; j<n_sf; j++){
               M(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
            }
         }
      }
      dof_map.constrain_element_matrix_and_vector (M, f, dof_indices);

      eigen_system.matrix->add_matrix (M, dof_indices);
      eigen_system.rhs->add_vector (f, dof_indices);

   } // end of element loop
         
   /**
   * All done!
   */
   return;
}

void assemble_Spherical(EquationSystems & es, const std::string & system_name){
   // Get a constant reference to the mesh object.
   const MeshBase& mesh = es.get_mesh();
   // The dimension that we are running.
   const unsigned int dim = mesh.mesh_dimension();
   
   // Get a reference to our system.
   LinearImplicitSystem & eigen_system = es.get_system<LinearImplicitSystem> (system_name);

   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = eigen_system.get_dof_map().variable_type(0);

   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));  // here, try AutoPtr instead...
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (dim, fe_type.default_quadrature_order());
   //QGauss qrule (dim, SIXTH);
   //QGauss qrule (dim, TWENTIETH);

   Real k=es.parameters.get<Real>("current frequency")*2*3.1415926;
   Real l=es.parameters.get<int>("L_guess");
   Real m=es.parameters.get<int>("M_guess");
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
      
   // A reference to the \p DofMap object for this system.  The \p DofMap
   // object handles the index translation from node and element numbers
   // to degree of freedom numbers.
   const DofMap& dof_map = eigen_system.get_dof_map();
      
   // The element mass matrix and Hamiltonian
   DenseMatrix<Number> M;
   DenseVector<Number> f;
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
      
   std::vector<Number> Value;
   int index=l*l+l+1+m; // this is the index of vector Value to evaluate.
      
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
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }
     
      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real>& JxW = cfe->get_JxW();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >& phi = cfe->get_phi();
      const std::vector<Point>& q_point = cfe->get_xyz();
      // get extra data needed for infinite elements
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D
   
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);
      M.resize (dof_indices.size(), dof_indices.size());
      f.resize (dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){
         // get the value of DO at q_point[qp]
         Value=evalSphWave(l, q_point[qp], k);
         
         // Now, get number of shape functions:
         unsigned int n_sf = cfe->n_shape_functions();
         for (unsigned int i=0; i<n_sf; i++){
            f(i)+=JxW[qp]*weight[qp]*phi[i][qp]*Value[index];
            for (unsigned int j=0; j<n_sf; j++){
               M(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
            }
         }
      }
      dof_map.constrain_element_matrix_and_vector (M, f, dof_indices);

      eigen_system.matrix->add_matrix (M, dof_indices);
      eigen_system.rhs->add_vector (f, dof_indices);

   } // end of element loop
         
   /**
   * All done!
   */
   return;
}
