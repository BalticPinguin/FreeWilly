#include <vector>
#include <math.h> // needed for sqrt function
#include <cstring> // for strlen.
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/linear_implicit_system.h"

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/point.h"
#include "libmesh/system.h"
#include "libmesh/elem.h"
#include "libmesh/fe_type.h"

// includes for calculate_norm, point_*
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_tools.h"

void getDyson(const char *filename, int namelength, std::vector<std::vector<double> >& do_j, std::vector<unsigned int>& l,std::vector<double>& alpha, double&  energy, double& normDO);
double evalDO(const std::vector<std::vector<double> >& do_j, const std::vector<unsigned int>& l, const std::vector<double>& alpha, const std::vector<libMesh::Node>& geometry, const libMesh::Point pt);
std::vector<libMesh::Node> getGeometry(std::string fname);

using namespace libMesh;

enum IntegralType: int{
   MU=0,
   OVERLAP
};

Number calculate_overlap(EquationSystems& eq_sys, const std::string sys1, int var1, const std::string sys2, int var2 , IntegralType int_type){
   //run at all processors at once:
   //parallel_object_only(); --> can not be used here.

   Number overlap=0;
   
   // this should be checked somehow as well:
   //libmesh_assert_not_equal_to(es1.comm(), es2.comm());
   //libmesh_assert_not_equal_to(es1.get_dof_map().variable_type(var1),
   //                            es2.get_dof_map().variable_type(var2));
   //CondensedEigenSystem& es1=equation_systems.get_system<CondensedEigenSystem> (sys1);
   //LinearImplicitSystem & es2 = equation_systems.get_system<LinearImplicitSystem> (sys2);
   System & es1 = eq_sys.get_system<System> (sys1);
   System & es2 = eq_sys.get_system<System> (sys2);

   // Localize the potentially parallel vectors
   UniquePtr<NumericVector<Number> > local_v1 = NumericVector<Number>::build(es1.comm());
   local_v1->init((*es1.solution).size(), true, SERIAL);
   (*es1.solution).localize (*local_v1, es1.get_dof_map().get_send_list());
   UniquePtr<NumericVector<Number> > local_v2 = NumericVector<Number>::build(es2.comm());
   local_v2->init((*es2.solution).size(), true, SERIAL);
   (*es2.solution).localize (*local_v2, es2.get_dof_map().get_send_list());

   const FEType & fe_type = es2.get_dof_map().variable_type(var2);
   // Allow space for dims 0-3, even if we don't use them all
   std::vector<FEBase *> fe_ptrs(4,libmesh_nullptr);
   std::vector<QBase *> q_rules(4,libmesh_nullptr);

   const std::set<unsigned char> & elem_dims = es1.get_mesh().elem_dimensions();
   // Prepare finite elements for each dimension present in the mesh
   for (std::set<unsigned char>::const_iterator d_it = elem_dims.begin(); d_it != elem_dims.end(); ++d_it){
      q_rules[*d_it] =fe_type.default_quadrature_rule (*d_it).release();

      // Construct finite element object

      fe_ptrs[*d_it] = FEBase::build(*d_it, fe_type).release();

      // Attach quadrature rule to FE object
      fe_ptrs[*d_it]->attach_quadrature_rule (q_rules[*d_it]);
   }

   std::vector<dof_id_type> dof_indices;
   FEBase * cfe = libmesh_nullptr;

   // Begin the loop over the elements
   MeshBase::const_element_iterator       el     = es1.get_mesh().active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = es1.get_mesh().active_local_elements_end();
   for ( ; el != end_el; ++el){
       const Elem * elem = *el;
       const unsigned int dim = elem->dim();

      //QGauss qrule (dim, FIFTH);
      QGauss qrule (dim, fe_type.default_quadrature_order());
      UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
      UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
      fe->attach_quadrature_rule (&qrule);
      inf_fe->attach_quadrature_rule (&qrule);
      
      if (elem->infinite())
         cfe = inf_fe.get();
      else
         cfe = fe.get();
      const std::vector<Real> &  JxW     = cfe->get_JxW();
      const std::vector<Point> & q_point = cfe->get_xyz();

      cfe->reinit(elem);

      es1.get_dof_map().dof_indices (elem, dof_indices, var1);

      //const unsigned int n_qp = qrule.n_points(); --> fails for infinite elements.
      unsigned int n_qp = cfe->n_quadrature_points();
      const unsigned int n_sf = cast_int<unsigned int>(dof_indices.size());

      // Begin the loop over the Quadrature points.
      for (unsigned int qp=0; qp<n_qp; qp++){
         Number u_h = 0.;
         Number v_h = 0.;
               
         Point mapped_qp = FEInterface::inverse_map(dim, fe_type, elem, q_point[qp], TOLERANCE, true); 
         FEComputeData fe_data(eq_sys, mapped_qp);
         FEInterface::compute_data(dim, fe_type, elem, fe_data);

         for (unsigned int i=0; i != n_sf; ++i){
            //Use FEComputeData because with infinite elements the value at q_point[i][qp]
            // is not just phi[i][qp].
            u_h += fe_data.shape[i] * (*local_v1)(dof_indices[i]);
            v_h += fe_data.shape[i] * (*local_v2)(dof_indices[i]);
            
         }
         //out<<"later: "<<q_point[qp];
         //out<<"  "<<u_h<<std::endl;

         //norm += JxW[qp] * TensorTools::norm_sq(u_h);
         if(int_type==MU)
            overlap += JxW[qp] * std::conj(u_h) * q_point[qp].norm() * v_h;
         else //if(int_type==OVERLAP)
            overlap += JxW[qp] * std::conj(u_h) * v_h;
      }
   }
   // Need to delete the FE and quadrature objects to prevent a memory leak
   for(unsigned int i=0; i<fe_ptrs.size(); i++){
       if(fe_ptrs[i])
         {
           delete fe_ptrs[i];
         }
     }
   for(unsigned int i=0; i<q_rules.size(); i++){
       if(q_rules[i])
         {
           delete q_rules[i];
         }
     }

   es1.comm().sum(overlap);
   //overlap = std::sqrt(overlap);

   return overlap;
}

Number overlap_DO(EquationSystems& eq_sys, const std::string sys1, int var1, IntegralType int_type){
   //run at all processors at once:
   //parallel_object_only(); --> can not be used here.

   Number overlap=0;
   
   // this should be checked somehow as well:
   //libmesh_assert_not_equal_to(es1.comm(), es2.comm());
   //libmesh_assert_not_equal_to(es1.get_dof_map().variable_type(var1),
   //                            es2.get_dof_map().variable_type(var2));
   //CondensedEigenSystem& es1=equation_systems.get_system<CondensedEigenSystem> (sys1);
   //LinearImplicitSystem & es2 = equation_systems.get_system<LinearImplicitSystem> (sys2);
   System & es1 = eq_sys.get_system<System> (sys1);

   // Localize the potentially parallel vectors
   UniquePtr<NumericVector<Number> > local_v1 = NumericVector<Number>::build(es1.comm());
   local_v1->init((*es1.solution).size(), true, SERIAL);
   (*es1.solution).localize (*local_v1, es1.get_dof_map().get_send_list());

   // get the dyson orbital:
   std::vector<std::vector<double> > do_j;
   std::vector<unsigned int> l;
   std::vector<double> alpha;
   std::vector<Node> geometry= getGeometry(eq_sys.parameters.get<std::string>("DO_file"));
   Real energy=0, normDO=0;
   const char* filename=eq_sys.parameters.get<std::string>("DO_file").c_str();
   int namelength=strlen(filename);
   getDyson(filename, namelength, do_j, l, alpha, energy, normDO);

   const FEType & fe_type = es1.get_dof_map().variable_type(var1);
   std::vector<dof_id_type> dof_indices;
   FEBase * cfe = libmesh_nullptr;

   // set correct k-vector
   eq_sys.parameters.set<Real>("current frequency")=sqrt(2.*
                     std::abs(eq_sys.parameters.get<Real>("energy")));

   // Begin the loop over the elements
   MeshBase::const_element_iterator       el     = es1.get_mesh().active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = es1.get_mesh().active_local_elements_end();
   for(; el != end_el; ++el){
       const Elem * elem = *el;
       const unsigned int dim = elem->dim();

      //QGauss qrule (dim, FIFTH);
      QGauss qrule (dim, fe_type.default_quadrature_order());
      UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
      UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
      fe->attach_quadrature_rule (&qrule);
      inf_fe->attach_quadrature_rule (&qrule);
      
      if (elem->infinite())
         cfe = inf_fe.get();
      else
         cfe = fe.get();
      const std::vector<Real> &  JxW     = cfe->get_JxW();
      const std::vector<Point> & q_point = cfe->get_xyz();

      cfe->reinit(elem);

      es1.get_dof_map().dof_indices (elem, dof_indices, var1);

      //const unsigned int n_qp = qrule.n_points(); --> fails for infinite elements.
      unsigned int n_qp = cfe->n_quadrature_points();
      const unsigned int n_sf = cast_int<unsigned int>(dof_indices.size());

      // Begin the loop over the Quadrature points.
      for (unsigned int qp=0; qp<n_qp; qp++){
         Number u_h = 0.;
         Number do_val=evalDO(do_j, l, alpha, geometry, q_point[qp]);
               
         Point mapped_qp = FEInterface::inverse_map(dim, fe_type, elem, q_point[qp], TOLERANCE, true); 
         FEComputeData fe_data(eq_sys, mapped_qp);
         FEInterface::compute_data(dim, fe_type, elem, fe_data);

         for (unsigned int i=0; i != n_sf; ++i){
            //Use FEComputeData because with infinite elements the value at q_point[i][qp]
            // is not just phi[i][qp].
            u_h += fe_data.shape[i] * (*local_v1)(dof_indices[i]);
            
         }
         //norm += JxW[qp] * TensorTools::norm_sq(u_h);
         if(int_type==MU)
            overlap += JxW[qp] * std::conj(u_h) * q_point[qp].norm() * do_val;
            //overlap += JxW[qp] * std::conj(do_val) * q_point[qp].norm() * do_val;
         else //if(int_type==OVERLAP)
            overlap += JxW[qp] * std::conj(u_h) * do_val;
            //overlap += JxW[qp] * std::conj(do_val) * do_val;
      }
   }

   es1.comm().sum(overlap);

   return overlap;
}

Number norm_DO(EquationSystems& eq_sys){
   //run at all processors at once:
   //parallel_object_only(); --> can not be used here.

   Number overlap=0;
   System & es = eq_sys.get_system<System> ("DO");

   // get the dyson orbital:
   std::vector<std::vector<double> > do_j;
   std::vector<unsigned int> l;
   std::vector<double> alpha;
   std::vector<Node> geometry= getGeometry(eq_sys.parameters.get<std::string>("DO_file"));
   Real energy=0, normDO=0;
   const char* filename=eq_sys.parameters.get<std::string>("DO_file").c_str();
   int namelength=strlen(filename);
   getDyson(filename, namelength, do_j, l, alpha, energy, normDO);

   const FEType & fe_type = es.get_dof_map().variable_type(0);
   FEBase * cfe = libmesh_nullptr;

   // Begin the loop over the elements
   MeshBase::const_element_iterator       el     = es.get_mesh().active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = es.get_mesh().active_local_elements_end();
   for(; el != end_el; ++el){
       const Elem * elem = *el;
       const unsigned int dim = elem->dim();

      QGauss qrule (dim, FIFTH);
      //QGauss qrule (dim, fe_type.default_quadrature_order());
      UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
      UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
      fe->attach_quadrature_rule (&qrule);
      inf_fe->attach_quadrature_rule (&qrule);
      
      if (elem->infinite())
         cfe = inf_fe.get();
      else
         cfe = fe.get();
      const std::vector<Real> &  JxW     = cfe->get_JxW();
      const std::vector<Point> & q_point = cfe->get_xyz();

      cfe->reinit(elem);

      //const unsigned int n_qp = qrule.n_points(); --> fails for infinite elements.
      unsigned int n_qp = cfe->n_quadrature_points();

      // Begin the loop over the Quadrature points.
      for (unsigned int qp=0; qp<n_qp; qp++){
         Number do_val=evalDO(do_j, l, alpha, geometry, q_point[qp]);
         overlap += JxW[qp] * TensorTools::norm_sq(do_val);

      }
   }

   es.comm().sum(overlap);

   return overlap;
}

Real normalise(EquationSystems& equation_systems, bool infel){
   CondensedEigenSystem& eigen_system=equation_systems.get_system<CondensedEigenSystem> ("EigenSE");
   LinearImplicitSystem & DO = equation_systems.get_system<LinearImplicitSystem> ("DO");
   //normalise eigen_system:
   Number norm_phi=0;
   //if (!infel)
   //   norm_phi=eigen_system.calculate_norm( *eigen_system.solution, 0, L2);
   //out<<"norm of phi: "<<norm_phi<<"   ";
   //out<< calculate_overlap(equation_systems, "EigenSE", 0, "EigenSE", 0, OVERLAP) <<std::endl;
   norm_phi=calculate_overlap(equation_systems, "EigenSE", 0, "EigenSE", 0, OVERLAP);

   //Real normDO = 0;
   //if (!infel) 
   //   normDO= DO.calculate_norm(*DO.solution, 0, L2);
   //out<<"norm of DO:   "<< normDO <<"  ";
   //out<< sqrt(calculate_overlap(equation_systems, "DO", 0, "DO", 0, OVERLAP))<<"  ";
   //out<< sqrt(norm_DO(equation_systems))<<std::endl;

   Number overlap=0;
   //compute <DO| mu |phi>
   //overlap=calculate_overlap(equation_systems, "DO", 0, "EigenSE", 0, MU);
   //out <<"solution overlap"<< overlap<<std::endl;
   overlap = overlap_DO(equation_systems, "EigenSE", 0, MU);
   //out <<"direct overlap"<< overlap<<std::endl;
   
   return abs(overlap)*abs(overlap)/(abs(norm_phi)*abs(norm_phi));
}
