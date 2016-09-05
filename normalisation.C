#include <vector>
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
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/quadrature.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_tools.h"


using namespace libMesh;

enum IntegralType: int{
   MU=0,
   OVERLAP
};

Number calculate_overlap(System& es1, const NumericVector<Number>& vect1, int var1, System& es2, const NumericVector<Number>& vect2, int var2, IntegralType int_type ){
   //run at all processors at once:
   //parallel_object_only(); --> can not be used here.

   Number overlap=0;
   if (int_type != MU){
      libmesh_not_implemented();
   }
   else{
      //libmesh_assert_not_equal_to(es1.comm(), es2.comm());
      //libmesh_assert_not_equal_to(es1.get_dof_map().variable_type(var1),
      //                            es2.get_dof_map().variable_type(var2));

      // Localize the potentially parallel vectors
      UniquePtr<NumericVector<Number> > local_v1 = NumericVector<Number>::build(es1.comm());
      local_v1->init(vect1.size(), true, SERIAL);
      vect1.localize (*local_v1, es1.get_dof_map().get_send_list());
      UniquePtr<NumericVector<Number> > local_v2 = NumericVector<Number>::build(es2.comm());
      local_v2->init(vect2.size(), true, SERIAL);
      vect2.localize (*local_v2, es2.get_dof_map().get_send_list());

      const FEType & fe_type = es2.get_dof_map().variable_type(var1);
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

      // Begin the loop over the elements
      MeshBase::const_element_iterator       el     =
            es1.get_mesh().active_local_elements_begin();
      const MeshBase::const_element_iterator end_el =
            es1.get_mesh().active_local_elements_end();
      for ( ; el != end_el; ++el){
          const Elem * elem = *el;
          const unsigned int dim = elem->dim();


         FEBase * fe = fe_ptrs[dim];
         QBase * qrule = q_rules[dim];
         libmesh_assert(fe);
         libmesh_assert(qrule);

         const std::vector<Real> &               JxW = fe->get_JxW();
         const std::vector<Point> &              q_point = fe->get_xyz();
         const std::vector<std::vector<Real> > * phi = libmesh_nullptr;
         //const std::vector<std::vector<Real> > * phi = fe->get_phi();

         fe->reinit (elem);
         es1.get_dof_map().dof_indices (elem, dof_indices, var1);

         const unsigned int n_qp = qrule->n_points();

         const unsigned int n_sf = cast_int<unsigned int>(dof_indices.size());

         // Begin the loop over the Quadrature points.
         for (unsigned int qp=0; qp<n_qp; qp++){
            Number u_h = 0.;
            Number v_h = 0.;
            for (unsigned int i=0; i != n_sf; ++i){
               u_h += (*phi)[i][qp] * (*local_v1)(dof_indices[i]);
               v_h += (*phi)[i][qp] * (*local_v2)(dof_indices[i]);
               
               //overlap += JxW[qp] * TensorTools::norm_sq(u_h);
               if(int_type==MU)
                  overlap += JxW[qp] * u_h*q_point[qp].norm()*v_h;
               else if(int_type==OVERLAP)
                  overlap += JxW[qp] * u_h*v_h;
            }

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
   }

   es1.comm().sum(overlap);
   overlap = std::sqrt(overlap);

   return overlap;
}

Real normalise(EquationSystems& equation_systems){
   CondensedEigenSystem& eigen_system=equation_systems.get_system<CondensedEigenSystem> ("EigenSE");
   LinearImplicitSystem & DO = equation_systems.get_system<LinearImplicitSystem> ("DO");
   //normalise eigen_system:
   double norm_phi=eigen_system.calculate_norm( *eigen_system.solution, 0, L2);
   out<<norm_phi<<"  ";
   out<<DO.calculate_norm(*DO.solution, 0, L2)<<std::endl;

   //compute <DO| mu |phi>
   Number overlap=calculate_overlap(eigen_system, *eigen_system.solution, 0,DO, *DO.solution, 0, MU);
   
   return abs(overlap)*abs(overlap)/(norm_phi*norm_phi);
}
