#include <cstring> // for strlen.

#include "read_DO.h"
#include "FreeWilly.h"

#include "fsu_soft/besselj.hpp"

// includes for calculate_norm, point_*
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_tools.h"

using namespace libMesh;

// Functions to get intensities and norms
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
   //UniquePtr<NumericVector<Number> > local_v1 = NumericVector<Number>::build(es1.comm());
   //local_v1->init((*es1.solution).size(), true, SERIAL);
   //(*es1.solution).localize (*local_v1, es1.get_dof_map().get_send_list());
   //UniquePtr<NumericVector<Number> > local_v2 = NumericVector<Number>::build(es2.comm());
   //local_v2->init((*es2.solution).size(), true, SERIAL);
   //(*es2.solution).localize (*local_v2, es2.get_dof_map().get_send_list());
   const std::string & formulation = eq_sys.parameters.get<std::string>("formulation");
   Real power=eq_sys.parameters.get<Real> ("power");

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
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D

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
	    if(formulation=="symmetric"){
       	       u_h += sqrt(weight[qp])* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	       v_h += sqrt(weight[qp])* fe_data.shape[i] * (*es2.solution)(dof_indices[i]);
	    }
	    else if(formulation=="root"){
       	       u_h += sqrt(sqrt(weight[qp]))* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	       v_h += sqrt(sqrt(weight[qp]))* fe_data.shape[i] * (*es2.solution)(dof_indices[i]);
            }
	    else if(formulation=="power"){
       	       u_h += pow(weight[qp],power)* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	       v_h += pow(weight[qp],power)* fe_data.shape[i] * (*es2.solution)(dof_indices[i]);
            }
            else{
       	       u_h += fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	       v_h += fe_data.shape[i] * (*es2.solution)(dof_indices[i]);
	    }
         }
         //out<<"later: "<<q_point[qp];
         //out<<"  "<<u_h<<std::endl;

         //norm += JxW[qp] * TensorTools::norm_sq(u_h);
         if(int_type==LENGTH)
            // this expression is wrong:
            overlap += JxW[qp] * std::conj(u_h) * q_point[qp].norm() * v_h* 0.;
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

Real overlap_DO(EquationSystems& eq_sys, const std::string sys1, int var1, IntegralType int_type, bool infinite){
   //run at all processors at once:
   //parallel_object_only(); --> can not be used here.

   Number overlap_x=0;
   Number overlap_y=0;
   Number overlap_z=0;
   Number overlap=0;
   
   // this should be checked somehow as well:
   //libmesh_assert_not_equal_to(es1.comm(), es2.comm());
   //libmesh_assert_not_equal_to(es1.get_dof_map().variable_type(var1),
   //                            es2.get_dof_map().variable_type(var2));
   //CondensedEigenSystem& es1=equation_systems.get_system<CondensedEigenSystem> (sys1);
   //LinearImplicitSystem & es2 = equation_systems.get_system<LinearImplicitSystem> (sys2);
   System & es1 = eq_sys.get_system<System> (sys1);

   // Localize the potentially parallel vectors
   //UniquePtr<NumericVector<Number> > local_v1 = NumericVector<Number>::build(es1.comm());
   //local_v1->init((*es1.solution).size(), true, SERIAL);
   //(*es1.solution).localize (*local_v1, es1.get_dof_map().get_send_list());

   // get the dyson orbital:
   DOrbit dyson (eq_sys.parameters.get<std::string>("DO_file"));
   const std::string & formulation = eq_sys.parameters.get<std::string>("formulation");
   Real power=eq_sys.parameters.get<Real> ("power");

   const FEType & fe_type = es1.get_dof_map().variable_type(var1);
   std::vector<dof_id_type> dof_indices;
   FEBase * cfe = libmesh_nullptr;

   // Begin the loop over the elements
   MeshBase::const_element_iterator       el     = es1.get_mesh().active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = es1.get_mesh().active_local_elements_end();
   for(; el != end_el; ++el){
       const Elem * elem = *el;
       const unsigned int dim = elem->dim();
   
      // skip infinite elements if inf
      if (!infinite && elem->infinite())
         continue;

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
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D

      cfe->reinit(elem);

      es1.get_dof_map().dof_indices (elem, dof_indices, var1);

      //const unsigned int n_qp = qrule.n_points(); --> fails for infinite elements.
      unsigned int n_qp = cfe->n_quadrature_points();
      const unsigned int n_sf = cast_int<unsigned int>(dof_indices.size());

      // Begin the loop over the Quadrature points.
      for (unsigned int qp=0; qp<n_qp; qp++){
         Number u_h = 0.;
         Number do_val=dyson.evalDO(q_point[qp]);
               
         Point mapped_qp = FEInterface::inverse_map(dim, fe_type, elem, q_point[qp], TOLERANCE, true); 
         FEComputeData fe_data(eq_sys, mapped_qp);
         FEInterface::compute_data(dim, fe_type, elem, fe_data);

         for (unsigned int i=0; i != n_sf; ++i){
            //Use FEComputeData because with infinite elements the value at q_point[i][qp]
            // is not just phi[i][qp].
	    if(formulation=="symmetric")
		u_h += sqrt(weight[qp])* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	    else if(formulation=="root")
       	       u_h += sqrt(sqrt(weight[qp]))* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	    else if(formulation=="power")
       	       u_h += pow(weight[qp],power)* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
            else
		u_h += fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
            
         }
         //norm += JxW[qp] * TensorTools::norm_sq(u_h);
         if(int_type==LENGTH)
         {
            overlap_x += JxW[qp] * std::conj(u_h) * q_point[qp](0) * do_val;
            overlap_y += JxW[qp] * std::conj(u_h) * q_point[qp](1) * do_val;
            overlap_z += JxW[qp] * std::conj(u_h) * q_point[qp](2) * do_val;
            //overlap += JxW[qp] * std::conj(do_val) * q_point[qp].norm() * do_val;
         }
         else //if(int_type==OVERLAP)
            overlap+= JxW[qp] * std::conj(u_h) * do_val;
            //overlap += JxW[qp] * std::conj(do_val) * do_val;
      }
   }
   
   if(int_type==LENGTH)
      {
      overlap=overlap_x*conj(overlap_x)+
              overlap_y*conj(overlap_y)+
              overlap_z*conj(overlap_z);
     }
   else if (int_type==VELOCITY)
      {
      overlap=overlap_x*conj(overlap_x)+
              overlap_y*conj(overlap_y)+
              overlap_z*conj(overlap_z) ;
     }
   else
      overlap*=conj(overlap);

   es1.comm().sum(overlap);
   
   // abs is needed here to avoid compiler errors.
   return std::abs(overlap);
}

Number norm_DO(EquationSystems& eq_sys, bool infinite){
   //run at all processors at once:
   //parallel_object_only(); --> can not be used here.

   Number overlap=0;
   System & es = eq_sys.get_system<System> ("EigenSE");

   // get the dyson orbital:
   DOrbit dyson (eq_sys.parameters.get<std::string>("DO_file"));

   const FEType & fe_type = es.get_dof_map().variable_type(0);
   FEBase * cfe = libmesh_nullptr;

   // Begin the loop over the elements
   MeshBase::const_element_iterator       el     = es.get_mesh().active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = es.get_mesh().active_local_elements_end();
   for(; el != end_el; ++el){
       const Elem * elem = *el;
       const unsigned int dim = elem->dim();
   
      // skip infinite elements if inf
      if (!infinite && elem->infinite())
         continue;

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
         Number do_val=dyson.evalDO(q_point[qp]);
         overlap += JxW[qp] * TensorTools::norm_sq(do_val);

      }
   }

   es.comm().sum(overlap);

   return overlap;
}

Real normalise(EquationSystems& equation_systems, bool infel){
   //normalise eigen_system:
   //Number norm_phi=calculate_overlap(equation_systems, "EigenSE", 0, "EigenSE", 0, OVERLAP);
   //out<<"norm of phi: ";
   //out<< sqrt(norm_phi) <<std::endl;

   //compute <DO| mu |phi>
   Real overlap = overlap_DO(equation_systems, "EigenSE", 0, LENGTH, infel);
   
   // abs needed for type conversion.
   return overlap; ///std::abs((norm_phi*conj(norm_phi)));
}

// functions used to project solution onto spherical waves
void Y_lm(Real x, Real y, Real z, int l_max, std::vector<Number>& solution){
   //http://www.ppsloan.org/publications/StupidSH36.pdf
   //12.5663706144=4*pi
   solution.resize((l_max+1)*(l_max+1));
   Real K_lm;
   
   Real** value= new Real*[l_max+1];

   // memory for value is allocated in pm_polynomial_value()!
   //for(int m=0; m<=l_max; m++)
   //   value[m]=new Real[l_max+1];

   double r=sqrt(x*x+y*y+z*z),
            phival=0;
   double theta[1];

   if( r<1e-12){
      theta[0]=0;
      phival=0;
   }
   else{
      theta[0]=z/r;
      phival=atan2(y,x); //-pi -> pi
   }

   //value = p_polynomial_value(1, l, theta );
   // value is not normalised!
   for(int m=0; m<=l_max; m++){
      // this is against the convention in QM: for m>0 && m%2==1
      // it should have the negative of it; I don't care at this point.
      value[m] = pm_polynomial_value ( 1, l_max, m, theta);
   }
   int k=0;
   for(int l=0; l<=l_max; l++){
      for(int m=-l; m<=l; m++){
         K_lm=sqrt((2*l+1)*factorial(l-abs(m))/(12.5663706144*factorial(l+abs(m))));
         solution[k]=value[abs(m)][l]*K_lm*Number(cos(m*phival),sin(m*phival));
         k++;
      }
   }

   for(int m=0; m<=l_max; m++)
      delete [] value[m];

   delete [] value;
}

void evalSphWave(int l_max, Point qp, Real k, std::vector<Number>& solution){
   int error;
   Number wave;
   Real* R = new Real[l_max+1];
   solution.resize((l_max+1)*(l_max+1));
   
   std::vector<Number> angular;
   //angular.resize((l_max+1)*(l_max+1));

   Real kr=qp.norm()*k;
   rjbesl(kr, 0 ,l_max+1, R, error);
   assert(error>=0);
   if (error!=l_max+1){
      err<<"The evaluation of bessel functions returned"<<std::endl;
      err<<"   "<<error<<std::endl;
      if ( 0 > qp.norm()*k )
         err<<"   qp.norm*k < 0"<<std::endl;
      if ( qp.norm()*k > 1e+04 )
         err<<"   qp.norm*k > x_max "<<std::endl;
   }

   Y_lm(qp(0), qp(1), qp(2), l_max, angular);

   int j=0;
   for(int l=0; l<=l_max; l++){
      for(int m=-l; m<=l; m++){
         if(error< l)
           //when the besselfunction did not converge.
           solution[j]=0;
         else
            solution[j]=R[l]*angular[j];
            if(solution[j]!=solution[j]){
              err<<"!! "<<R[l]<<"  "<<angular[j]<<"  "<<j<<std::endl;
            }
         j++;
      }
   }
   delete [] R;
}

void projection(EquationSystems& es, const std::string sys, int l_max, std::vector<Number>& overlap ){
   // this should be checked somehow as well:
   //libmesh_assert_not_equal_to(es1.comm(), es2.comm());
   //libmesh_assert_not_equal_to(es1.get_dof_map().variable_type(var1),
   //                            es2.get_dof_map().variable_type(var2));
   //CondensedEigenSystem& es1=equation_systems.get_system<CondensedEigenSystem> (sys1);
   //LinearImplicitSystem & es2 = equation_systems.get_system<LinearImplicitSystem> (sys2);
   System & es1 = es.get_system<System> (sys);

   std::vector<Number> spherical_qp;
   //std::vector<Number> overlap;
   std::vector<Number> normSphere;
   overlap.resize((l_max+1)*(l_max+1));
   normSphere.resize((l_max+1)*(l_max+1));

   Number norm=0;
   //Number norm2=0;
    
   Real k = std::abs(es.parameters.get<Number>("momentum"));
   const std::string & formulation = es.parameters.get<std::string>("formulation");
   Real power=es.parameters.get<Real> ("power");

   // Localize the potentially parallel vectors
   //UniquePtr<NumericVector<Number> > local_v1 = NumericVector<Number>::build(es1.comm());
   //local_v1->init((*es1.solution).size(), true, SERIAL);
   //(*es1.solution).localize (*local_v1, es1.get_dof_map().get_send_list());

   const FEType & fe_type = es1.get_dof_map().variable_type(0);
   std::vector<dof_id_type> dof_indices;
   FEBase * cfe = libmesh_nullptr;

   // Begin the loop over the elements
   MeshBase::const_element_iterator       el     = es1.get_mesh().active_local_elements_begin();
   const MeshBase::const_element_iterator end_el = es1.get_mesh().active_local_elements_end();
   for(; el != end_el; ++el){
       const Elem * elem = *el;
       const unsigned int dim = elem->dim();

      // skip infinite elements
      if (elem->infinite())
         continue;

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
      const std::vector<Real>& weight = cfe->get_Sobolev_weight(); // in publication called D

      cfe->reinit(elem);

      es1.get_dof_map().dof_indices (elem, dof_indices, 0);

      //const unsigned int n_qp = qrule.n_points(); --> fails for infinite elements.
      unsigned int n_qp = cfe->n_quadrature_points();
      const unsigned int n_sf = cast_int<unsigned int>(dof_indices.size());

      // Begin the loop over the Quadrature points.
      for (unsigned int qp=0; qp<n_qp; qp++){
         Number u_h = 0.;
         evalSphWave(l_max, q_point[qp], k, spherical_qp);
         
         Point mapped_qp = FEInterface::inverse_map(dim, fe_type, elem, q_point[qp], TOLERANCE, true); 
         FEComputeData fe_data(es, mapped_qp);
         FEInterface::compute_data(dim, fe_type, elem, fe_data);

         for (unsigned int i=0; i != n_sf; ++i){
            //Use FEComputeData because with infinite elements the value at q_point[i][qp]
            // is not just phi[i][qp].
	    if(formulation=="symmetric")
		u_h += sqrt(weight[qp])* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	    else if(formulation=="root")
       	       u_h += sqrt(sqrt(weight[qp]))* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
	    else if(formulation=="power")
		u_h += pow(weight[qp],power)* fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
            else
		u_h += fe_data.shape[i] * (*es1.solution)(dof_indices[i]);
         }
         int j=0;
         for(int l=0; l<=l_max; l++){
            for(int m=-l; m<=l; m++){
               overlap[j]+= JxW[qp] * std::conj(u_h) * spherical_qp[j];
               normSphere[j]+=JxW[qp] * std::conj(spherical_qp[j]) * spherical_qp[j];
               j++;
            }
         }
         norm += JxW[qp] * TensorTools::norm_sq(u_h);
      }
   }
   int j=0;
   for(int l=0; l<=l_max; l++){
      for(int m=-l; m<=l; m++){
         overlap[j]=overlap[j]/(sqrt(normSphere[j])*sqrt(norm));
         es1.comm().sum(overlap[j]);
         j++;
      }
   }
   //return overlap;
}

void ProjectSphericals (EquationSystems& es, int l_max, int /*i*/){
   Number norm_phi=calculate_overlap(es, "EigenSE", 0, "EigenSE", 0, OVERLAP);
   norm_phi=sqrt(norm_phi*conj(norm_phi));
   Number tot_proj=0, this_proj;
   int m;
   libmesh_assert_greater(l_max,0);
   std::vector<Number> contributions;
   projection(es,"EigenSE", l_max, contributions);
   // make some nice output:
   out<<"====================================";
   out<<std::endl;
   int k=0;
   for( int l=0; l<=l_max; l++){
      for(m=-l; m<=l; m++){
         this_proj=contributions[k];
         tot_proj+=this_proj*conj(this_proj);
         out<<"|     l = "<<l<<"       ";
         out<<"        m = "<<m<<"       ";
         k++;
         out<<"  \t"<<this_proj<<" ";
         out<<"\t|"<<std::endl;
      }
      out<<std::endl;
   }
   out<<"||  Total:  ";
   out<<sqrt(tot_proj)<<"   |"<<std::endl;
   out<<"====================================";
   out<<std::endl<<std::endl;
}
