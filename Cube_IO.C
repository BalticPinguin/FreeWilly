#include <complex.h>
#include <iostream>
#include <fstream>
// libMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/tree.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/fe.h"
#include "libmesh/dof_map.h"
#include "libmesh/inf_fe.h"
#include "libmesh/fe_interface.h" 
#include "libmesh/fe_compute_data.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void cube_io(EquationSystems& es, std::vector<Node> geom, std::string output){
   CondensedEigenSystem & system = es.get_system<CondensedEigenSystem> ("EigenSE"); // --> how to generalise??
   const MeshBase & mesh = es.get_mesh();
   const DofMap & dof_map = system.get_dof_map();
   
   UniquePtr<NumericVector<Number> > solution_vect = 
        NumericVector<Number>::build(es.comm());

   solution_vect->init((*system.solution).size(), true, SERIAL);
   (*system.solution).localize(* solution_vect);
   
   const FEType & fe_type = dof_map.variable_type(0);
   UniquePtr<FEBase> fe (FEBase::build(3, fe_type));
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(3, fe_type));
   FEBase * cfe = libmesh_nullptr;
   QGauss qrule (3, SECOND);
   std::vector<dof_id_type> dof_indices;
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);

   // set output to filename
   std::ofstream out(output);
   out<<"EigenSE"<<std::endl<<std::endl; // print first two lines: comments
   out<<std::setw(5)<<"  "<<geom.size();
   // where do I start?
   Point start(-5.5,-5.5,-5.5);
   out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   // print # points per axis and step in Cartesian Coordinates:

   Real dx=0.283450;
   Real dy=0.283450;
   Real dz=0.283450;
   unsigned int nx=40;
   unsigned int ny=40;
   unsigned int nz=40;

   out<<std::setw(5)<<nx;
   out<<std::setw(12)<<std::setprecision(6)<<" \t"<<dx<<" \t 0.00000 \t0.00000"<<std::endl;
   out<<std::setw(5)<<ny;
   out<<std::setw(12)<<std::setprecision(6)<<" \t 0.00000"<<dy<<" \t 0.00000"<<std::endl;
   out<<std::setw(5)<<nz;
   out<<std::setw(12)<<std::setprecision(6)<<" \t 0.00000 \t0.00000 \t "<<dz<<std::endl;

   for(unsigned int i=0; i<geom.size(); i++){
      //out<<std::setw(5)<<geom[i].id()<<"\t";
      out<<std::setw(5)<<"1.00000"<<"\t";
      out<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      out<<std::setw(12)<<std::setprecision(6)<<geom[i](0)<<"\t";
      out<<std::setw(12)<<std::setprecision(6)<<geom[i](1)<<"\t";
      out<<std::setw(12)<<std::setprecision(6)<<geom[i](2)<<"\n";
   }

   unsigned int ix, iy, iz;
   PointLocatorTree pt_lctr(mesh);
   //pt_lctr.enable_out_of_mesh_mode();
   //pt_lctr.init(); 
   for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
         for (iz=0;iz<nz;iz++) {

            Point q_point(start(0)+(Real)ix*dx,
                          start(1)+(Real)iy*dy,
                          start(2)+(Real)iz*dz);
            
            const Elem * elem=pt_lctr(q_point);
            if(elem==NULL){
               out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<0.0;
            }
            else{

               dof_map.dof_indices (elem, dof_indices);
   
               Point map_point=FEInterface::inverse_map(3, fe_type, elem, q_point, TOLERANCE, true); 
               FEComputeData data(es, map_point); 
               if (elem->infinite()){
                  out<<"point: "<<q_point<<"   ";
                  out<<map_point(2)<<std::endl;
               }
               FEInterface::compute_data(3, fe_type, elem, data);
            
               //compute solution value at that point.
               Number soln=0;
               if (elem->infinite())
                  cfe = inf_fe.get();
               else
                  cfe = fe.get();
               cfe->reinit(elem);
               unsigned int n_sf= cfe->n_shape_functions();
               for (unsigned int i=0; i<n_sf; i++){
                  soln+=(*solution_vect)(dof_indices[i])*data.shape[i];
               }
               out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::real(soln);
            }

            if (iz % 6 == 5)
               out<<std::endl;
         }
         out<<std::endl;
      }
   }
}
