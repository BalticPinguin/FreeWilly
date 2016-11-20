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

using namespace libMesh;

Number evalSphWave(int l, int m, Point qp, Real k);
void cube_sphere(EquationSystems& es, std::string output, int l, int m);

void PlotSphericals (EquationSystems& es, int l_max){
   Number tot_proj=0, this_proj;
   int m;
   libmesh_assert_greater(l_max,0);
   for( int l=0; l<=l_max; l++){
      for(m=-l; m<=l; m++){
         std::ostringstream output;
         output<<"-"<<l<<"_"<<m<<"_.cube";
         cube_sphere(es, output.str(), l, m);
      }
   }
}

void cube_sphere(EquationSystems& es, std::string output, int l, int m){
   //CondensedEigenSystem & system = es.get_system<CondensedEigenSystem> ("EigenSE"); // --> how to generalise??
   System & system = es.get_system<System> ("EigenSE");
   const DofMap & dof_map = system.get_dof_map();
   
   const FEType & fe_type = dof_map.variable_type(0);
   UniquePtr<FEBase> fe    (FEBase::build(3, fe_type));
   UniquePtr<FEBase> inf_fe(FEBase::build_InfFE(3, fe_type));
   QGauss qrule (3, fe_type.default_quadrature_order());
   std::vector<dof_id_type> dof_indices;
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);

   // set output to filename
   std::ostringstream re_output;
   re_output<<"re_"<<output;
   std::ostringstream im_output;
   im_output<<"im_"<<output;
   std::ostringstream abs_output;
   abs_output<<"abs_"<<output;
   std::ofstream re_out(re_output.str());
   std::ofstream im_out(im_output.str());
   std::ofstream abs_out(abs_output.str());
   re_out<<"deine Mudda"<<std::endl<<std::endl; // print first two lines: comments
   im_out<<"deine Mudda"<<std::endl<<std::endl; 
   abs_out<<"deine Mudda"<<std::endl<<std::endl;

   re_out<<std::setw(5)<<"  "<<1;
   im_out<<std::setw(5)<<"  "<<1;
   abs_out<<std::setw(5)<<"  "<<1;
   // where do I start?
   Point mol_center(0,0,0);

   Real r = 2.*es.parameters.get<Real>("radius");
   Real lambda = 1./es.parameters.get<Real>("current frequency");

   Real dx=lambda/6.;
   Real dy=lambda/6.;
   Real dz=lambda/6.;
   unsigned int nx=(2*r)/dx;
   unsigned int ny=(2*r)/dy;
   unsigned int nz=(2*r)/dz;

   Point start(mol_center(0)-dx*nx/2.,
               mol_center(1)-dy*ny/2.,
               mol_center(2)-dz*nz/2.);

   re_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   re_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   re_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   im_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   im_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   im_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   abs_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(0);
   abs_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(1);
   abs_out<<std::setw(12)<<std::setprecision(6)<<"   "<<start(2)<<std::endl;
   // print # points per axis and step in Cartesian Coordinates:

   re_out<<std::setw(5)<<nx;
   re_out<<std::setw(12)<<std::setprecision(5)<<" \t "<<dx<<" \t\t 0.00000 \t 0.00000"<<std::endl;
   re_out<<std::setw(5)<<ny;
   re_out<<std::setw(12)<<std::setprecision(5)<<" \t\t 0.00000 \t "<<dy<<" \t\t 0.00000"<<std::endl;
   re_out<<std::setw(5)<<nz;
   re_out<<std::setw(12)<<std::setprecision(5)<<" \t\t 0.00000 \t 0.00000 \t "<<dz<<std::endl;
   im_out<<std::setw(5)<<nx;
   im_out<<std::setw(12)<<std::setprecision(5)<<" \t "<<dx<<" \t\t 0.00000 \t 0.00000"<<std::endl;
   im_out<<std::setw(5)<<ny;
   im_out<<std::setw(12)<<std::setprecision(5)<<" \t\t 0.00000 \t "<<dy<<" \t\t 0.00000"<<std::endl;
   im_out<<std::setw(5)<<nz;
   im_out<<std::setw(12)<<std::setprecision(5)<<" \t\t 0.00000 \t 0.00000 \t "<<dz<<std::endl;
   abs_out<<std::setw(5)<<nx;
   abs_out<<std::setw(12)<<std::setprecision(5)<<" \t "<<dx<<" \t\t 0.00000 \t 0.00000"<<std::endl;
   abs_out<<std::setw(5)<<ny;
   abs_out<<std::setw(12)<<std::setprecision(5)<<" \t\t 0.00000 \t "<<dy<<" \t\t 0.00000"<<std::endl;
   abs_out<<std::setw(5)<<nz;
   abs_out<<std::setw(12)<<std::setprecision(5)<<" \t\t 0.00000 \t 0.00000 \t "<<dz<<std::endl;

   {
      re_out<<" "<<std::setw(5)<<1<<"\t";
      im_out<<" "<<std::setw(5)<<1<<"\t";
      abs_out<<" "<<std::setw(5)<<1<<"\t";
      //out<<std::setw(5)<<"1.00000"<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      re_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\n";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      im_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\n";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\t";
      abs_out<<" "<<std::setw(12)<<std::setprecision(6)<<"0.00000"<<"\n";
   }

   unsigned int ix, iy, iz;
   unsigned int num_line=0;
   Number soln;
   Real k = es.parameters.get<Real>("current frequency")*2.*pi;
   for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
         for (iz=0;iz<nz;iz++) {

            num_line++;
            Point q_point(start(0)+(Real)ix*dx,
                          start(1)+(Real)iy*dy,
                          start(2)+(Real)iz*dz);
            soln=evalSphWave(l,  m, q_point, k);

            re_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::real(soln);
            im_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::imag(soln);
            abs_out<<" "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::abs(soln);
         }
      }
   }
}
