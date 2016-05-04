#include <vector>
#include <fstream>
#include "libmesh/point.h" // to have access to "point"-object
#include "libmesh/mesh.h"  // for 'mesh'
#include "libmesh/elem.h"
#include "libmesh/cell_hex8.h"
// for GetPotential:
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/explicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h" // Define Gauss quadrature rules.
#include <cmath> // needed for abs().
// for printing the equation system
#include "libmesh/exodusII_io.h"

#include <libmesh/mesh_function.h>

using namespace libMesh;

//function to read the data from file and feed a mesh with respective nodes. After this, a .rhs is taken and the values are added:

struct ESP{
   // these vectors store the points and potential (at the points) as given.
   std::vector<Point> node;
   std::vector<double> potential;
   // k-th node opens an element, if there are 7 elements with distance <=3^1/2*min_dist
   // whose components are >= the components of k-th node.
   std::vector<std::vector<unsigned int>> neighbour;
   std::vector<bool> used; 
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

void FindNeighbours(ESP& esp){
   // idea of the function: find neighbours, always looking back (smaller coordinate values)
   // the strategy in different directions is different:
   // x-direction: if previous element is close enough: accept it as neighbour.
   // y-direction: try to find element in last row at same position and check elements 
   //        with indices 1 larger and 1 smaller as well.
   // z-direction: try to find element that is below the current one (or closest to it)
   //        and test all of its neighbours as own neighbours.

   unsigned int y=0, z=0, y_new=0, z_new=0, elem_below;
   int offset_y, offset_z;
   esp.neighbour.resize(esp.size);
   // this will fail when the first two elements are not neighours.
   // I am quite sure this will never happen but this can cause trouble!
   double THRESHOLD=0.0001;
   double step=esp.node[1](0)-esp.node[0](0)+THRESHOLD; // add some tolerance
   for(unsigned int i=1; i<esp.size; i++){
      if(esp.node[i](0)-esp.node[i-1](0)<=step && esp.node[i](0)-esp.node[i-1](0)> 0){
         esp.neighbour[i].push_back(i-1);
         esp.neighbour[i-1].push_back(i);
      }
      if (esp.node[i](1)!=esp.node[i-1](1)){
         offset_y=0; 
         // y is the start-value of the last row.
         y=y_new;
         // y_new is the start-value of this row.
         y_new=i;
      }
      if (esp.node[i](2)!=esp.node[i-1](2)){
         offset_z=0;
         // z is the start-value of the last row.
         z=z_new;
         // z_new is the start-value of this row.
         z_new=i;
      }

      // there is no neighbour with lower y-value in this row if I just changed the
      // z-value (and hence are in the row with smallest y-values!)
      if (y_new!=z_new){
         // there might be a whole! it might be in this or in previous layer...
         if (std::abs(esp.node[i](0)-esp.node[y+offset_y+(i-y_new)](0))>THRESHOLD){
            // so adjust offset_y accordingly and move on!
            // the last line started more left of it:
            while (esp.node[i](0)>esp.node[y+offset_y+i-y_new](0) 
                     && esp.node[i](1)>esp.node[y+offset_y+i-y_new+1](1) +THRESHOLD){
               offset_y++;
               //out<<"      "<<offset_y<<std::endl;
            }
           // while ((esp.node[i](0)<esp.node[y+offset_y+i-y_new](0) && (i-y_new+offset_y)>0)
           //             || esp.node[i](1)-esp.node[y+offset_y+i-y_new-1](1)<THRESHOLD){
            while ( (i-y_new+offset_y)>0
                   && (esp.node[i]-esp.node[y+offset_y+i-y_new]).size()
                   >(esp.node[i]-esp.node[y+offset_y+i-y_new-1]).size()){
               offset_y--;
               //out<<"             "<<offset_y<<std::endl;
            }
            // if there is still an offset: there is no nearest neighbour in this direction.
            if (std::abs(esp.node[i](0)-esp.node[y+offset_y+(i-y_new)](0))>step)
               continue;
         }
         if (std::abs(esp.node[i](0)-esp.node[y+offset_y+(i-y_new)](0))<=step
                && std::abs(esp.node[i](1)-esp.node[y+offset_y+(i-y_new)](1))<=step){
            // the stepping in y-direction is asured to be exactly 'step'.
            esp.neighbour[i].push_back(y+offset_y+i-y_new);
            esp.neighbour[y+offset_y+i-y_new].push_back(i);
         }
         // get the diagonal elements as well
         if (i-y_new+offset_y >0 && std::abs(esp.node[i](0)-esp.node[y+offset_y+i-y_new-1](0))<=step
                       && std::abs(esp.node[i](1)-esp.node[y+offset_y+i-y_new-1](1))<=step){
            // the correct step-size in y- and z-direction is assured.
            esp.neighbour[i].push_back(y+offset_y-1+i-y_new);
            esp.neighbour[y+offset_y-1+i-y_new].push_back(i);
         }
         if (std::abs(esp.node[i](0)-esp.node[y+offset_y+i-y_new+1](0))<=step
                       && std::abs(esp.node[i](1)-esp.node[y+offset_y+i-y_new+1](1))<=step
                       && y+offset_y+i-y_new <y_new ){
            // the correct step-size in y- and z-direction is assured.
            esp.neighbour[i].push_back(y+offset_y+1+i-y_new);
            esp.neighbour[y+offset_y+1+i-y_new].push_back(i);
         }
      }
      if (z_new>0){
         // even in 3D, only the variable that changes fastest needs to be considered explicitly.
         if (std::abs(esp.node[i](0)-esp.node[z+offset_z+i-z_new](0))>THRESHOLD
               || std::abs(esp.node[i](1)-esp.node[z+offset_z+i-z_new](1))>THRESHOLD){
            // for efficiency: check first, if the previous element fits better:
            if (z+offset_z+i-z_new>0 
                  && (esp.node[i]-esp.node[z+offset_z+i-z_new]).size()
                    >(esp.node[i]-esp.node[z+offset_z+i-z_new-1]).size()+THRESHOLD)
               offset_z--;
            // so adjust offset_z accordingly and move on. Therefore, search neighbours
            // that are closer to the current element.
            if (std::abs(esp.node[i](0)-esp.node[z+offset_z+i-z_new](0))>THRESHOLD
                  || std::abs(esp.node[i](1)-esp.node[z+offset_z+i-z_new](1))>THRESHOLD){
               double oldDist=(esp.node[i]-esp.node[z+offset_z+i-z_new]).size();
               for (unsigned int layer_below=0; 
                        layer_below<esp.neighbour[z+offset_z+i-z_new].size(); 
                        layer_below++){
                  if (oldDist>(esp.node[i]-esp.node[esp.neighbour[z+offset_z+i-z_new][layer_below]]).size() 
                     &&(esp.node[i](2)-esp.node[esp.neighbour[z+offset_z+i-z_new][layer_below]](2)>THRESHOLD)){

                     // readjust offset_z and start the search in its neighbours from the beginning.
                     offset_z=esp.neighbour[z+offset_z+i-z_new][layer_below]-z-i+z_new;
                     oldDist=(esp.node[i]-esp.node[z+offset_z+i-z_new]).size();
                     layer_below=0;
                  }
               }
            }
         }
         // if this is reached now:
         elem_below=z+offset_z+i-z_new;
         if (std::abs(esp.node[i](0)-esp.node[elem_below](0))<=step
               && std::abs(esp.node[i](1)-esp.node[elem_below](1))<=step
               && std::abs(esp.node[i](2)-esp.node[elem_below](2))<=step){
            // now, search for other neighbours in the layer below: therefore, search all neigbours
            // of the element below this one:
            for(unsigned int layer_below=0; 
               layer_below<esp.neighbour[elem_below].size(); 
               layer_below++){
               // the first two conditions might be tautologies.
               if (std::abs(esp.node[i](0)-esp.node[esp.neighbour[elem_below][layer_below]](0))<=step 
                       && std::abs(esp.node[i](1)-esp.node[esp.neighbour[elem_below][layer_below]](1))<=step
                       && std::abs(esp.node[elem_below](2)-esp.node[esp.neighbour[elem_below][layer_below]](2))
                             <=THRESHOLD){
                  esp.neighbour[i].push_back(esp.neighbour[elem_below][layer_below]);
                  esp.neighbour[esp.neighbour[elem_below][layer_below]].push_back(i);
                  //out<<esp.neighbour[elem_below][layer_below]<<std::endl;
               }
            }
            // the correct step-size in y- and z-direction is assured.
            esp.neighbour[i].push_back(elem_below);
            esp.neighbour[elem_below].push_back(i);
         }
      }
   }
}

void MakeMesh(ESP & esp, libMesh::UnstructuredMesh& mesh){
   for (unsigned int i=0; i<esp.size; i++){
      mesh.add_point(esp.node[i], i);
   }

   BoundaryInfo& boundary_info = mesh.get_boundary_info();

   // Add sideset names to boundary info (Z axis out of the screen)
   boundary_info.sideset_name(0) = "back";
   boundary_info.sideset_name(1) = "bottom";
   boundary_info.sideset_name(2) = "right";
   boundary_info.sideset_name(3) = "top";
   boundary_info.sideset_name(4) = "left";
   boundary_info.sideset_name(5) = "front";

   // Add nodeset names to boundary info
   boundary_info.nodeset_name(0) = "back";
   boundary_info.nodeset_name(1) = "bottom";
   boundary_info.nodeset_name(2) = "right";
   boundary_info.nodeset_name(3) = "top";
   boundary_info.nodeset_name(4) = "left";
   boundary_info.nodeset_name(5) = "front";

   //find all elements that are left bottom front corner of a HEX-element
   double THRESHOLD=0.0001;
   double step=esp.node[1](0)-esp.node[0](0)+THRESHOLD; // add some tolerance
   std::vector<unsigned int> span;
   esp.used.resize(esp.size);
   bool back, bottom, left, front, top, right;
   // first set up all elements:
   for(unsigned int k=0; k<esp.size; k++){
      // in other case this index will not be able to span an element.
      if (esp.neighbour[k].size()>=7){
         //now check for nodes with coordinates with smaller values:
         span.clear();
         for(unsigned int i=0; i<esp.neighbour[k].size(); i++){
            if (esp.node[k](0)>=esp.node[esp.neighbour[k][i]](0)
                  && esp.node[k](1)>=esp.node[esp.neighbour[k][i]](1)
                  && esp.node[k](2)>=esp.node[esp.neighbour[k][i]](2)){
               span.push_back(esp.neighbour[k][i]);
               if (span.size()==7)
                  continue;
            }
         }
         // the element itself is also part of the element!
         span.push_back(k);
         //sort the span-vector:
         std::sort (span.begin(), span.end());
         if (span.size()>=8){ // '>' should never be true!
            if (span.size()>8){
               out<<"   too many neighbours!!"<<std::endl;
            }
            back= bottom= left= front= top= right=true;
            Elem* elem=mesh.add_elem(new Hex8);
            for (unsigned int j=0; j<8; j++){
                // this is needed to fit the enumeration of nodes
                // in the Hex8-element as used in libmesh.
                if (j==2)
                   elem->set_node(3)=mesh.node_ptr(span[j]);
                else if (j==3)                           
                   elem->set_node(2)=mesh.node_ptr(span[j]);
                else if (j==6)                           
                   elem->set_node(7)=mesh.node_ptr(span[j]);
                else if (j==7)                           
                   elem->set_node(6)=mesh.node_ptr(span[j]);
               else
                  elem->set_node(j)=mesh.node_ptr(span[j]);
               esp.used[span[j]]=true;
               // now, check if this element has a surface to outer region:
               if(esp.neighbour[span[j]].size()>=26){
                  // this element has all possible neighbours
                  // and hence has no surface.
                  continue;
               }
               //
               for(unsigned int i=0; i<esp.neighbour[span[j]].size(); i++){
                  //check, if it is a nearest neighbour:
                  if ((esp.node[esp.neighbour[span[j]][i]]-esp.node[span[j]]).size_sq()<step*step){
                     if (esp.node[esp.neighbour[span[j]][i]](0)!=esp.node[k](0)){
                        if (esp.node[esp.neighbour[span[j]][i]](0)>esp.node[k](0)){
                           right=false;
                           continue;
                        }
                        else if (esp.node[esp.neighbour[span[j]][i]](0)<esp.node[k](0)-step){
                           left=false;
                           continue;
                        }
                     }
                     if (esp.node[esp.neighbour[span[j]][i]](1)!=esp.node[k](1)){
                        if (esp.node[esp.neighbour[span[j]][i]](1)>esp.node[k](1)){
                           back=false;
                           continue;
                        }
                        else if (esp.node[esp.neighbour[span[j]][i]](1)<esp.node[k](1)-step){
                           front=false;
                           continue;
                        }
                     }
                     if (esp.node[esp.neighbour[span[j]][i]](2)!=esp.node[k](2)){
                        if (esp.node[esp.neighbour[span[j]][i]](2)>esp.node[k](2))
                           top=false;
                        else if (esp.node[esp.neighbour[span[j]][i]](2)<esp.node[k](2)-step)
                           bottom=false;
                     }
                  }
               }
            }
            if (bottom)
               boundary_info.add_side(elem, 1, 1);
            if (top)
               boundary_info.add_side(elem, 3, 3);
            if (left)
               boundary_info.add_side(elem, 4, 4);
            if (right)
               boundary_info.add_side(elem, 2, 2);
            if (front)
               boundary_info.add_side(elem, 5, 5);
            if (back)
               boundary_info.add_side(elem, 0, 0);
         }
      }
   }
}

void GetPotential1(ESP& esp, EquationSystems& equation_systems){
   const MeshBase& mesh = equation_systems.get_mesh();

   // create an explicit system to load the solution into:
   ExplicitSystem & pot = equation_systems.get_system<ExplicitSystem> ("esp");
   
   MeshBase::const_node_iterator           nd = mesh.nodes_begin();
   const MeshBase::const_node_iterator nd_end = mesh.nodes_end();
   unsigned int dn=0;
   for (; nd != nd_end; ++nd)
     {

      if (dn >=pot.rhs->size()){
         out<<"error with the sizes"<<std::endl;
         break;
      }
      //this may work, if the elements are not renumbered already.
      pot.rhs->set(dn, esp.potential[dn]);
      pot.solution->set(dn, esp.potential[dn]);
      dn++;
     }
}

unsigned int minind(Point q_point, ESP& esp, unsigned int guessind){
   int mini=-20;
   Real dist=(q_point-esp.node[guessind]).norm_sq();
   //out<<" ._.  "<<dist<<std::endl;
   Real distTest;
   for(unsigned int i=0; i<esp.neighbour[guessind].size(); i++){
      distTest=(q_point - esp.node[esp.neighbour[guessind][i]]).norm_sq();
      //out<<" ...  "<<distTest<<std::endl;
      if (distTest< dist){
         dist=distTest;
         mini=esp.neighbour[guessind][i];
      }
   }
   //out<<" ._.  "<<dist<<std::endl;
   if (mini> -10)
      return mini;
   else
      return guessind;
}

unsigned int findIndex(Point q_point, unsigned int  guessind, ESP& esp){
   unsigned int oldguess=guessind+2;
   //out<<" ->   "<<guessind<<std::endl;
   //while (!q_point.absolute_fuzzy_equals(esp.node[guessind], 0.005)){ //--> deltax + deltay + deltaz < tolerance.
   int  plus=0;
   //while ((q_point - esp.node[guessind]).norm_sq() > 0.005){
   while ((q_point - esp.node[guessind]).norm_sq() > 0.03571*0.866025){ // 0.03571 is the square of distance between NN
      oldguess=guessind;
      guessind=minind(q_point, esp, guessind);
      //if I didn't make progress, try to start at an other point:
      if(guessind==oldguess){
         if (plus==0){
            guessind-=2; // just lets try this!
            plus=1;
         }
         else if (plus==1){
            guessind+=2; // just lets try this!
            plus=2;
         }
         else{
            // I have no solution to this problem.
            out<<"stuck in loop!   ";
            out<<guessind<<"  "<<q_point;
            out<<"   "<<(q_point - esp.node[guessind]).norm_sq()<<std::endl;
            break;
         }
      }
   }
   return guessind;
}
         
double interpolate(ESP& esp, unsigned int dn, Point qp){
   double pot=esp.potential[dn]/(qp - esp.node[dn]).norm_sq();
   double norm=1/(qp - esp.node[dn]).norm_sq();
 // // a) weighted mean value:
 //  for(unsigned int i=0; i<esp.neighbour[dn].size(); i++){
 //     //-> interpolate using some interpolation formula:
 //     pot+=esp.potential[esp.neighbour[dn][i]]/(qp - esp.node[esp.neighbour[dn][i]]).norm_sq();
 //     norm+=1/(qp - esp.node[esp.neighbour[dn][i]]).norm_sq();
 //  }
 // // b) linear interpolation
 // // c) 1/r interpolation!?
   return pot/norm;
}

void GetPotential2(ESP& esp, EquationSystems& equation_systems){
   const MeshBase& mesh = equation_systems.get_mesh();

   // create an explicit system to load the solution into:
   ExplicitSystem & pot = equation_systems.get_system<ExplicitSystem> ("esp");
   
   MeshBase::const_node_iterator           nd = mesh.nodes_begin();
   const MeshBase::const_node_iterator nd_end = mesh.nodes_end();
   DofMap& dof_map=pot.get_dof_map();

   // This vector will hold the degree of freedom indices for
   // the element.  These define where in the global system
   // the element degrees of freedom get mapped.
   std::vector<dof_id_type> dof_indices;
   
   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = pot.get_dof_map().variable_type(0);
      
   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(3, fe_type));  // here, try AutoPtr instead...
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(3, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (3, fe_type.default_quadrature_order());
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);
   
   unsigned int dn;
   double espot;

   MeshBase::const_element_iterator el= mesh.elements_begin();
   const MeshBase::const_element_iterator end_el =  mesh.elements_end();
      
   for ( ; el != end_el; ++el){
      const Elem* elem = *el;
      
      dof_map.dof_indices (elem, dof_indices);
      
      // unifyging finite and infinite elements
      FEBase * cfe = libmesh_nullptr;

      if (elem->infinite()){
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }
      
      const std::vector<Point>& q_point = cfe->get_xyz();
      cfe->reinit (elem);

      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){

         unsigned int n_sf = cfe->n_shape_functions();
         dn=findIndex(q_point[qp], dn+1, esp);
         espot=interpolate(esp, dn, q_point[qp]);
         for (unsigned int i=0; i<n_sf; i++){
            //dn=findIndex(q_point[qp], dof_indices[i], esp);
            //pot.rhs->set(dof_indices[i], esp.potential[dn]);
            //pot.solution->set(dof_indices[i], esp.potential[dn]);
            pot.rhs->set(dof_indices[i], espot);
            pot.solution->set(dof_indices[i], espot);
         }
      }
   }
}

void GetPotential3(ESP& esp, EquationSystems& equation_systems){
   const MeshBase& mesh = equation_systems.get_mesh();

   // create an explicit system to load the solution into:
   ExplicitSystem & pot = equation_systems.get_system<ExplicitSystem> ("esp");
   
   MeshBase::const_node_iterator           nd = mesh.nodes_begin();
   const MeshBase::const_node_iterator nd_end = mesh.nodes_end();
   DofMap& dofs=pot.get_dof_map();
   out<<mesh.n_nodes()<<"   "<<esp.size<<"    "<<pot.rhs->size()<<std::endl;
   out<<dofs.n_dofs()<<std::endl;
   MeshBase::const_element_iterator el= mesh.elements_begin();
   const MeshBase::const_element_iterator end_el =  mesh.elements_end();
   //bool same;
   std::vector<dof_id_type> id;
   for ( ; el != end_el; ++el){
      const Elem* elem = *el;
      dofs.dof_indices(elem, id);
      for (unsigned int i=0; i<id.size(); i++){
         //Node* node=elem->get_node(i);
         //same=node->relative_fuzzy_equals(esp.node[id[i]], 0.001);
         //if (same){
         pot.rhs->set(id[i], esp.potential[id[i]]);
         pot.solution->set(id[i], esp.potential[id[i]]);
         //}
         //else
         //   out<<"error finding the correct node"<<std::endl;
      }
   }
}

void writeESP(EquationSystems & equation_systems){
   const MeshBase& mesh = equation_systems.get_mesh();
   const unsigned int n_nodes = mesh.n_nodes();
   // print the head line. The number here will be smaller because not all
   // points contribute to the mesh...
   out<< n_nodes<<"      "<<1.000000<<std::endl;
   ExplicitSystem & pot = equation_systems.get_system<ExplicitSystem> ("esp");
     
   MeshBase::const_element_iterator el= mesh.elements_begin();
   const MeshBase::const_element_iterator end_el =  mesh.elements_end();
   for ( ; el != end_el; ++el){
      const Elem* elem = *el;
      for(unsigned int i=0; i<elem->n_nodes(); i++){
         Node* node=elem->get_node(i);
         out<< node->id()<<"     ";
         out<<node->operator() (0)<<"  ";
         out<<node->operator() (1)<<"  ";
         out<<node->operator() (2)<<"  ";
         //out<<std::endl;
         out<< pot.rhs->el(node->id())<<std::endl;
      }
   }
}

void writeESP2(EquationSystems & equation_systems){
   const MeshBase& mesh = equation_systems.get_mesh();

   // create an explicit system to load the solution into:
   ExplicitSystem & pot = equation_systems.get_system<ExplicitSystem> ("esp");
   
   MeshBase::const_node_iterator           nd = mesh.nodes_begin();
   const MeshBase::const_node_iterator nd_end = mesh.nodes_end();
   DofMap& dof_map=pot.get_dof_map();

   // This vector will hold the degree of freedom indices for
   // the element.  These define where in the global system
   // the element degrees of freedom get mapped.
   std::vector<dof_id_type> dof_indices;
   
   // Get a constant reference to the Finite Element type
   // for the first (and only) variable in the system.
   FEType fe_type = pot.get_dof_map().variable_type(0);
      
   // Build a Finite Element object of the specified type.  Since the
   // \p FEBase::build() member dynamically creates memory we will
   // store the object as an \p UniquePtr<FEBase>.  This can be thought
   // of as a pointer that will clean up after itself.
   UniquePtr<FEBase> fe (FEBase::build(3, fe_type));  // here, try AutoPtr instead...
   UniquePtr<FEBase> inf_fe (FEBase::build_InfFE(3, fe_type));
   
   // A  Gauss quadrature rule for numerical integration.
   // Use the default quadrature order.
   QGauss qrule (3, fe_type.default_quadrature_order());
      
   // Tell the finite element object to use our quadrature rule.
   fe->attach_quadrature_rule (&qrule);
   inf_fe->attach_quadrature_rule (&qrule);

   MeshBase::const_element_iterator el= mesh.elements_begin();
   const MeshBase::const_element_iterator end_el =  mesh.elements_end();
      
   for ( ; el != end_el; ++el){
      const Elem* elem = *el;
      
      dof_map.dof_indices (elem, dof_indices);
      
      // unifyging finite and infinite elements
      FEBase * cfe = libmesh_nullptr;

      if (elem->infinite()){
         cfe = inf_fe.get();
      }
      else{
        cfe = fe.get();
      }

      const std::vector<Point>& q_point = cfe->get_xyz();
      cfe->reinit (elem);

      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++){

         unsigned int n_sf = cfe->n_shape_functions();
         for (unsigned int i=0; i<n_sf; i++){
            out<<dof_indices[i]<<"   ";
            out<<q_point[qp]<<"   ";
            out<<pot.rhs->operator()(dof_indices[i])<<std::endl;
         }
      }
   }
}

EquationSystems & InsertPot(std::string potfile, Mesh& pot_mesh, EquationSystems & equation_systems){
   struct ESP esp;
   Read(esp, potfile);
   FindNeighbours(esp);
   MakeMesh(esp, pot_mesh);

   // find neighbours and other things:
   pot_mesh.prepare_for_use(true, false);
   // create an explicit system to load the solution into:
   ExplicitSystem& pot = equation_systems.add_system<ExplicitSystem> ("esp");

   // Create an FEType describing the approximation
   // characteristics of the InfFE object.  Note that
   // the constructor automatically defaults to some
   // sensible values.  But use \p FIRST order
   // approximation.
   FEType fe_type(FIRST);
   pot.add_variable("potential", fe_type);
   // Initialize the data structures for the equation system.
   equation_systems.init();
   equation_systems.print_info();
  
   GetPotential2(esp,equation_systems);
   ExodusII_IO (pot_mesh).write_equation_systems("potential.e", equation_systems);
   // for checking, if everything worked,
   // try to reconstruct K_esp.grid, just using the pot_mesh:
   //writeESP(equation_systems);
   writeESP2(equation_systems);
   return equation_systems;
}
