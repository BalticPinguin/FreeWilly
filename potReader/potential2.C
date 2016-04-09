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

using namespace libMesh;

//function to read the data from file and feed a mesh with respective nodes. After this, a .rhs is taken and the values are added:

struct ESP{
   // these vectors store the points and potential (at the points) as given.
   std::vector<Point> node;
   std::vector<double> potential;
   // k-th node opens an element, if there are 7 elements with distance <=3^1/2*min_dist
   // whose components are >= the components of k-th node.
   std::vector<std::vector<unsigned int>> neighbour; // fill this with nanoflann!
   std::vector<std::vector<unsigned int>> setup_elem; 
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
      esp.potential[i]=v;
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
   double THRESHOLD=0.00001;
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
   unsigned int node_id = 0;

   for (unsigned int i=0; i<esp.size; i++){
      mesh.add_point(esp.node[i], node_id++);
   }
   // how does the system know which element belongs to which node?

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
   double THRESHOLD=0.00001;
   double step=esp.node[1](0)-esp.node[0](0)+THRESHOLD; // add some tolerance
   std::vector<unsigned int> span;
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

void GetPotential(ESP& esp, EquationSystems& equation_systems){
   const MeshBase& mesh = equation_systems.get_mesh();

   // create an explicit system to load the solution into:
   ExplicitSystem & pot = equation_systems.get_system<ExplicitSystem> ("esp");

   {
      MeshBase::const_node_iterator           nd = mesh.local_nodes_begin();
      const MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
      unsigned int dn=0;
      for (; nd != nd_end; ++nd)
        {
        // this may work, if the elements are not renumbered already.
         pot.solution->set (dn, esp.potential[dn]);
         dn++;
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
  // MeshBase::const_node_iterator           nd = mesh.nodes_begin();
  // const MeshBase::const_node_iterator nd_end = mesh.nodes_end();
  // for (; nd != nd_end; ++nd){
  // }

}

int main (int argc, char** argv){
   struct ESP esp;
   Read(esp, "K_esp.grid");
   FindNeighbours(esp);
   LibMeshInit init (argc, argv);
   Mesh mesh(init.comm(), 3);
   MakeMesh(esp, mesh);

   // find neighbours, renumber nodes and other things:
   // do NOT renumber!
   mesh.prepare_for_use(true, false);
   EquationSystems equation_systems (mesh);
   // create an explicit system to load the solution into:
   //ExplicitSystem & pot = equation_systems.add_system<ExplicitSystem> ("esp");
   ExplicitSystem& pot = equation_systems.add_system<ExplicitSystem> ("esp");

   // Create an FEType describing the approximation
   // characteristics of the InfFE object.  Note that
   // the constructor automatically defaults to some
   // sensible values.  But use \p FIRST order
   // approximation.
   FEType fe_type(FIRST);
   pot.add_variable("p", fe_type);
   // Initialize the data structures for the equation system.
   equation_systems.init();
   equation_systems.print_info();
  
   GetPotential(esp,equation_systems);
   mesh.write("foobar.e");
   ExodusII_IO (mesh).write_equation_systems( "system.e", equation_systems);
   // for checking, if everything worked,
   // try to reconstruct K_esp.grid, just using the mesh:
   writeESP(equation_systems);
   return 0;
}

// b potential2.C:59
