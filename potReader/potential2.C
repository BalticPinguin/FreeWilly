#include <vector>
#include <fstream>
#include "libmesh/point.h" // to have access to "point"-object
#include "libmesh/mesh.h"  // for 'mesh'
#include "libmesh/elem.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex27.h"  // this is only dummy and will be removed later.
#include "libmesh/nanoflann.hpp"  // for FindeNeighors()
#include <cmath> // needed for abs().

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
   // idea of the function:

   unsigned int y=0, z=0, y_new=0, z_new=0, elem_below;
   int overhang_y, overhang_z;
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
         overhang_y=0; if (esp.node[i](0)>esp.node[y+i-y_new](0)+THRESHOLD){ // the last line started more left of it: while(esp.node[i](0)>esp.node[y+overhang_y+i-y_new](0)+THRESHOLD) overhang_y++; overhang_y=-overhang_y; } else if (esp.node[i](0)<esp.node[y+i-y_new](0)-THRESHOLD){ // the last line started more right of it: while (esp.node[i+overhang_y](0)>esp.node[y+i-y_new](0)-THRESHOLD)
               overhang_y++;
         }
         // y is the start-value of the last row.
         y=y_new;
         // y_new is the start-value of this row.
         y_new=i;
      }
      if (esp.node[i](2)!=esp.node[i-1](2)){
         overhang_z=0;
         // adjust myself with respect to x
         // alignment wrt. to y is assumed to exist, since 
         // otherwise there would be an empty line in y...
         if (esp.node[i](0)>esp.node[z+(i-z_new)](0)){
            // the last line started more left of it:
            while (esp.node[i](0)>esp.node[z+overhang_z+(i-z_new)](0))
               overhang_z++;
            overhang_z=-overhang_z;
         }
         else if (esp.node[i](0)<esp.node[z+i-z_new](0)){
            // the last line started more right of it:
            while (esp.node[i+overhang_z](0)>esp.node[z+i-z_new](0))
               overhang_z++;
         }
         // z is the start-value of the last row.
         z=z_new;
         // z_new is the start-value of this row.
         z_new=i;
      }

      if (y_new>0){
         // there might be a whole! it might be in this or in previous layer...
         if (std::abs(esp.node[i](0)-esp.node[y+overhang_y+(i-y_new)](0))>step
                       || std::abs(esp.node[i](1)-esp.node[y+overhang_y+i-y_new+1](1))>step){
            // so adjust overhang_y accordingly and move on!
            if (esp.node[i](0)>esp.node[y+i-y_new](0)){
               // the last line started more left of it:
               while (esp.node[i](0)>esp.node[y+overhang_y+i-y_new](0))
                  overhang_y++;
            }
            else if (esp.node[i](0)<esp.node[y+i-y_new](0)){
               // the last line started more right of it:
               while (esp.node[i+overhang_y](0)>esp.node[y+i-y_new](0))
                  overhang_y++;
               overhang_y=-overhang_y;
            }
         }
         if (std::abs(esp.node[i](0)-esp.node[y+overhang_y+(i-y_new)](0))<=step
                && std::abs(esp.node[i](1)-esp.node[y+overhang_y+(i-y_new)](1))<=step){
            // the stepping in y-direction is asured to be exactly 'step'.
            esp.neighbour[i].push_back(y+overhang_y+i-y_new);
            esp.neighbour[y+overhang_y+i-y_new].push_back(i);
         }
         // get the diagonal elements as well
         if (i-y_new>0 && std::abs(esp.node[i](0)-esp.node[y+overhang_y+i-y_new-1](0))<=step
                       && std::abs(esp.node[i](1)-esp.node[y+overhang_y+i-y_new-1](1))<=step){
            // the correct step-size in y- and z-direction is assured.
            esp.neighbour[i].push_back(y+overhang_y-1+i-y_new);
            esp.neighbour[y+overhang_y-1+i-y_new].push_back(i);
         }
         if (std::abs(esp.node[i](0)-esp.node[y+overhang_y+i-y_new+1](0))<=step
                       && std::abs(esp.node[i](1)-esp.node[y+overhang_y+i-y_new+1](1))<=step){
            // the correct step-size in y- and z-direction is assured.
            esp.neighbour[i].push_back(y+overhang_y+1+i-y_new);
            esp.neighbour[y+overhang_y+1+i-y_new].push_back(i);
         }
      }
      if (z_new>0){
         // even in 3D, only the variable that changes fastest needs to be considered explicitly.
         if (std::abs(esp.node[i](0)-esp.node[z+overhang_z+i-z_new](0))>step
               || std::abs(esp.node[i](1)-esp.node[z+overhang_z+i-z_new](1))>step){
            // so adjust overhang_y accordingly and move on!
            if (esp.node[i](0)>esp.node[z+i-z_new](0)){
               // the last line started more left of it:
               out<<"adjust overhang_z+"<<std::endl;
               while (esp.node[i](0)>esp.node[z+overhang_z+i-z_new](0))
                  overhang_z++;
               overhang_z=-overhang_z;
            }
            else if (esp.node[i](0)<esp.node[z+i-z_new](0)){
               // the last line started more right of it:
               out<<"adjust overhang_z-"<<std::endl;
               while (esp.node[i+overhang_z](0)>esp.node[z+i-z_new](0))
                  overhang_z++;
            }
         }
         elem_below=z+overhang_z+i-z_new;
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
               }
            }
            // the correct step-size in y- and z-direction is assured.
            esp.neighbour[i].push_back(elem_below);
            esp.neighbour[elem_below].push_back(i);
         }
      }
   }
}

void MakeMesh(ESP & esp, libMesh::UnstructuredMesh& mesh, const ElemType Eltype=HEX8){
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

   switch (Eltype){
      case HEX8:{
         //find all elements that are left bottom front corner of a HEX-element
         double THRESHOLD=0.00001;
         double step=esp.node[1](0)-esp.node[0](0)+THRESHOLD; // add some tolerance
         std::vector<unsigned int> span;
         bool back, bottom, left, front, top, right;
         // first set up all elements:
         for(unsigned int k=0; k<esp.size; k++){
            // in other case this index will not be able to span an element.
            if (esp.neighbour[k].size()>8){
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
         break;
      }
      default:
         libmesh_error_msg("only HEX8 is available now.");
   }
}

int main (int argc, char** argv){
   struct ESP esp;
   Read(esp, "K_esp.grid");
   //for(unsigned int i=0; i<esp.size; i++){
   //   out<<esp.node[i](0)<<"  ";
   //   out<<esp.node[i](1)<<"  ";
   //   out<<esp.node[i](2)<<"  "<<std::endl;
   //}
   FindNeighbours(esp);
   //for(unsigned int i=0; i<esp.size; i++){
   //   out<<esp.node[i](0)<<"  ";
   //   out<<esp.node[i](1)<<"  ";
   //   out<<esp.node[i](2)<<"  "<<std::endl;
   //   out<<"   neighbours:\n    ";
   //   for (unsigned int j=0;j<esp.neighbour[i].size(); j++){
   //      out<<esp.neighbour[i][j]<<"  ";
   //   }
   //   out<<std::endl;
   //}

   LibMeshInit init (argc, argv);
   Mesh mesh(init.comm(), 3);
   MakeMesh(esp, mesh);
   mesh.write("foobar.e");
   return 0;
}

// b potential2.C:59
