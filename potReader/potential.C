#include <vector>
#include <fstream>
#include "libmesh/point.h" // to have access to "point"-object
#include "libmesh/mesh.h"  // for 'mesh'
#include "libmesh/elem.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex27.h"  // this is only dummy and will be removed later.
#include "libmesh/nanoflann.hpp"  // for FindeNeighors()

using namespace libMesh;

//function to read the data from file and feed a mesh with respective nodes. After this, a .rhs is taken and the values are added:

struct ESP{
   std::vector<Point> node;
   // k-th node opens an element, if there are 7 elements with distance <=3^1/2*min_dist
   // whose components are >= the components of k-th node.
   std::vector<std::vector<unsigned int>> neighbour; // fill this with nanoflann!
   std::vector<std::vector<unsigned int>> setup_elem; 
   std::vector<double> potential;
   unsigned int size;
   double dist;

   template <class BBOX>
   bool kdtree_get_bbox(BBOX&) const {return false;}
   
   inline unsigned int kdtree_get_point_count() const { return this->size;}

   inline double kdtree_get_pt(const unsigned int  idx, int dim) const{
         if (dim==0) return node[idx](0);
         else if (dim==1) return node[idx](1);
         else return node[idx](2);
    }

   // return the distance between a point and a second point with index idx.
   inline double kdtree_distance(const Point point, const unsigned int idx) const {
      double x,y,z;
      x=point(0)-this->node[idx](0);
      y=point(1)-this->node[idx](1);
      z=point(2)-this->node[idx](2);
      return x*x+y*y+z*z;
   }
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
   esp.dist=esp.node[1](0)-esp.node[0](0);
   esp_file.close();
}

void FindNeighbors(ESP & esp){
   const unsigned int search_radius = esp.dist*1.733; //sqrt(3)
   std::vector<std::pair<long unsigned int, double> > ret_matches;

   nanoflann::SearchParams params;
   params.sorted = false;
   esp.neighbour.resize(esp.size);
   esp.setup_elem.resize(esp.size);

   // construct a kd-tree index:
   typedef nanoflann::KDTreeSingleIndexAdaptor<
                      nanoflann::L2_Simple_Adaptor<Real, ESP > ,
                      ESP, 3/* dim */ > esp_kd_tree_t;
   esp_kd_tree_t index(3 /*dim*/, esp , nanoflann::KDTreeSingleIndexAdaptorParams(27 /* max leaf */ ) );
   index.buildIndex();

   for(unsigned int i=0; i<esp.size; i++){
      //search for points that share a common cube with the i-th point:
      Real* query_pt[3];
      *query_pt[0]=esp.node[i](0);
      *query_pt[1]=esp.node[i](1);
      *query_pt[2]=esp.node[i](2);
      const size_t nMatches = index.radiusSearch(query_pt[0], search_radius, ret_matches, params);
      for (unsigned int j=0; j<nMatches; j++){
         esp.neighbour[i].push_back(ret_matches[j].first);
         if ((      esp.node[ret_matches[j].first](0)>=esp.node[i](0))
                && (esp.node[ret_matches[j].first](1)>=esp.node[i](1))
                && (esp.node[ret_matches[j].first](2)>=esp.node[i](2))){
            esp.setup_elem[i].push_back(ret_matches[j].first);
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
      case HEX27:{
         // the '-1' below is wrong, but maybe leads to not having the error!
         for(unsigned int k=0; k<esp.size; k++){
            Elem* elem=mesh.add_elem(new Hex8);
            for (unsigned int j=0; j<8; j++){
               out<<k<<"  "<<j<<"  "<<k*8+j<<std::endl;
               //elem->set_node(j)=mesh.node_ptr(complicated function of);
               elem->set_node(j)=mesh.node_ptr(k*8+j);
               // do I need to enumerate tham as above??
            }
            if (k==0){
               boundary_info.add_side(elem, 0, 0);
               boundary_info.add_side(elem, 1, 1);
               boundary_info.add_side(elem, 4, 4);
            }
            else if (k==esp.size-1){
               boundary_info.add_side(elem, 5, 5);
               boundary_info.add_side(elem, 3, 3);
               boundary_info.add_side(elem, 2, 2);
            }
            else{
               //I am at a boarder, when y and/or z change from this to the previous element.
               if (esp.node[k](1)>esp.node[k+1](1))
                  boundary_info.add_side(elem, 2, 2);
               if (esp.node[k](1)>esp.node[k+1](1))
                  boundary_info.add_side(elem, 3, 3);
               //if (esp[k].node(2)>esp[k+1].node(2)) in z-direction the values always grow.
               //   boundary_info.add_side(elem, 5, 5);
               if( esp.node[k](0)<esp.node[k-1](0))
                  boundary_info.add_side(elem, 4, 4);
               if( esp.node[k](1)<esp.node[k-1](1))
                  boundary_info.add_side(elem, 1, 1);
               if( esp.node[k](2)<esp.node[k-1](2))
                  boundary_info.add_side(elem, 0, 0);
            }
         }
      }
      case HEX8:{
         //find all elements that are left bottom front corner of a HEX-element
         FindNeighbors( esp);
         bool x_min, x_max, y_min, y_max, z_min, z_max;
         // first set up all elements:
         for(unsigned int k=0; k<esp.size; k++){
            if (esp.setup_elem[k].size()==7){
               Elem* elem=mesh.add_elem(new Hex8);
               for (unsigned int j=0; j<8; j++){
                  elem->set_node(j)=mesh.node_ptr(esp.neighbour[k][j]);
               }
               if (esp.neighbour[k].size()<27){
                  x_min=false;
                  y_min=false;
                  z_min=false;
                  x_max=false;
                  y_max=false;
                  z_max=false;
                  for(unsigned int j=0; j<esp.neighbour[k].size(); j++){
                     if (esp.node[esp.neighbour[k][j]](0)<esp.node[k](0)){
                        x_min=true;
                     }
                     if (esp.node[esp.neighbour[k][j]](1)<esp.node[k](1)){
                        y_min=true;
                     }
                     if (esp.node[esp.neighbour[k][j]](2)<esp.node[k](2)){
                        z_min=true;
                     }
                     if (esp.node[esp.neighbour[k][j]](0)>esp.node[k](0)){
                        x_max=true;
                     }
                     if (esp.node[esp.neighbour[k][j]](1)>esp.node[k](1)){
                        y_max=true;
                     }
                     if (esp.node[esp.neighbour[k][j]](2)>esp.node[k](2)){
                        z_max=true;
                     }
                  }
                  // no neighbour has larger y-value
                  if (!y_max)
                     boundary_info.add_side(elem, 0, 0);
                  // no neighbour has smaller z-value
                  if (!z_min)
                     boundary_info.add_side(elem, 1, 1);
                  // no neighbour has larger x-value
                  if (!x_max)
                     boundary_info.add_side(elem, 2, 2);
                  // no neighbour has larger z-value
                  if (!z_max)
                     boundary_info.add_side(elem, 3, 3);
                  // no neighbour has smaller x-value
                  if (!x_min)
                     boundary_info.add_side(elem, 4, 4);
                  // no neighbour has smaller y-value
                  if (!y_min)
                     boundary_info.add_side(elem, 5, 5);
            
               }
            }
         }
      }
      default:
         libmesh_error_msg("only HEX8 is available now.");
   }
}

int main (int argc, char** argv){
   struct ESP esp;
   Read(esp, "K_esp.grid");
   LibMeshInit init (argc, argv);
   Mesh mesh(init.comm(), 3);
   MakeMesh(esp, mesh);
   mesh.write("foobar.e");
   return 0;
}

