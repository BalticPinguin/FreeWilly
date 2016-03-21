#include <vector>
#include <fstream>
#include "libmesh/point.h" // to have access to "point"-object
#include "libmesh/mesh.h"  // for 'mesh'
#include "libmesh/elem.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/face_quad4.h"  // this is only dummy and will be removed later.

using namespace libMesh;

//function to read the data from file and feed a mesh with respective nodes. After this, a .rhs is taken and the values are added:

struct ESP{
   std::vector<Point> node;
   std::vector<double> potential;
   unsigned int size;
};

ESP Read(std::string input_file){
   double noonecares;
   double x,y,z,v;
   std::ifstream esp_file;
   struct ESP esp;

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
   return esp;
}

unsigned int idx(const ElemType type, const unsigned int nx,
                const unsigned int ny, const unsigned int i,
                const unsigned int j, const unsigned int k){
   // this function is just a shortened copy of a function of same name in namespace 
   // MeshTools::Generation::Private.
   switch (type){
      case HEX8:{
         return i+(nx+1)*(j+k*(ny+1));
      }
      default:
      libmesh_error_msg("ERROR: Unrecognized element type.");
   }
}

void MakeMesh(const ESP esp, libMesh::UnstructuredMesh& mesh, const ElemType Eltype=HEX8){
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
      case QUAD4:{ // this is hex8 actually as well.
         const unsigned int nx=2,ny=2,nz=2;
         for(unsigned int k=0; k<nz; k++){
            for(unsigned int j=0; j<ny; j++){
               for(unsigned int i=0; i<nx; i++){
                  Elem* elem=mesh.add_elem(new Quad4);
                  elem->set_node(0)=mesh.node_ptr(idx(Eltype,nx,ny,i  ,j  ,k  ) );
                  elem->set_node(1)=mesh.node_ptr(idx(Eltype,nx,ny,i+1,j  ,k  ) );
                  elem->set_node(2)=mesh.node_ptr(idx(Eltype,nx,ny,i+1,j+1,k  ) );
                  elem->set_node(3)=mesh.node_ptr(idx(Eltype,nx,ny,i  ,j+1,k  ) );
                  elem->set_node(4)=mesh.node_ptr(idx(Eltype,nx,ny,i  ,j  ,k+1) );
                  elem->set_node(5)=mesh.node_ptr(idx(Eltype,nx,ny,i+1,j  ,k+1) );
                  elem->set_node(6)=mesh.node_ptr(idx(Eltype,nx,ny,i+1,j+1,k+1) );
                  elem->set_node(7)=mesh.node_ptr(idx(Eltype,nx,ny,i  ,j+1,k+1) );
                  if (k == 0)
                    boundary_info.add_side(elem, 0, 0);

                  if (k == (nz-1))
                    boundary_info.add_side(elem, 5, 5);

                  if (j == 0)
                    boundary_info.add_side(elem, 1, 1);

                  if (j == (ny-1))
                    boundary_info.add_side(elem, 3, 3);

                  if (i == 0)
                    boundary_info.add_side(elem, 4, 4);

                  if (i == (nx-1))
                    boundary_info.add_side(elem, 2, 2);
               }
            }
         }
         break;
      }
      case HEX8:{
         for(unsigned int k=0; k<esp.size; k++){
            Elem* elem=mesh.add_elem(new Hex8);
            out<<k<<std::endl;
            for (unsigned int j=0; j<8; j++){
               elem->set_node(j)=mesh.node_ptr(k*8+j);
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
      default:
         libmesh_error_msg("only HEX8 is available now.");
   }
}

int main (int argc, char** argv){
   struct ESP esp;
   esp=Read("K_esp.grid");
   LibMeshInit init (argc, argv);
   Mesh mesh(init.comm(), 3);
   MakeMesh(esp, mesh);
   mesh.write("foobar.e");
   return 0;
}

