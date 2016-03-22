#include <vector>
#include <fstream>
#include "libmesh/point.h" // to have access to "point"-object

using namespace libMesh;

//function to read the data from file and feed a mesh with respective nodes. After this, a .rhs is taken and the values are added:

std::vector<Point> Read(std::string input_file){
   unsigned int num_points, noonecares;
   std::ifstream esp_file;
   std::vector<Point> nodes;
   std::vector<double> ESP;

   esp_file.open(input_file.c_str(),ios::in);
   if (!esp_file.open()){
      //give an error message.
   }

   esp_file>>num_points>> noonecares;
   nodes.resize(num_points);
   ESP.resize(num_points);
   for(std::string line, int i=0; getline(esp_file, line); i++){
      esp_file >> x>>y>>z>>v;
      nodes[i]=Point(x,y,z);
      ESP[i]=v;
   }
   esp_file.close();
   return nodes;

}

void main(void){
   std::vector<Point> foo;
   foo=Read();
   out<<foo<<std::endl;
}
