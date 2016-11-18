// Local header files
#include "FreeWilly.h"


// Bring in everything from the libMesh namespace using namespace libMesh;
using namespace libMesh;

unsigned int closest(std::vector<Node> geom, Point pt){
   unsigned int molsize=geom.size();
   //most outer loop: nuclear sites
   unsigned int closest=0;
   Point cl= geom[0]-pt;
   Point ac;
   for (unsigned int i=1; i<molsize; i++ ){
      ac= geom[i]-pt;
      if ( ac.size()<cl.size() ){
         closest = i;
         cl = ac;
      }
   }
   return closest;
}

void tetrahedralise_sphere(UnstructuredMesh& mesh, std::vector<Node> geometry, std::string creator, Real r, std::string scheme, Real p, Real VolConst, Real L, unsigned int N){
   #ifdef LIBMESH_HAVE_TETGEN
   
   add_sphere_convex_hull_to_mesh(mesh, r, scheme, p, geometry, creator, L, N);
   
   // 3.) Update neighbor information so that TetGen can verify there is a convex hull.
   mesh.find_neighbors();
   
   // 5.) Set parameters and tetrahedralize the domain
   
   // 0 means "use TetGen default value"
   Real quality_constraint = 1.5;
   //Real quality_constraint = 420; // practically no constraint
   
   // The volume constraint determines the max-allowed tetrahedral
   // volume in the Mesh.  TetGen will split cells which are larger than
   // this size
   Real volume_constraint = VolConst; 
   
   // Construct the Delaunay tetrahedralization
   TetGenMeshInterface t(mesh);
   
   t.pointset_convexhull();
   mesh.find_neighbors();
   t.set_switches("V");
   t.triangulate_conformingDelaunayMesh(quality_constraint, volume_constraint);
   
   
   // Find neighbors, etc in preparation for writing out the Mesh
   mesh.prepare_for_use();
   
   // Finally, write out the result --> is done in main already.
   //mesh.write("sphere_3D.e");
   #else
   // Avoid compiler warnings
   libmesh_ignore(mesh.comm);
   #endif // LIBMESH_HAVE_TETGEN
}

std::pair<int, Real> chose_scheme(int circle, int N, Real L, Real p, Real r_max, std::string scheme){
   Real scale;
   int pts_circle;
   if (scheme=="son"){
      scale=L*circle/(N-circle+L*N/r_max);
      pts_circle=(int)(12.5/ (1.- (circle-1.)/circle*(N-circle+ L*N/r_max)/(N-circle+1.+L*N/r_max) )
                           / (1.- (circle-1.)/circle*(N-circle+ L*N/r_max)/(N-circle+1.+L*N/r_max)));
   }
   else if( scheme=="tm"){
      Real power=pow((Real)N/circle,p);
      scale=L*circle/( power* (L*N/r_max-1.)+1. );
      pts_circle=(int)(12.5/ (1.- (circle-1.)/circle* (power*(N*L/r_max-1.)+1.)
                                                   /(pow(N/(circle-1.),p)*(N*L/r_max-1.)+1.) )
                           / (1.- (circle-1.)/circle* (power*(N*L/r_max-1.)+1.)
                                                   /(pow(N/(circle-1.),p)*(N*L/r_max-1.)+1.)));
   }
   else if( scheme=="tm_300"){
      Real power=pow((Real)N/circle,p);
      scale=L*circle/(power* (L*N/r_max-1.)+1.);
      pts_circle=std::min((int)(12.5/ (1.- (circle-1.)/circle* (power*(N*L/r_max-1.)+1.)
                                                   /(pow(N/(circle-1.),p)*(N*L/r_max-1.)+1.) )
                           / (1.- (circle-1.)/circle* (power*(N*L/r_max-1.)+1.)
                                                   /(pow(N/(circle-1.),p)*(N*L/r_max-1.)+1.)))
                     ,300);
   }
   else if (scheme == "const"){
      scale=pow(L,circle)*r_max/pow(L,N);
      pts_circle = std::min((int)(12.5*L*L/((1-L)*(1-L))), 300);
   }
   else if (scheme == "const_tm"){
      scale=L*circle/( pow((Real)N/circle,p)* (L*N/r_max-1.)+1. );
      pts_circle = std::min((int)(12.5*L*L/((1-L)*(1-L))), 300);
   }
   else if (scheme == "quadr" ){
      scale=circle/N*r_max;
      pts_circle=(int)(12.5*(circle*r_max*circle*r_max/(N*N)));
   }
   else if (scheme == "sqrt_tm" ){
      scale=L*circle/( pow((Real)N/circle,p)* (L*N/r_max-1.)+1. );
      pts_circle = (int)(40*sqrt(scale)); 
   } 

   return std::pair<int, Real>(pts_circle,scale);
}

std::vector<Point> fibonacci(unsigned int points_on_sphere){
   std::vector<Point> points(points_on_sphere);
   double x,y,z,r,phi;
   double pi=3.1415926;
   for (unsigned int i=0; i<points_on_sphere; i++){
      y=2.*i/points_on_sphere-1.;
      r=sqrt(1.-y*y);
      phi=i*pi*(3.-sqrt(5.));
      z=r*sin(phi);
      x=r*cos(phi);
      points[i]=Node(x,y,z);
      //out<<points[i]<<std::endl;
   }
   //out<<std::endl<<std::endl;
   return points;
}

std::vector<Point> archimedes(unsigned int points_on_sphere){
   int n=floor(sqrt(points_on_sphere));
   std::vector<Point> points(n*n);
   double theta, phi;
   for(int i=0; i<n; i++){
      theta=i*6.2831852/n;
      for(int j=1; j<=n; j++){
         phi=acos(2.*j/n -1.);
         points[i*n+j-1]=Node(cos(theta)*sin(phi),
                   sin(theta)*sin(phi), 
                   2.*j/(n+1.)-1.);
      //out<<points[i*n+j-1]<<std::endl;
      }
   }
   return points;
}

std::vector<Point> spiral(unsigned int points_on_sphere){
   std::vector<Point> points(points_on_sphere);
   double pi=3.1415926;
   double dl=pi*(3.-sqrt(5.));
   double dz=2./(points_on_sphere +1.);
   double z, l, r;
   for(unsigned int i=1; i<=points_on_sphere; i++){
      z=i*dz-1.;
      l=i*dl;
      r=sqrt(1.-z*z);
      points[i-1]=Node(r*cos(l),r*sin(l),z);
   }
   return points;
}

void add_sphere_convex_hull_to_mesh(MeshBase& mesh, libMesh::Real r_max, std::string scheme, Real p, std::vector<Node> geometry, std::string creator, const Real L, const unsigned int N){
   #ifdef LIBMESH_HAVE_TETGEN
   std::vector<Point> point;
   double scale;
   // For each point in the map, insert it into the input mesh and copy it to all nuclear sites.
   // keep track of the ID assigned.
   unsigned int molsize=geometry.size();
   unsigned int pts_circle;

   // add one node each at the nuclear position.
   for(unsigned int i=0; i<molsize; i++){
      // add one point each at the nuclear position:
      mesh.add_point( geometry[i]);
   }

   //now, add further points to the mesh: The nodes are located on spheres (circle)
   // around the molecules; each sphere has a different radius (scale) in (0,r_max].

   std::pair<int, Real> ball;
   //outermost loop: over different spheres
   for(unsigned int circle=1; circle<=N; circle++){
      // chose the number of points on the sphere and the radius according to the
      // scheme requested by the user:
      ball=chose_scheme(circle, N, L, p, r_max, scheme);
      pts_circle=ball.first;
      scale=ball.second;

      libmesh_assert(pts_circle>0);
      libmesh_assert(scale>0.);

      // out<<" sphere nr "<< circle; // out<<"  "<<scale;
      // out<<"  "<<pts_circle<<std::endl;
      if (creator=="fibonacci")
         point=fibonacci(pts_circle);
      else if (creator=="archimedes")
         point=archimedes(pts_circle);
      else if (creator=="spiral")
         point=spiral(pts_circle);
      else if (creator=="Sdesign"){
         int rule, num_pts;
         // pts_circle approx num_pts, get resp. rule.
         rule = design_size(pts_circle, &num_pts);
         double** x;
         x= new double*[3];
         x[0]=new double[num_pts];
         x[1]=new double[num_pts];
         x[2]=new double[num_pts];
         design_points(rule, num_pts, x);
         point.resize(num_pts);

         out<<"order: "<<num_pts<<"  ";
         out<<pts_circle<<std::endl;
         for(int i=0; i<num_pts; i++){
            point[i]=Point(x[0][i],x[1][i],x[2][i]);
         }
         delete [] x[0];
         delete [] x[1];
         delete [] x[2];
         delete[] x;
      }
      else{
         // in these schemes, no 'rule' is needed.
         int num_pts;
         //set the order:
         if (creator=="lebedev"){
            avail_pts(pts_circle, &num_pts);
         }
         // these functions are disabled until I get them somehow working.
         //else if (creator=="geodesic4"){
         //   rule=point_size(4, pts_circle, rule);
         //}
         //else if (creator=="geodesic6"){
         //   rule=point_size(6, pts_circle, rule);
         //}
         else{ // some Womersley
            Wom_precision_table(pts_circle, &num_pts);
         }
         // not all orders are implemented.
         out<<"order: "<<num_pts<<"  ";
         out<<pts_circle<<std::endl;
         double *x;
         double *y;
         double *z;
         double *w;
         x= new double[num_pts];
         y= new double[num_pts];
         z= new double[num_pts];
         w= new double[num_pts];
         if (creator=="lebedev")
            ld_by_order (num_pts, x, y, z, w);
         else if (creator=="womEV")
            Wom_points (1, num_pts, x, y, z, w);
         else if (creator=="womMD")
            Wom_points (2, num_pts, x, y, z, w);
         else if (creator=="womMN")
            Wom_points (3, num_pts, x, y, z, w);
         else if (creator=="fliME") 
            Wom_points (4, num_pts, x, y, z, w);
         //else if (creator=="geodesic4")
         //   gen_grid(x, y, z, rule, 4);
         //else if (creator=="geodesic6")
         //   gen_grid(x, y, z, rule, 6);
         else
            libmesh_error_msg("no valid creator specified.\n");
         // this is not needed at all.
         delete [] w; 

         // now, convert the arrays to a Point-vector.
         point.resize(num_pts);
         for(int i=0; i<num_pts; i++){
            point[i]=Point(x[i],y[i],z[i]);
         }
         delete [] x;
         delete [] y;
         delete [] z;
      }
      // loop over atoms: each of them has a shpere around it
      for(unsigned int i=0; i<molsize; i++){
         // loop over the points on the sphere and add it respectively.
         for (unsigned int it=0; it<point.size(); it++){
            // shift and scale the coordinate as needed
            point[it]*=scale;
            point[it]+=geometry[i];
            if (i==closest(geometry, point[it])){
              // out<<"mesh:  "<<point[it](0)<<"  ";
              // out<<point[it](1)<<"  ";
              // out<<point[it](2)<<std::endl;
               mesh.add_point( point[it] );
            }
            // shift and scale back for use in the next round.
            point[it]-= geometry[i];
            point[it]/= scale;
         }
      }
   }
   // now, the point set is in mesh and we can call triangulate_pointset().
   #else
   // Avoid compiler warnings
   libmesh_ignore(mesh);
   libmesh_ignore(lower_limit);
   libmesh_ignore(upper_limit);
   #endif // LIBMESH_HAVE_TETGEN
}
