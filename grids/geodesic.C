#include <math.h>
#include <iostream>

unsigned int start_4(double* x, double* y, double* z, bool** neighbours){
   x[0]= 0.0000000000000000; y[0]= 0.0000000000000000; z[0]= 1.00000000000000000;
   x[1]= 0.0000000000000000; y[1]= 0.9428090415820634; z[1]=-0.33333333333333333;
   x[2]= 0.4714045207910317; y[2]=-0.4714045207910317; z[2]=-0.33333333333333333;
   x[3]=-0.4714045207910317; y[3]=-0.4714045207910317; z[3]=-0.33333333333333333;
   // is point i with j a nearest neighbour?
   neighbours[0][1]=true;
   neighbours[0][2]=true;
   neighbours[0][3]=true;

   neighbours[1][0]=true;
   neighbours[1][2]=true;
   neighbours[1][3]=true;

   neighbours[2][0]=true;
   neighbours[2][1]=true;
   neighbours[2][3]=true;

   neighbours[3][0]=true;
   neighbours[3][1]=true;
   neighbours[3][2]=true;
   return 4;
}

unsigned int start_6(double* x, double* y, double* z, bool** neighbours){
   x[0]= 1.000; y[0]=0.0; z[0]=0.0;
   x[1]=-1.000; y[1]=0.0; z[1]=0.0;
   x[2]= 0.000; y[2]=1.0; z[2]=0.0;
   x[3]= 0.000; y[3]=-1.0;z[3]=0.0;
   x[4]= 0.000; y[2]=0.0; z[2]=1.0;
   x[5]= 0.000; y[3]=0.0; z[3]=-1.0;
   // is point i with j a nearest neighbour?
   neighbours[0][2]=true;
   neighbours[0][3]=true;
   neighbours[0][4]=true;
   neighbours[0][5]=true;

   neighbours[1][2]=true;
   neighbours[1][3]=true;
   neighbours[1][4]=true;
   neighbours[1][5]=true;

   neighbours[2][0]=true;
   neighbours[2][1]=true;
   neighbours[2][4]=true;
   neighbours[2][5]=true;

   neighbours[3][0]=true;
   neighbours[3][1]=true;
   neighbours[3][4]=true;
   neighbours[3][5]=true;

   neighbours[4][0]=true;
   neighbours[4][1]=true;
   neighbours[4][2]=true;
   neighbours[4][3]=true;

   neighbours[4][0]=true;
   neighbours[4][1]=true;
   neighbours[4][2]=true;
   neighbours[4][3]=true;
   return 4;
}

//unsgned int start_12(double* x, double* y, double* z, double* neighbours, int* n){
// return 5;
//}

//unsigned start_18(double* x, double* y, double* z, double* neighbours, int* n){
// return ?;
//}

unsigned int num_neighbour(int n){
   switch (n){
      case 4:
         return 4;
      case 6:
         return 4;
      case 12:
         return 5;       
      //case 18
      //   m=start_18(x,y,z,neighbours, n);
      default:
         return 0;
   }
}


void iterate(double*x, double*y, double*z, bool** neighbours, unsigned int* n, unsigned int m){
   bool** newneighbour;
   double* newx;
   double* newy;
   double* newz;
   newx = new double[(*n)*(m+1)];
   newy = new double[(*n)*(m+1)];
   newz = new double[(*n)*(m+1)];
   newneighbour = new bool*[(*n)*(m+1)];
   unsigned int k=0;
   double r, maxdist=0;
   double mindist=(x[0]-x[1])*(x[0]-x[1])+
                  (y[0]-y[1])*(y[0]-y[1])+
                  (z[0]-z[1])*(z[0]-z[1]);
   // first create all points:
   for(unsigned int i=0; i<(*n); i++){
      newx[k]=x[i];
      newy[k]=y[i];
      newz[k]=z[i];
      newneighbour[k] = new bool[m];
      for(unsigned int j=0; j<m; j++){
         newneighbour[k][j] = false;
      }
      k++;
      for(unsigned int j=0; j<(*n); j++){
         if (neighbours[i][j]){
            newx[k]=(x[i]+x[j])/2.;
            newy[k]=(y[i]+y[j])/2.;
            newz[k]=(z[i]+z[j])/2.;
            //normalise new point:
            r= sqrt(newx[k]*newx[k]+
                    newy[k]*newy[k]+
                    newz[k]*newz[k]);
            newx[k]*=r;
            newy[k]*=r;
            newz[k]*=r;

            //check the distances to the nearest neighbours:
            r=(newx[k]-x[i])*(newx[k]-x[i])+
              (newy[k]-y[i])*(newy[k]-y[i])+
              (newz[k]-z[i])*(newz[k]-z[i]);
            if(r<mindist)
               mindist=r;
            if(r>maxdist)
               maxdist=r;
            r=(newx[k]-x[j])*(newx[k]-x[j])+
              (newy[k]-y[j])*(newy[k]-y[j])+
              (newz[k]-z[j])*(newz[k]-z[j]);
            if(r<mindist)
               mindist=r;
            if(r>maxdist)
               maxdist=r;
            k++;
         }
      }
   }
   std::cout<<"extremal distances: "<<mindist<<"  "<<maxdist<<std::endl;
   //now, find their neighbours via the distance.
   for(unsigned int i=0; i<k; i++){
      for(unsigned int j=0; j<k; j++){
         r=(newx[j]-newx[i])*(newx[j]-newx[i])+
           (newy[j]-newy[i])*(newy[j]-newy[i])+
           (newz[j]-newz[i])*(newz[j]-newz[i]);
         if( r<maxdist && r>mindist)
            newneighbour[i][j]=true;
            newneighbour[j][i]=true;
      }
   }
   int sum;
   //Finally, check that each element has m neighbours again.
   for(unsigned int i=0; i<k; i++){
      sum=0;
      for(unsigned int j=0; j<k; j++){
         if (newneighbour[i][j])
            sum++;
      }
      if(sum+=m)
         std::cout<<"does not have correct number of neighours: ";
         std::cout<<i<<"  "<<sum<<"  "<<m<<std::endl;
   }
   *n=k;
   x=newx;
   y=newy;
   z=newz;
   neighbours=newneighbour;
}

void gen_grid(double* x, double* y, double* z, unsigned int iterations, unsigned int * n){
   unsigned int m;
   const unsigned nval=*n;
   bool** neighbours= new bool*[num_neighbour(nval)];
   for(unsigned int i=0; i<num_neighbour(nval); i++){
      neighbours[i] = new bool[nval];
      for(unsigned int j=0; j<nval; j++){
         neighbours[i][j]=false;
      }
   }
   switch (*n){
      // in this part, the meaning of n is changed.
      case 4:
         m= start_4(x,y,z,neighbours);
      case 6:
         m= start_6(x,y,z,neighbours);
      //case 12
      //   m=start_12(x,y,z,neighbours, n);
      //case 18
      //   m=start_18(x,y,z,neighbours, n);
      default:
         return;
   }
   for(unsigned int i=0; i<iterations; i++){
      iterate(x,y,z,neighbours,n, m);
   }
   return;
}


unsigned int point_size(unsigned int iterations, int n){
   int num_points=0;
   unsigned int m=0;
   switch (n){
      // in this part, the meaning of n is changed.
      case 4:
      {
         num_points=4;
         m=4;
      }
      case 6:
      {
         num_points=6;
         m=4;
      }
      case 12:
      {
         num_points=12;
         m=5;
      }
      case 18:
      {
         num_points=18;
         m=5;
      }
      default:
         return -1;
   }
   // in each iteration, per point m additional points are added:
   return num_points*pow(m,iterations);
}

//void iterate2(double*x, double*y, double*z, bool** neighbours, unsigned int* n, unsigned int m){
//   bool** newneighbour;
//   double* newx;
//   double* newy;
//   double* newz;
//   newx = new double[(*n)*(m+1)];
//   newy = new double[(*n)*(m+1)];
//   newz = new double[(*n)*(m+1)];
//   newneighbour = new bool*[m*(*n)*(m+1)];
//   unsigned int k=0, l=0;
//   double r, neighbourdist;
//   for(unsigned int i=0; i<*n; i++){
//      newx[k]=x[i];
//      newy[k]=y[i];
//      newz[k]=z[i];
//      newneighbour[k] = new bool[m];
//      k++;
//      for(unsigned int j=0; j<*n; j++){
//         if (neighbours[i][j]){
//            // add point of intersection:
//            newx[k]=(x[i]+x[j])/2.;
//            newy[k]=(y[i]+y[j])/2.;
//            newz[k]=(z[i]+z[j])/2.;
//            neighbourdist=((newx[k]-x[i])*(newx[k]-x[i])+
//                           (newy[k]-z[i])*(newy[k]-z[i])+
//                           (newz[k]-y[i])*(newz[k]-y[i]));
//            //normalise new point:
//            r= sqrt(newx[k]*newx[k]+
//                    newy[k]*newy[k]+
//                    newz[k]*newz[k]);
//            newx[k]*=r;
//            newy[k]*=r;
//            newz[k]*=r;
//            // find my neighbours:
//            {
//               // first add the 'parent' nodes
//               newneighbour[k][l]=true;
//               newneighbour[l][k]=true;
//               // now search tham as the nn of them:
//               for(unsigned int nn=0; nn<k; nn++){
//                  if(neighbours[i][nn] && 
//                    (newx[k]-newx[i])*(newx[k]-newx[i])+
//                    (newy[k]-newy[i])*(newy[k]-newz[i])+
//                    (newy[k]-newy[i])*(newy[k]-newz[i])+ <= neighbourdist)
//                  {
//                     newneighbour[k][nn]=true;
//                     newneighbour[nn][k]=true;
//                  }
//               }
//            }
//            k++;
//         }
//      }
//      l++;
//   }
//   x=newx;
//   y=newy;
//   z=newz;
//   neighbours=newneighbours;
//   *n=k;
//   return;
//}
