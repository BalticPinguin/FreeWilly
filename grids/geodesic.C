#include "geodesic.h"

unsigned int start_4(double* x, double* y, double* z, int** neighbours){
   //x[0]= 0.0000000000000000; y[0]= 0.0000000000000000; z[0]= 1.00000000000000000;
   //x[1]= 0.0000000000000000; y[1]= 0.9428090415820634; z[1]=-0.33333333333333333;
   //x[2]= 0.4714045207910317; y[2]=-0.4714045207910317; z[2]=-0.33333333333333333;
   //x[3]=-0.4714045207910317; y[3]=-0.4714045207910317; z[3]=-0.33333333333333333;
   x[0]= 0.577350269189626; y[0]= 0.577350269189626; z[0]= 0.577350269189626;
   x[1]= 0.577350269189626; y[1]=-0.577350269189626; z[1]=-0.577350269189626;
   x[2]=-0.577350269189626; y[2]= 0.577350269189626; z[2]=-0.577350269189626;
   x[3]=-0.577350269189626; y[3]=-0.577350269189626; z[3]= 0.577350269189626;

   // is point i with j a nearest neighbour?
   neighbours[0][0]=1;
   neighbours[0][1]=2;
   neighbours[0][2]=3;

   neighbours[1][0]=0;
   neighbours[1][1]=2;
   neighbours[1][2]=3;
                    
   neighbours[2][0]=0;
   neighbours[2][1]=1;
   neighbours[2][2]=3;
                    
   neighbours[3][0]=0;
   neighbours[3][1]=1;
   neighbours[3][2]=2;
   return 3;
}

unsigned int start_6(double* x, double* y, double* z, int** neighbours){
   x[0]= 1.000; y[0]=0.0; z[0]=0.0;
   x[1]=-1.000; y[1]=0.0; z[1]=0.0;
   x[2]= 0.000; y[2]=1.0; z[2]=0.0;
   x[3]= 0.000; y[3]=-1.0;z[3]=0.0;
   x[4]= 0.000; y[2]=0.0; z[2]=1.0;
   x[5]= 0.000; y[3]=0.0; z[3]=-1.0;
   // is point i with j a nearest neighbour?
   neighbours[0][0]=2;
   neighbours[0][1]=3;
   neighbours[0][2]=4;
   neighbours[0][3]=5;
                    
   neighbours[1][0]=2;
   neighbours[1][1]=3;
   neighbours[1][2]=4;
   neighbours[1][3]=5;
                    
   neighbours[2][0]=0;
   neighbours[2][1]=1;
   neighbours[2][2]=4;
   neighbours[2][3]=5;
                    
   neighbours[3][0]=0;
   neighbours[3][1]=1;
   neighbours[3][2]=4;
   neighbours[3][3]=5;
                    
   neighbours[4][0]=0;
   neighbours[4][1]=1;
   neighbours[4][2]=2;
   neighbours[4][3]=3;
                    
   neighbours[4][0]=0;
   neighbours[4][1]=1;
   neighbours[4][2]=2;
   neighbours[4][3]=3;
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
         return 3;
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


void iterate(double*x, double*y, double*z, int** neighbours, unsigned int* n, unsigned int m){
   int ** newneighbour;
   unsigned int totlen=(*n)*(m+1);
   double* newx;
   double* newy;
   double* newz;
   newx = new double[totlen];
   newy = new double[totlen];
   newz = new double[totlen];
   newneighbour = new int*[totlen];
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
      newneighbour[k] = new int[m];
      for(unsigned int j=0; j<m; j++){
         newneighbour[k][j] = -1;
      }
      k++;
      for(unsigned int j=0; j<i; j++){
         if (neighbours[i][j]>=0){
            newx[k]=(x[i]+x[j])/2.;
            newy[k]=(y[i]+y[j])/2.;
            newz[k]=(z[i]+z[j])/2.;
            //normalise new point:
            r= sqrt(newx[k]*newx[k]+
                    newy[k]*newy[k]+
                    newz[k]*newz[k]);
            newx[k]/=r;
            newy[k]/=r;
            newz[k]/=r;

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
            newneighbour[k] = new int[m];
            for(unsigned int l=0; l<m; l++){
               newneighbour[k][l] = -1;
            }
            k++;
         }
      }
   }
   std::cout<<"extremal distances: "<<mindist<<"  "<<maxdist<<std::endl;
   //now, find their neighbours via the distance.
   for(unsigned int i=0; i<totlen; i++){
      k=0;
      std::cout<<newx[i]<<"  "<<newy[i]<<"  "<<newz[i]<<std::endl;
      for(unsigned int j=0; j<totlen; j++){
         if ( i==j) continue;
         r=(newx[j]-newx[i])*(newx[j]-newx[i])+
           (newy[j]-newy[i])*(newy[j]-newy[i])+
           (newz[j]-newz[i])*(newz[j]-newz[i]);
         if( r<=maxdist && r>=mindist){
            newneighbour[i][k]=j;
            k++;
         }
         //std::cout<<i<<"  "<<j<<"  "<<r<<"  "<<newneighbour[j][k]<<std::endl;
      }
      //Finally, check that each element has m neighbours again.
      for(unsigned int j=0; j<m; j++){
         if (newneighbour[i][j]<0){
            //std::cout<<"does not have correct number of neighours: ";
            //std::cout<<i<<"  "<<j<<"  "<<m<<std::endl;
         }
      }
   }

   *n=totlen;
   x=newx;
   y=newy;
   z=newz;
   neighbours=newneighbour;
   
   delete [] newx;
   delete [] newy;
   delete [] newz;
   for(unsigned int i=0; i<totlen; i++)
      delete [] newneighbour[i];
   delete [] newneighbour;
}

void gen_grid(double* x, double* y, double* z, const unsigned int iterations, const unsigned int nval){
   unsigned int m;
   unsigned int *n;
   *n=nval;
   int** neighbours= new int*[nval];
   for(unsigned int i=0; i<nval; i++){
      neighbours[i] = new int[num_neighbour(nval)];
      for(unsigned int j=0; j<num_neighbour(nval); j++){
         neighbours[i][j]=-1;
      }
   }
   switch (nval){
      // in this part, the meaning of n is changed.
      case 4:
         m= start_4(x,y,z,neighbours);
         break;
      case 6:
         m= start_6(x,y,z,neighbours);
         break;
      //case 12
      //   m=start_12(x,y,z,neighbours, n);
      //   break;
      //case 18
      //   m=start_18(x,y,z,neighbours, n);
      //   break;
      default:
         return;
   }
   for(unsigned int i=0; i<iterations; i++){
      iterate(x,y,z,neighbours,n, m);
   }
   for(unsigned int i=0; i<num_neighbour(nval); i++)
      delete [] neighbours[i];
   delete neighbours;
   return;
}


unsigned int point_size(int n, int iterations){
   int num_points=0;
   unsigned int m=0;
   switch (n){
      // in this part, the meaning of n is changed.
      case 4:
      {
         num_points=4;
         m=4;
         break;
      }
      case 6:
      {
         num_points=6;
         m=4;
         break;
      }
      case 12:
      {
         num_points=12;
         m=5;
         break;
      }
      case 18:
      {
         num_points=18;
         m=5;
         break;
      }
      default:
         return -1;
   }
   // in each iteration, per point m additional points are added:
   return num_points*pow(m,iterations);
}
