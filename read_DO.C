#include "read_DO.h"
/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
! Author:  Tobias Moehle, University of Rostock                       !
! Date:    30.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/


void getDyson(const char *filename, int namelength, std::vector<std::vector<double> >& do_j, std::vector<unsigned int>& l,std::vector<double>& alpha, double&  energy, double& normDO){
   //define the constants for accessing the atom, basis and env arrays:
   //atom array
   const int atom_slots = 6;
   //const int charge_of = 0;
   //const int ptr_coord = 1;
   //const int nuc_mod = 2; //not used by my code
   //const int ptr_zeta = 3; // not used by my code
   //basis array
   const int basis_slots = 8;
   //const int atom_of = 0;
   //const int ang_of = 1;
   //const int nprim_of = 2;
   //const int nctr_of = 3;
   //const int kappa_of = 4; // not used by my code
   //const int ptr_exp = 5;
   //const int ptr_coeff = 6;
   int* basis;
   int* atom;
   double* env;
   double* DyOr;
   
   // get the filename from command line:
   // check the number of command line arguments:
   
   // count the number of characters:
   int n = namelength;
   
   // call the fortran routine pass_parameters() from the parse_inp module to get the parameters
   int natom, nbas, nbasf, env_bas_dim, ptr_env_start;
   pass_parameters(filename, &(n), &(natom), &(nbas), &(nbasf), &(env_bas_dim), &(ptr_env_start));
   
   //fprintf(stdout,"natom = %d nbas =  %d, nbasf = %d, env_bas_dim = %d, ptr_env_start = %d\n", natom, nbas,nbasf, env_bas_dim, ptr_env_start);
   // allocate the environmental arrays:
   basis=new int[nbas*basis_slots];
   // basis= (a*8 +b) 
   //  a:   specifies the # of current basis
   //  b=0  index of atom
   //  b=1  angular mementum
   //  b=2  #primitive functions
   //  b=3  #contr. functions functions
   //  b=4  Kappa-offset --> not needed
   //  b=5  *env: start alpha
   //  b=6  *env: start contractions
   //  b=7  not used ?
   atom= new int[natom*atom_slots];
   // atom= (a*6+b)
   // a: number of atoms
   // b=0  charge
   // b=1  *env: coordinates
   // b=2  nuc. model (not used)
   // b=3  *zeto (not used)
   // b=4  ?? unused
   // b=5  ?? unused
   env=new double[ptr_env_start+3*natom+env_bas_dim];
   // #Block          length      content
   // 1                20           ???
   // 2                3*Natom      geometry env[atom*3+coord]
   // 3 basis  (atom, ang.mom.)
   //      a)      #primitives  alpha_j
   //      b)      #contraced   contr.-coefficients 
   // 4
   DyOr=new double[nbasf];
   
   // call the fortran routine pass_arrays() from the parse_inp module to get the arrays
   pass_arrays(filename, &(n), atom, basis, env);
   // normalize the contraction coefficients
   c_norm_contr_coeff(atom, basis, env);
   
   // read the coefficients of the Dyson Orbital:
   pass_dyor(filename, &n, DyOr, &energy, &normDO); 
   //int expoff=ptr_env_start+3*natom;
   
   makeDO_j(do_j, l, alpha, nbas, basis, env, DyOr);
   delete[] basis;
   delete[] atom;
   delete[] env;
   delete[] DyOr;
}

double eval_DO(const std::vector<std::vector<double> >& do_j, const std::vector<unsigned int>& l, const std::vector<double>& alpha, const std::vector<libMesh::Node>& geometry, const libMesh::Point pt){
   // evaluates the value of dyson orbital at x,y,z. I would consider 
   // this as quite efficient way.
   double angular, do_xyz=0;
   unsigned int i=0;
   double diff_x, diff_y, diff_z;
   for(unsigned int k=0; k<do_j.size(); k++){
      angular=0;
      diff_x=pt(0)-geometry[i](0);
      diff_y=pt(1)-geometry[i](1);
      diff_z=pt(2)-geometry[i](2);
      for(unsigned int m=0; m<2*l[k]+1; m++){
         angular+=solHar2(diff_x,diff_y,diff_z,l[k],m)*do_j[k][m]; 
      }
      do_xyz+=exp(-alpha[k]*(diff_x*diff_x+diff_y*diff_y+diff_z*diff_z))*angular;
      if (k>0 && l[k]<l[k-1]) i++;  // than next atom is considered.
   }
   return do_xyz;
}

void makeDO_j(std::vector<std::vector<double> >& do_j, std::vector<unsigned int>& l, std::vector<double>& alpha, int nbas, int* basis,double* env, double* DyOr){
   /////////////////////////////////////////////////////////////////
   /* This function fills the 2-D vector do_j with values          *
    *  DO(m,k(j))=sum_n do(n,m,l(j))* c(n,j)                       *
    *  where k runs over j for each bases after an other.          *
    * to have easier access to the coefficient times Dyson orbital.*
    * With this, the Dyson orbital can be expressed as             *
    * sum_k exp(-alpha_k*r^2) * sum_(m=-l)^l Y_lm DO[k][m]         */
    /////////////////////////////////////////////////////////////////
    // The goal is to make the expr. for indices easier and cancel 1 loop already.
   unsigned int k=0, NumBas;
   unsigned int atome=0; // the name is misleading.
   unsigned int basis_slots=8;
   do_j.clear();
   l.clear();
   alpha.clear();
   NumBas=0;
   for(unsigned int b=0; b<(unsigned)nbas; b++){
      NumBas+=basis[b*basis_slots+2];
   }
   do_j.resize(NumBas);
   l.resize(NumBas);
   alpha.resize(NumBas);
   for(unsigned int b=0; b<(unsigned)nbas; b++){
      for(unsigned int j=0; j<(unsigned)basis[b*basis_slots+2]; j++){
         // count k separately since j goes over different many steps.
         l[k]=basis[b*basis_slots+1];
         alpha[k]=env[basis[b*basis_slots+5]+j]; // bring alpha to continuous array.
         do_j[k].resize(2*l[k]+1);
         for(unsigned int m=0; m<2*l[k]+1; m++){
            // now, put different ns together into do_j.
            do_j[k][m]=0;
            for(unsigned int n=0; n<basis[b*basis_slots+3]; n++){ //quant. number: n+l
               do_j[k][m]+=DyOr[n+m*basis[b*basis_slots+3]+ atome]*env[basis[b*basis_slots+6]+n*basis[b*basis_slots+2]+j];
               // this follows my expectations, I expect...
               // printf("%d  %d  %d  %d  %d   %g ",
               //              basis[b*basis_slots],n+l[k]+1,l[k], m-l[k], n+m*basis[b*basis_slots+3]+atome,
               //              DyOr[n+m*basis[b*basis_slots+3]+ atome]);
               // printf(" %f  %g \n", env[basis[b*basis_slots+6]+n*basis[b*basis_slots+2]+j] , alpha[k]);
            }
         }
         k++;
      }
      atome+=(basis[b*basis_slots+3])*(2*l[k-1]+1);
   }
}

int factorial(unsigned int n){
   if(n>1)
      return n*factorial(n-1);
   return 1;
}

double solHar2(double x,double y,double z, unsigned int l, unsigned int mpl){
   //http://www.ppsloan.org/publications/StupidSH36.pdf
   int m=(int)(mpl-l); // m=m+l-l.     //12.5663706144=4*pi
   double K_lm=sqrt((2.*l+1.)*factorial(l-abs(m))/(12.5663706144*factorial(l+abs(m))));
   double* value;
   double r=sqrt(x*x+y*y+z*z), phival;
   double cos_theta[1];

   if( r<1e-12){
      cos_theta[0]=0;
      phival=0;
   }
   else{
      //theta[0]=acos ( z/r ); //0-> pi
      cos_theta[0]=z/r;  // cosTheta
      phival=atan2(y,x); //-pi -> pi
   }
   // I am not interested in spherical Harmonics but solid Harmonics:
   //  correct by factor r^l.
   if (l!=0){
      K_lm=K_lm*pow(r,l);
   }
   //theta[0]=cos(theta[0]); // make cos(theta) out of it.

   //evaluate associated legendre polynomial;
   // it is defined only for m>=0 here, therefore need to distinguish
   // three cases (following the convention in QM)
   // value is allocated in pm_polynomial_value()!
   if (m<0)
      value = pm_polynomial_value ( 1, l, -m, cos_theta);

   else
      value = pm_polynomial_value ( 1, l, m, cos_theta);
   double val=value[l];
   delete [] value;

   if (m>0){
      return 1.424214 *val*K_lm*cos(m*phival);
      }
   if (m<0){
      // negative sign due to convention in QM
      // for m<0 this doesn't come up.
      return 1.424214*val*K_lm*sin(-m*phival);
      }
   return val*K_lm;
}

double solHar(double x,double y,double z, unsigned int l, unsigned int m){
   // this uses self-written functions; maybe more efficient, but might contain bugs.
   //////////////////////////////////////////////////////////////////////////////
   /* please find the definition of the real spherical harmonics in             *
   *  https://en.wikipedia.org/wiki/Table_of_spherical_harmonics                *
   *  where the normalisation and r^l term change since I use solid harmonics:  *
   *  https://en.wikipedia.org/wiki/Solid_harmonics                            */
   //////////////////////////////////////////////////////////////////////////////
   // check the normalisation!! Something seems to be left here...
   // the quantum number m corresponds to m-l here.
   switch(l){
   case 0:
   {
      return 1;
   }
   case 1:
   {
      switch(m){
      case 0:
         return y;
      case 1:
         return z;
      case 2:
         return x;
      default:
         return -300e30;
      }
   }
   case 2:
   {
      switch(m){
      case 0:
         // *sqrt(3)
         return x*y*1.73205080756888;
      case 1:
         return y*z*1.73205080756888;
      case 2:
         return (-x*x-y*y+2*z*z);
      case 3:
         return z*x*1.73205080756888;
      case 4:
         // *sqrt(3)/2
         return (x*x-y*y)*0.866025403784439;
      default:
         return -300e30;
      }
   }
   case 3:
   {
      switch(m){
      case 0:
         // sqrt(5/8)
         return y*(3*x*x-y*y)*0.590043589926644;
      case 1:
         // sqrt(15)
         return x*y*z*3.87298334620742;
      case 2:
         // sqrt(3/8)
         return y*(4*z*z-x*x-y*y)*0.612372435695794;
      case 3:
         return z*(2*z*z-3*x*x-3*y*y)*0.5;
      case 4:
         return x*(4*z*z-x*x-y*y)*0.612372435695794;
      case 5:
         return z*(x*x-y*y)*3.87298334620742;
      case 6:
         return x*(x*x-3*y*y)*0.590043589926644;
      default:
         return -300e30;
      }
   }
   case 4:
   {
      switch(m){
      case 0:
         // sqrt(35/4)
         return x*y*(x*x-y)*2.95803989154981;
      case 1:
         // sqrt(35/8)
         return y*z*(3*x*x-y*y)*2.09165006633519;
      case 2:
         // sqrt(5/4)
         return x*y*(6*z*z-x*x-y*y)*1.11803398874989;
      case 3:
         // * sqrt(5/8)
         return y*z*(4*z*z-3*x*x-3*y*y)*.790569415042095;
      case 4:
         // 35z^4-30z^2 r^2+3r^4
         // *1/8
         return 8*z*z*z*z-30*z*z*(y*y+x*x)+3*(x*x+y*y)*.125;
      case 5:
         // * sqrt(5/8)
         return x*z*(4*z*z-x*x-y*y)*.790569415042095;
      case 6:
         // * sqrt(5/16)
         return (x*x-y*y)*(6*z*z-x*x-y*y)*0.559016994374947;
      case 7:
         // sqrt(35/8)
         return x*z*(x*x-3*y*y)*2.09165006633519;
      case 8:
         // sqrt(35)/8
         return x*x*(x*x-3*y*y)-y*y*(3*x*x-y*y)*0.739509972887452;
      default:
         return -300e30;
      }
   }
   default:
      return -300e30;
   }
   // never will be reached.
   return -300e30;
}
   
std::vector<libMesh::Node> getGeometry(std::string fname){
   const int atom_slots = 6;
   const int basis_slots = 8;
   const char* filename=fname.c_str();
   int n=strlen(filename);
   
   // call the fortran routine pass_parameters() from the parse_inp module to get the parameters
   int natom, nbas, nbasf, env_bas_dim, ptr_env_start;
   pass_parameters(filename, &(n), &(natom), &(nbas), &(nbasf), &(env_bas_dim), &(ptr_env_start));
   int* atom;
   int* basis;
   double* env;
   atom= new int [natom*atom_slots];
   basis= new int[nbas*basis_slots];
   env=new double[ptr_env_start+3*natom+env_bas_dim];
   
   // call the fortran routine pass_arrays() from the parse_inp module to get the arrays
   pass_arrays(filename, &(n), atom, basis, env);
   
   // write charge, geometry to vector:
   std::vector<libMesh::Node> geometry;
   for(int i=0; i<natom; i++){
      libMesh::Node tmpnd(env[atom[i*6+1]], env[atom[i*6+1]+1], env[atom[i*6+1]+2], atom[i*6]);
      geometry.push_back(tmpnd);
   }
   
   delete[] atom;
   delete[] basis;
   delete[] env;

   return geometry;
}
