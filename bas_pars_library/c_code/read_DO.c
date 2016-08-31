/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bas_pars v.0.1 - fork of the AUGER program                          !
!                                                                     !
! code example that shows how to get the ATOM, BASIS and ENV arrays   !
!      using the Bas_pars library.                                    !
!                                                                     !
! Author:  Gilbert Grell, University of Rostock                       !
! Date:    09.08.2016                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <stdio.h>
#include <math.h> // need only exp() from it.
#include <vector>

extern void pass_arrays(char *file, int *length, int *atom, int *basis, double *env);
extern void pass_parameters(char *file, int *length, int *natom, int *nbas, int *nbasf, int *env_bas_dim, int *ptr_env_start);
extern char ang_label(int * ang);
extern void c_norm_contr_coeff(int *atom, int *basis, double *env);
extern void global_par_allocate();
extern void global_par_deallocate();
extern void pass_dyor(char *file, int *length, double* dyor); 
// skip the character arrays for the moment. due to principle incompatibility with c.

double DOval(double x, double y, double z, double* DyOr,double *env,int* basis,int* atom, int nbas);
double solHar(double x,double y,double z, unsigned int l, int m);
void printDO(double* DyOr, int nbasf);

int main( int argc, char **argv)
{
   //define the constants for accessing the atom, basis and env arrays:
   //atom array
   const int atom_slots = 6;
   const int charge_of = 0;
   const int ptr_coord = 1;
   //  const int nuc_mod = 2; //not used by my code
   //  const int ptr_zeta = 3; // not used by my code
   //basis array
   const int basis_slots = 8;
   const int atom_of = 0;
   const int ang_of = 1;
   const int nprim_of = 2;
   const int nctr_of = 3;
   //  const int kappa_of = 4; // not used by my code
   const int ptr_exp = 5;
   const int ptr_coeff = 6;
   
   // get the filename from command line:
   // check the number of command line arguments:
   
   if(argc != 2) //(the first argument is the invocation of the program itself)
   {
      fprintf(stdout, "please give the name of the inputfile as the only commandline argument!\n");
      fprintf(stdout, "Program aborts now!\n");
      return 0;
   }
   
   // count the number of characters:
   int n = 0;
   for(int i = 0; argv[1][i] != '\0'; i++)
   {
      n++;
   }
   
   // call the fortran routine pass_parameters() from the parse_inp module to get the parameters
   int natom, nbas, nbasf, env_bas_dim, ptr_env_start;
   pass_parameters(argv[1], &(n), &(natom), &(nbas), &(nbasf), &(env_bas_dim), &(ptr_env_start));
   
   //fprintf(stdout,"natom = %d nbas =  %d, nbasf = %d, env_bas_dim = %d, ptr_env_start = %d\n", natom, nbas,nbasf, env_bas_dim, ptr_env_start);
   // allocate the environmental arrays:
   int basis[nbas*basis_slots];
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
   int atom[natom*atom_slots];
   // atom= (a*6+b)
   // a: number of atoms
   // b=0  charge
   // b=1  *env: coordinates
   // b=2  nuc. model (not used)
   // b=3  *zeto (not used)
   // b=4  ?? unused
   // b=5  ?? unused
   double env[ptr_env_start+3*natom+env_bas_dim];
   // #Block          length      content
   // 1                20           ???
   // 2                3*Natom      geometry env[atom*3+coord]
   // 3 a) basis  (atom, ang.mom.)
   //         i)      #primitives  alpha_j
   //        ii)      #contraced   contr.-coefficients 
   // 4
   double DyOr[nbasf];
   double energy; //--> get from pass_dyor
   double normDO; //--> get from pass_dyor
   double x,y,z;
   
   // call the fortran routine pass_arrays() from the parse_inp module to get the arrays
   pass_arrays(argv[1], &(n), atom, basis, env);
   // normalize the contraction coefficients
   //c_norm_contr_coeff(atom, basis, env); TODO
   
   // read the coefficients of the Dyson Orbital:
   pass_dyor(argv[1], &n, DyOr); 
   //printDO(DyOr, nbasf);
   int oldatom=7;
   int expoff=ptr_env_start+3*natom;
   int k=0, l;
   double c_j;

   std::vector<double> do_j;
   std::

   MakeVectors();

   //setup a test-grid to evaluate the function values at
   for(unsigned int i=0; i<1; i++){
      for(unsigned int j=0; j<1; j++){
         for(unsigned int k=0; k<1; k++){
            x=-5.+i/10.;
            y=-5.+j/10.;
            z=-5.+k/10.;
            // evaluate data at grid and printh the values.
            //printf("%f   %f   %f   %g\n",x,y,z,DOval(x, y, z, DyOr, env, basis, atom, nbas));
            DOval(x, y, z, DyOr, env, basis, atom, nbas);
         }
      }
   }
   
   return 0;
}

double DOval(double x, double y, double z, double* DyOr,double *env,int* basis,int* atom, int nbas){
   double do_val=0;
   double radial, coeff, ylm;
   int k=0, b,j,n, l,m;
   const int basis_slots = 8;
   for(b=0; b<nbas;b++){                      // loop over bases
      l=basis[b*basis_slots+1];                    // each basis has fixed ang. momentum
      for(m=-l; m<=l; m++){                         // loop over m
         radial=0;
         ylm=solHar(x,y,z, l, m);
         for(j=0; j<basis[b*basis_slots+2]; j++ ){ // loop over primitives
            coeff=0;
            for(n=0; n<basis[b*basis_slots+3]; n++){  // quantum number=n+l
               coeff+=DyOr[k]*env[basis[b*basis_slots+6]+n*basis[b*basis_slots+2]+n];
               printf(" %d   %d   %d   %g   %e \n",n+l, l, m, DyOr[k], env[basis[b*basis_slots+6]+n*basis[b*basis_slots+2]+n]);
               k++;
            }
            coeff;
            radial+=coeff*exp(-env[basis[b*basis_slots+5]+j]*(x*x+y*y+z*z));
         }
         do_val+=radial*ylm;
      }
   }
   printf("\n\n");
   return do_val;
}
   //  b=0  index of atom
   //  b=1  angular mementum
   //  b=2  #primitive functions
   //  b=3  #contr. functions functions
   //  b=4  Kappa-offset --> not needed
   //  b=5  *env: start alpha
   //  b=6  *env: start contractions
   //  b=7  not used ?

void printDO(double* DyOr, int nbasf){
   for(unsigned int i=0; i<nbasf; i++){
      printf("%d  %f \n",i,DyOr[i]);
   }
}

double solHar(double x,double y,double z, unsigned int l, int m){
   //////////////////////////////////////////////////////////////////////////////
   /* please find the definition of the real spherical harmonics in             *
   *  https://en.wikipedia.org/wiki/Table_of_spherical_harmonics                *
   *  where the normalisation and r^l term change since I use solid harmonics:  *
   *  https://en.wikipedia.org/wiki/Solid_harmonics                             */
   //////////////////////////////////////////////////////////////////////////////
   // check the normalisation!! Something seems to be left here...
   switch(l){
   case 0:
   {
      return 1;
   }
   case 1:
   {
      switch(m){
      case -1:
         return y;
      case 0:
         return z;
      case 1:
         return x;
      default:
         return -300e30;
      }
   }
   case 2:
   {
      switch(m){
      case -2:
         // *sqrt(3)
         return x*y*1.73205080756888;
      case -1:
         return y*z*1.73205080756888;
      case 0:
         return (-x*x-y*y+2*z*z);
      case 1:
         return z*x*1.73205080756888;
      case 2:
         // *sqrt(3)/2
         return (x*x-y*y)*0.866025403784439;
      default:
         return -300e30;
      }
   }
   case 3:
   {
      switch(m){
      case -3:
         // sqrt(5/8)
         return y*(3*x*x-y*y)*0.590043589926644;
      case -2:
         // sqrt(15)
         return x*y*z*3.87298334620742;
      case -1:
         // sqrt(3/8)
         return y*(4*z*z-x*x-y*y)*0.612372435695794;
      case 0:
         return z*(2*z*z-3*x*x-3*y*y)*0.5;
      case 1:
         return x*(4*z*z-x*x-y*y)*0.612372435695794;
      case 2:
         return z*(x*x-y*y)*3.87298334620742;
      case 3:
         return x*(x*x-3*y*y)*0.590043589926644;
      default:
         return -300e30;
      }
   }
   case 4:
   {
      switch(m){
      case -4:
         // sqrt(35/4)
         return x*y*(x*x-y)*2.95803989154981;
      case -3:
         // sqrt(35/8)
         return y*z*(3*x*x-y*y)*2.09165006633519;
      case -2:
         // sqrt(5/4)
         return x*y*(6*z*z-x*x-y*y)*1.11803398874989;
      case -1:
         // * sqrt(5/8)
         return y*z*(4*z*z-3*x*x-3*y*y)*.790569415042095;
      case 0:
         // 35z^4-30z^2 r^2+3r^4
         // *1/8
         return 8*z*z*z*z-30*z*z*(y*y+x*x)+3*(x*x+y*y)*.125;
      case 1:
         // * sqrt(5/8)
         return x*z*(4*z*z-x*x-y*y)*.790569415042095;
      case 2:
         // * sqrt(5/16)
         return (x*x-y*y)*(6*z*z-x*x-y*y)*0.559016994374947;
      case 3:
         // sqrt(35/8)
         return x*z*(x*x-3*y*y)*2.09165006633519;
      case 4:
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
