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

extern void pass_arrays(char *file, int *length, int *atom, int *basis, double *env);
extern void pass_parameters(char *file, int *length, int *natom, int *nbas, int *nbasf, int *env_bas_dim, int *ptr_env_start);
extern char ang_label(int * ang);
extern void c_norm_contr_coeff(int *atom, int *basis, double *env);
extern void  global_par_allocate();
extern void  global_par_deallocate();
// skip the character arrays for the moment. due to principle incompatibility with c.

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
  
  fprintf(stdout,"natom = %d nbas =  %d, nbasf = %d, env_bas_dim = %d, ptr_env_start = %d\n", natom, nbas,nbasf, env_bas_dim, ptr_env_start);
// allocate the environmental arrays:
  int basis[nbas*basis_slots];
  int atom[natom*atom_slots];
  double env[ptr_env_start+3*natom+env_bas_dim];
// allocate the fortran global arrays:
//  global_par_allocate();

// call the fortran routine pass_arrays() from the parse_inp module to get the arrays
  pass_arrays(argv[1], &(n), atom, basis, env);

// print information on the geometry, basis array and basis set:
  fprintf(stdout, "Geometry as read and converted to a.u.:\n");
  fprintf(stdout,"charge      X               Y              Z\n");
  for (int i = 0; i<natom; i++)
  {
    fprintf(stdout,"%2d   % 13.9f  % 13.9f  % 13.9f\n",atom[i*atom_slots+charge_of], env[atom[i*atom_slots + ptr_coord] + 0], env[atom[i*atom_slots + ptr_coord] + 1], env[atom[i*atom_slots + ptr_coord] + 2]);
  }
  fprintf(stdout, "\n Basis set as read and assigned to the atoms:\n");
  for(int i = 0; i<nbas; i++)
  {
    if((i>0) && (basis[i * basis_slots + atom_of] != basis[(i-1) * basis_slots + atom_of]))
    {
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "%3d (%2d -> %2d) %c\n", atom[basis[i*basis_slots + atom_of]*atom_slots + charge_of], basis[i*basis_slots + nprim_of], basis[i*basis_slots + nctr_of], ang_label(&(basis[i*basis_slots + ang_of])) );
    for(int j = 0; j< basis[i*basis_slots + nprim_of]; j++)
    {
      fprintf(stdout, "% 16.8f ", env[basis[i*basis_slots + ptr_exp] + j]);
      for ( int k = 0; k < basis[i*basis_slots + nctr_of]; k++)
      {
        fprintf(stdout, " % 16.8f ", env[basis[i*basis_slots + ptr_coeff] + j + k*basis[i*basis_slots + nprim_of]]);
      }
      fprintf(stdout, "\n");
    }
  }
// normalize the contraction coefficients
  c_norm_contr_coeff(atom, basis, env);

  fprintf(stdout, "\n Basis set with normalized coefficients:\n");
  for(int i = 0; i<nbas; i++)
  {
    if((i>0) && (basis[i * basis_slots + atom_of] != basis[(i-1) * basis_slots + atom_of]))
    {
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "%3d (%2d -> %2d) %c\n", atom[basis[i*basis_slots + atom_of]*atom_slots + charge_of], basis[i*basis_slots + nprim_of], basis[i*basis_slots + nctr_of], ang_label(&(basis[i*basis_slots + ang_of])) );
    for(int j = 0; j< basis[i*basis_slots + nprim_of]; j++)
    {
      fprintf(stdout, "% 16.8f ", env[basis[i*basis_slots + ptr_exp] + j]);
      for ( int k = 0; k < basis[i*basis_slots + nctr_of]; k++)
      {
        fprintf(stdout, " % 16.8f ", env[basis[i*basis_slots + ptr_coeff] + j + k*basis[i*basis_slots + nprim_of]]);
      }
      fprintf(stdout, "\n");
    }
  }

// deallocate the fortran global arrays:
//  global_par_deallocate();

  return 0;
}
