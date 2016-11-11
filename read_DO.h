#ifndef READ_DO_H
#define READ_DO_H

#include <iostream>
#include <stdio.h>
#include <math.h> // need only exp() from it.
#include <string.h>
#include <vector>

#include "fsu_soft/legendre_polynomial.hpp"
#include "fsu_soft/sphere_lebedev_rule.hpp" // compute theta,phi from x,y,z.

#include "libmesh/node.h"
#include "libmesh/point.h"
#ifdef __cplusplus
extern "C"{
#endif 

extern void pass_arrays(const char* , int*, int*, int*, double*);
extern void pass_parameters(const char*, int*, int*, int*, int*, int*, int*);
extern char ang_label(int * ang);
extern void c_norm_contr_coeff(int *atom, int *basis, double *env);
extern void global_par_allocate();
extern void global_par_deallocate();
extern void pass_dyor(const char*, int*, double*, double*, double*); 
// skip the character arrays for the moment. due to principle incompatibility with c.

#ifdef __cplusplus
}
#endif

void getDyson(const char *filename, int namelength, std::vector<std::vector<double> >& do_j, std::vector<unsigned int>& l, std::vector<double>& alpha, double&  energy, double& normDO);
void makeDO_j(std::vector<std::vector<double> >& do_j, std::vector<unsigned int>& l, std::vector<double>& alpha, int nbas,int* basis,double* env, double* DyOr);
double solHar(double x,double y,double z, unsigned int l, unsigned int m);
double solHar2(double x,double y,double z, unsigned int l, unsigned int m);
double eval_DO(const std::vector<std::vector<double> >& do_j, const std::vector<unsigned int>& l, const std::vector<double>& alpha, const std::vector<libMesh::Node>& geometry, const libMesh::Point pt);
std::vector<libMesh::Node> getGeometry(std::string fname);

int factorial(unsigned int n);


class DOrbit{
public:
   DOrbit(std::string DO_file)
   {
     geometry= getGeometry(DO_file);
     //es.parameters.get<std::string>("DO_file");
     const char* filename=DO_file.c_str();
     int namelength=strlen(filename);
     getDyson(filename, namelength, do_j, l, alpha, energy, normDO);
   }
   
   ~DOrbit()
   {}

   libMesh::Real evalDO(libMesh::Point q_point){
      return eval_DO(do_j, l, alpha, geometry, q_point);
   }

   libMesh::Real get_norm()
   {return normDO;}

   libMesh::Real get_energy()
   {return energy;}

   std::vector<libMesh::Node> geometry;

private:
   std::vector<std::vector<double> > do_j;
   std::vector<unsigned int> l;
   std::vector<double> alpha;
   libMesh::Real energy=0, normDO=0;
};

#endif //READ_DO_H
