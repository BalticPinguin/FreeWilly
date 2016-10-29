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

void makeDO_j(std::vector<std::vector<double> >& do_j, std::vector<unsigned int>& l, std::vector<double>& alpha, int nbas,int* basis,double* env, double* DyOr);
double solHar(double x,double y,double z, unsigned int l, unsigned int m);
double solHar2(double x,double y,double z, unsigned int l, unsigned int m);
double evalDO(const std::vector<std::vector<double> >& do_j, const std::vector<unsigned int>& l, const std::vector<double>& alpha, const std::vector<libMesh::Node>& geometry, const libMesh::Point pt);

int factorial(unsigned int n);

#endif //READ_DO_H
