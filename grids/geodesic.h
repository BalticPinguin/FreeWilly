#include <math.h>
#include <iostream>

unsigned int start_4(double* x, double* y, double* z, bool** neighbours);
unsigned int start_6(double* x, double* y, double* z, bool** neighbours);
//unsgned int start_12(double* x, double* y, double* z, double* neighbours, int* n);
//unsigned start_18(double* x, double* y, double* z, double* neighbours, int* n);
unsigned int num_neighbour(int n);
void iterate(double*x, double*y, double*z, bool** neighbours, unsigned int* n, unsigned int m);
void gen_grid(double* x, double* y, double* z, const unsigned int iterations, const unsigned int nval);
unsigned int point_size(int n, int iterations);
//void iterate2(double*x, double*y, double*z, bool** neighbours, unsigned int* n, unsigned int m);
