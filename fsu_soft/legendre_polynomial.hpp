#ifndef LEGENDRE_H
#define LEDENGRE_H
char digit_to_ch ( int i );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
std::string i4_to_s ( int i );
void imtqlx ( int n, double d[], double e[], double z[] );
double *p_exponential_product ( int p, double b );
double p_integral ( int n );
double *p_polynomial_coefficients ( int n );
double *p_polynomial_prime ( int m, int n, double x[] );
double *p_polynomial_prime2 ( int m, int n, double x[] );
double *p_polynomial_value ( int m, int n, double x[] );
void p_polynomial_values ( int &n_data, int &n, double &x, double &fx );
double *p_polynomial_zeros ( int nt );
double *p_power_product ( int p, int e );
void p_quadrature_rule ( int nt, double t[], double wts[] );
double *pm_polynomial_value ( int mm, int n, int m, double x[] );
void pm_polynomial_values ( int &n_data, int &n, int &m, double &x, double &fx );
double *pmn_polynomial_value ( int mm, int n, int m, double x[] );
void pmn_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx );
double *pmns_polynomial_value ( int mm, int n, int m, double x[] );
void pmns_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx );
double *pn_pair_product ( int p );
double *pn_polynomial_coefficients ( int n );
double *pn_polynomial_value ( int m, int n, double x[] );
void pn_polynomial_values ( int &n_data, int &n, double &x, 
  double &fx );

#endif
