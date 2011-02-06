#ifndef HYPERGEOM_H
#define HYPERGEOM_H
#include <inttypes.h>

#define Integer int64_t

double alnfac ( int n );
double alngam ( double xvalue, int *ifault );
double alnorm ( double x, bool upper );
double chyper ( bool point, int kk, int ll, int mm, int nn, int *ifault );
void hypergeometric_cdf_values ( int *n_data, int *sam, int *suc, int *pop, 
  int *n, double *fx );
void hypergeometric_pdf_values ( int *n_data, int *sam, int *suc, int *pop, 
  int *n, double *fx );
int i4_min ( int i1, int i2 );
double r8_abs ( double x );
double r8_max ( double x, double y );
void timestamp ( void );

#define MIN(x,y) ((x)<(y)?(x):(y))

#endif

