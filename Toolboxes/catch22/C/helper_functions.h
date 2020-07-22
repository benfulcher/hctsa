#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "stats.h"

#if __cplusplus
#   include <complex>
typedef std::complex< double > cplx;
#else
#   include <complex.h>
#if defined(__GNUC__) || defined(__GNUG__)
typedef double complex cplx;
#elif defined(_MSC_VER)
typedef _Dcomplex cplx;
#endif
#endif

extern void linspace(double start, double end, int num_groups, double out[]);
extern double quantile(const double y[], const int size, const double quant);
extern void sort(double y[], int size);
extern void binarize(const double a[], const int size, int b[], const char how[]);
extern double f_entropy(const double a[], const int size);
extern void subset(const int a[], int b[], const int start, const int end);

extern cplx _Cminuscc(const cplx x, const cplx y);
extern cplx _Caddcc(const cplx x, const cplx y);
extern cplx _Cdivcc(const cplx x, const cplx y);
#if defined(__GNUC__) || defined(__GNUG__)
extern cplx _Cmulcc(const cplx x, const cplx y);
#endif

#endif
