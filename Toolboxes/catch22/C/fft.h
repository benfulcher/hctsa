#ifndef FFT_H
#define FFT_H
//#include <complex.h>

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

#include <stdlib.h>
#include <string.h>
#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Complex_I * (double)(y)))
#endif
extern void twiddles(cplx a[], int size);
// extern void _fft(cplx a[], cplx out[], int size, int step, cplx tw[]);
extern void fft(cplx a[], int size, cplx tw[]);
extern void ifft(cplx a[], int size, cplx tw[]);
#endif
