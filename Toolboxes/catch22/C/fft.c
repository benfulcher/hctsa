#include <stdlib.h>
#include <string.h>
#include <math.h>

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

#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Imaginary_I * (double)(y)))
#endif

#include "helper_functions.h"

void twiddles(cplx a[], int size)
{

    double PI = 3.14159265359;

    for (int i = 0; i < size; i++) {
        cplx tmp = { 0, -PI * i / size };
        a[i] = cexp(tmp);
        //a[i] = cexp(-I * M_PI * i / size);
    }
}

static void _fft(cplx a[], cplx out[], int size, int step, cplx tw[])
{   
    if (step < size) {
        _fft(out, a, size, step * 2, tw);
        _fft(out + step, a + step, size, step * 2, tw);

        for (int i = 0; i < size; i += 2 * step) {
            //cplx t = tw[i] * out[i + step];
            cplx t = _Cmulcc(tw[i], out[i + step]);
            a[i / 2] = _Caddcc(out[i], t);
            a[(i + size) / 2] = _Cminuscc(out[i], t);
        }
    }
}

void fft(cplx a[], int size, cplx tw[])
{
    cplx * out = malloc(size * sizeof(cplx));
    memcpy(out, a, size * sizeof(cplx));
    _fft(a, out, size, 1, tw);
    free(out);
}
