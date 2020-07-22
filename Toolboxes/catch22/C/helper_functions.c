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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "stats.h"

// compare function for qsort, for array of doubles
static int compare (const void * a, const void * b)
{
    if (*(double*)a < *(double*)b) {
        return -1;
    } else if (*(double*)a > *(double*)b) {
        return 1;
    } else {
        return 0;
    }
}

// wrapper for qsort for array of doubles. Sorts in-place
void sort(double y[], int size)
{
    qsort(y, size, sizeof(*y), compare);
}

// linearly spaced vector
void linspace(double start, double end, int num_groups, double out[])
{
    double step_size = (end - start) / (num_groups - 1);
    for (int i = 0; i < num_groups; i++) {
        out[i] = start;
        start += step_size;
    }
    return;
}

double quantile(const double y[], const int size, const double quant)
{   
    double quant_idx, q, value;
    int idx_left, idx_right;
    double * tmp = malloc(size * sizeof(*y));
    memcpy(tmp, y, size * sizeof(*y));
    sort(tmp, size);
    
    /*
    for(int i=0; i < size; i++){
        printf("y[%i]=%1.4f\n", i, y[i]);
    }
    for(int i=0; i < size; i++){
        printf("sorted[%i]=%1.4f\n", i, tmp[i]);
    }
     */
    
    // out of range limit?
    q = 0.5 / size;
    if (quant < q) {
        return tmp[0]; // min value
    } else if (quant > (1 - q)) {
        return tmp[size - 1]; // max value
    }
    
    quant_idx = size * quant - 0.5;
    idx_left = (int)floor(quant_idx);
    idx_right = (int)ceil(quant_idx);
    value = tmp[idx_left] + (quant_idx - idx_left) * (tmp[idx_right] - tmp[idx_left]) / (idx_right - idx_left);
    free(tmp);
    return value;
}

void binarize(const double a[], const int size, int b[], const char how[])
{   
    double m = 0.0;
    if (strcmp(how, "mean") == 0) {
        m = mean(a, size);
    } else if (strcmp(how, "median") == 0) {
        m = median(a, size);
    }
    for (int i = 0; i < size; i++) {
        b[i] = (a[i] > m) ? 1 : 0;
    }  
    return;
}

double f_entropy(const double a[], const int size)
{
    double f = 0.0;
    for (int i = 0; i < size; i++) {
        if (a[i] > 0) {
            f += a[i] * log(a[i]);
        }
    }
    return -1 * f;
}

void subset(const int a[], int b[], const int start, const int end)
{
    int j = 0;
    for (int i = start; i < end; i++) {
        b[j++] = a[i];
    }
    return;
}

#if defined(__GNUC__) || defined(__GNUG__)
    cplx _Cmulcc(const cplx x, const cplx y) {
        /*double a = x._Val[0];
        double b = x._Val[1];

        double c = y._Val[0];
        double d = y._Val[1];

        cplx result = { (a * c - b * d), (a * d + c * b) };
         */
        return x*y;
    }

    cplx _Cminuscc(const cplx x, const cplx y) {
        //cplx result = { x._Val[0] - y._Val[0], x._Val[1] - y._Val[1] };
        return x - y;
    }

    cplx _Caddcc(const cplx x, const cplx y) {
        // cplx result = { x._Val[0] + y._Val[0], x._Val[1] + y._Val[1] };
        return x + y;
    }

    cplx _Cdivcc(const cplx x, const cplx y) {
        /*
        double a = x._Val[0];
        double b = x._Val[1];

        double c = y._Val[0];
        double d = y._Val[1];

        cplx result = { (a*c + b*d) / (c*c + d*d), (b*c - a*d)/(c*c + d*d)};
         */
         
        return x / y;
    }

#elif defined(_MSC_VER)
    cplx _Cminuscc(const cplx x, const cplx y) {
        cplx result = { x._Val[0] - y._Val[0], x._Val[1] - y._Val[1] };
        return result;
    }

    cplx _Caddcc(const cplx x, const cplx y) {
        cplx result = { x._Val[0] + y._Val[0], x._Val[1] + y._Val[1] };
        return result;
    }

    cplx _Cdivcc(const cplx x, const cplx y) {
        double a = x._Val[0];
        double b = x._Val[1];

        double c = y._Val[0];
        double d = y._Val[1];

        cplx result = { (a*c + b*d) / (c*c + d*d), (b*c - a*d)/(c*c + d*d)};

        return result;
    }
#endif
