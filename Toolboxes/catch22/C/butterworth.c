//
//  butterworth.c
//  
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#include <stdlib.h>
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

#include "helper_functions.h"
#include "butterworth.h"

#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Imaginary_I * (double)(y)))
#endif

void poly(cplx x[], int size, cplx out[])
{
    /* Convert roots x to polynomial coefficients */
    
    // initialise
    #if defined(__GNUC__) || defined(__GNUG__)
        out[0] = 1;
        for(int i=1; i<size+1; i++){
            out[i] = 0;
        }
    #elif defined(_MSC_VER)
        cplx c1 = { 1, 0 };
        out[0] = c1;//1;
        for(int i=1; i<size+1; i++){
            c1._Val[0] = 0;
            out[i] = c1;
        }
    #endif
    
    
    cplx * outTemp = malloc((size+1)* sizeof(cplx));
    
    for(int i=1; i<size+1; i++){
        
        // save old out to not reuse some already changed values
        //(can be done more efficiently by only saving one old value, but who cares, pole number is usually ~4)
        for(int j=0; j<size+1; j++){
            outTemp[j] = out[j];
        }
        
        for(int j=1; j<i+1; j++){
            
            cplx temp1 = _Cmulcc(x[i - 1], outTemp[j - 1]);//x[i - 1] * outTemp[j - 1];
            cplx temp2 = out[j];
            out[j] = _Cminuscc(temp2, temp1);//x[i - 1] * outTemp[j - 1];
            
        }
        
    }
    
}

void filt(double y[], int size, double a[], double b[], int nCoeffs, double out[]){
    
    /* Filter a signal y with the filter coefficients a and b, output to array out.*/
    
    double offset = y[0];
    
    for(int i = 0; i < size; i++){
        out[i] = 0;
        for(int j = 0; j < nCoeffs; j++){
            if(i - j >= 0)
            {
                out[i] += b[j]*(y[i-j]-offset);
                out[i] -= a[j]*out[i-j];
            }
            else{
                out[i] += 0; //b[j]*offset; // 'padding'
                out[i] -= 0; //a[j]*offset;
            }
        }
    }
    
    for(int i = 0; i < size; i++){
        out[i] += offset;
    }
}

void reverse_array(double a[], int size){
    
    /* Reverse the order of the elements in an array. Write back into the input array.*/
    
    double temp;
    for(int i = 0; i < size/2; i++){
        temp = a[i];
        a[i] = a[size-i-1];
        a[size-1-i] = temp;
        /*
        printf("indFrom = %i, indTo = %i\n", i, size-1-i);
        for(int i=0; i < size; i++){
            printf("reversed[%i]=%1.3f\n", i, a[i]);
        }
         */
    }
}

void filt_reverse(double y[], int size, double a[], double b[], int nCoeffs, double out[]){
    
    /* Filter a signal y with the filter coefficients a and b _in reverse order_, output to array out.*/
    
    double * yTemp = malloc(size * sizeof(double));
    for(int i = 0; i < size; i++){
        yTemp[i] = y[i];
    }
    
    /*
    for(int i=0; i < size; i++){
        printf("yTemp[%i]=%1.3f\n", i, yTemp[i]);
    }
     */
    
    reverse_array(yTemp, size);
    
    /*
    for(int i=0; i < size; i++){
        printf("reversed[%i]=%1.3f\n", i, yTemp[i]);
    }
     */
    
    double offset = yTemp[0];
    
    for(int i = 0; i < size; i++){
        out[i] = 0;
        for(int j = 0; j < nCoeffs; j++){
            if(i - j >= 0)
            {
                out[i] += b[j]*(yTemp[i-j]-offset);
                out[i] -= a[j]*out[i-j];
            }
            else{
                out[i] += 0; //b[j]*offset; // 'padding'
                out[i] -= 0; //a[j]*offset;
            }
        }
    }
    
    for(int i = 0; i < size; i++){
        out[i] += offset;
    }
    
    reverse_array(out, size);
    
    free(yTemp);
    
}

/*
void butterworthFilter(const double y[], int size, const int nPoles, const double W, double out[]){
    
    double PI = 3.14159265359;

    double V = tan(W * PI/2);
    cplx * Q = malloc(nPoles * sizeof(cplx));
    
    for(int i = 0; i<nPoles; i++){
        cplx tmp1 = { 0, PI / 2 };
        cplx tmp2 = { nPoles, 0 };
        Q[i] =  conj(cexp((_Cdivcc(tmp1, tmp2)) * ((2 + nPoles - 1) + 2*i)));
        //printf("Q[%i]= real %1.3f imag %1.3f\n", i, creal(Q[i]), cimag(Q[i]));
        
    }
    
    double Sg = pow(V,nPoles);
    cplx * Sp = malloc(nPoles * sizeof(cplx));
    
    for(int i = 0; i<nPoles; i++){
        Sp[i] = V * Q[i];
        //printf("Sp[%i]= real %1.3f imag %1.3f\n", i, creal(Sp[i]), cimag(Sp[i]));
    }
    
    cplx * P = malloc(nPoles * sizeof(cplx));
    cplx * Z = malloc(nPoles * sizeof(cplx));
    
    cplx prod1mSp = 1; // %1 - Sp[0];
    
    // Bilinear transform for poles, fill zeros, compute products
    for(int i = 0; i<nPoles; i++){
        P[i] = (1 + Sp[i]) / (1 - Sp[i]);
        Z[i] = -1;
        
        // printf("P[%i]= real %1.3f imag %1.3f\n", i, creal(P[i]), cimag(P[i]));
        // printf("Z[%i]= real %1.3f imag %1.3f\n", i, creal(Z[i]), cimag(Z[i]));
        
        prod1mSp *= (1 - Sp[i]);
        
        //if(i > 0){
        //    prod1mSp *= (1 - Sp[i]);
        //}
         
    }

    double G = creal(Sg / prod1mSp);
    
    cplx * Zpoly = malloc((nPoles+1) * sizeof(cplx));
    cplx * Ppoly = malloc((nPoles+1) * sizeof(cplx));
    
    // polynomial coefficients from poles and zeros for filtering
    poly(Z, nPoles, Zpoly);
    
    //for(int i = 0; i < nPoles+1; i++){
    //    printf("Zpoly[%i]= %1.3f + %1.3f i\n", i, creal(Zpoly[i]), cimag(Zpoly[i]));
    //}
    
    poly(P, nPoles, Ppoly);
    
    //for(int i = 0; i < nPoles+1; i++){
    //    printf("Ppoly[%i]= %1.3f + %1.3f i\n", i, creal(Ppoly[i]), cimag(Ppoly[i]));
    //}
    
    
    // coeffs for filtering
    double * b = malloc((nPoles+1) * sizeof(double)); // zeros
    double * a = malloc((nPoles+1) * sizeof(double)); // poles
    
    for(int i = 0; i<nPoles+1; i++){
        b[i] = G * creal(Zpoly[i]);
        a[i] = creal(Ppoly[i]);
        //printf("a[%i]=%1.8f, b[%i]=%1.8f\n", i, a[i], i, b[i]);
    }
    
    
    // this is the simple solution, one forward filtering, phase not preserved.
    //filt(y, size, a, b, nPoles, out);
    
    // from here, better way of filtering. Forward and backward.
    // pad to both sides to avoid end-transients
    int nfact = 3*nPoles;
    double * yPadded = malloc((size + 2*nfact) * sizeof(double));
    
    for(int i = 0; i < nfact; i ++)
    {
        yPadded[i] = 2*y[0] - y[nfact-i];
        yPadded[nfact + size + i] = 2*y[size-1] - y[size-2-i];
        
    }
    for(int i = 0; i < size; i ++)
    {
        yPadded[nfact+i] = y[i];
    }
    
    // filter in both directions
    double * outPadded = malloc((size + 2*nfact) * sizeof(double));
    filt(yPadded, (size + 2*nfact), a, b, nPoles, outPadded);
    
    //for(int i=0; i < (size + 2*nfact); i++){
    //    printf("filtPadded[%i]=%1.3f\n", i, outPadded[i]);
    //}
     
    filt_reverse(outPadded, (size + 2*nfact), a, b, nPoles, outPadded);
    
    
    // for(int i=0; i < (size + 2*nfact); i++){
    //    printf("filtfiltPadded[%i]=%1.3f\n", i, outPadded[i]);
    //}
     
     
    for(int i = 0; i < size; i ++){
        out[i] = outPadded[nfact + i];
    }
    
    
    //for(int i=0; i < size; i++){
    //    printf("out[%i]=%1.3f\n", i, out[i]);
    //}
     
    
    free(Q);
    free(Sp);
    free(P);
    free(Z);
    free(a);
    free(b);
    free(Zpoly);
    free(Ppoly);
    free(outPadded);
    
}

*/