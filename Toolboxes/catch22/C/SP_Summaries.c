//
//  SP_Summaries.c
//  
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#include "SP_Summaries.h"
#include "CO_AutoCorr.h"

int welch(const double y[], const int size, const int NFFT, const double Fs, const double window[], const int windowWidth, double ** Pxx, double ** f){
    
    double dt = 1.0/Fs;
    double df = 1.0/(nextpow2(windowWidth))/dt;
    double m = mean(y, size);
    
    // number of windows, should be 1
    int k = floor((double)size/((double)windowWidth/2.0))-1;
    
    // normalising scale factor
    double KMU = k * pow(norm_(window, windowWidth),2);
    
    double * P = malloc(NFFT * sizeof(double));
    for(int i = 0; i < NFFT; i++){
        P[i] = 0;
    }
    
    // fft variables
    cplx * F = malloc(NFFT * sizeof *F);
    cplx * tw = malloc(NFFT * sizeof *tw);
    twiddles(tw, NFFT);
    
    double * xw = malloc(windowWidth * sizeof(double));
    for(int i = 0; i<k; i++){
        
        // apply window
        for(int j = 0; j<windowWidth; j++){
            xw[j] = window[j]*y[j + (int)(i*(double)windowWidth/2.0)];
        }
        
        // initialise F (
        for (int i = 0; i < windowWidth; i++) {
            cplx tmp = { xw[i] - m, 0.0 };
            F[i] = tmp; // CMPLX(xw[i] - m, 0.0);
        }
        for (int i = windowWidth; i < NFFT; i++) {
            // F[i] = CMPLX(0.0, 0.0);
            cplx tmp = { 0.0, 0.0 };
            F[i] = tmp;
        }
        
        fft(F, NFFT, tw);
        /*
        for(int i = 0; i < NFFT; i++){
            printf("F1[%i] real: %1.3f, imag: %1.3f\n", i, creal(F[i]), cimag(F[i]));
        }
         */
        
        for(int l = 0; l < NFFT; l++){
            P[l] += pow(cabs(F[l]),2);
        }
        /*
        for(int i = 0; i < NFFT; i++){
            printf("P[%i]: %1.3f\n", i, P[i]);
        }
         */
        
    }
    
    int Nout = (NFFT/2+1);
    *Pxx = malloc(Nout * sizeof(double));
    for(int i = 0; i < Nout; i++){
        (*Pxx)[i] = P[i]/KMU*dt;
        if(i>0 & i < Nout-1){
            (*Pxx)[i] *= 2;
        }
    }
    /*
    for(int i = 0; i < Nout; i++){
        printf("Pxx[%i]: %1.3f\n", i, Pxx[i]);
    }
     */
    
    *f = malloc(Nout * sizeof(double));
    for(int i = 0; i < Nout; i++){
        (*f)[i] = (double)i*df;
    }
    /*
    for(int i = 0; i < Nout; i++){
        printf("f[%i]: %1.3f\n", i, (*f)[i]);
    }
     */
    
    free(P);
    free(F);
    free(tw);
    free(xw);
    
    return Nout;
}

double SP_Summaries_welch_rect(const double y[], const int size, const char what[])
{
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    // rectangular window for Welch-spectrum
    double * window = malloc(size * sizeof(double));
    for(int i = 0; i < size; i++){
        window[i] = 1;
    }
    
    double Fs = 1.0; // sampling frequency
    int N = nextpow2(size);
    
    double * S;
    double * f;
    
    // compute Welch-power
    int nWelch = welch(y, size, N, Fs, window, size, &S, &f);
    free(window);
    
    // angualr frequency and spectrum on that
    double * w = malloc(nWelch * sizeof(double));
    double * Sw = malloc(nWelch * sizeof(double));
    
    double PI = 3.14159265359;
    for(int i = 0; i < nWelch; i++){
        w[i] = 2*PI*f[i];
        Sw[i] = S[i]/(2*PI);
        //printf("w[%i]=%1.3f, Sw[%i]=%1.3f\n", i, w[i], i, Sw[i]);
        if(isinf(Sw[i]) | isinf(-Sw[i])){
            return 0;
        }
    }
    
    double dw = w[1] - w[0];
    
    double * csS = malloc(nWelch * sizeof(double));
    cumsum(Sw, nWelch, csS);
    /*
    for(int i=0; i<nWelch; i++)
    {
        printf("csS[%i]=%1.3f\n", i, csS[i]);
    }
     */
    
    double output = 0;
    
    if(strcmp(what, "centroid") == 0){
        
        double csSThres = csS[nWelch-1]*0.5;
        double centroid = 0;
        for(int i = 0; i < nWelch; i ++){
            if(csS[i] > csSThres){
                centroid = w[i];
                break;
            }
        }
        
        output = centroid;
        
    }
    else if(strcmp(what, "area_5_1") == 0){
        double area_5_1 = 0;;
        for(int i=0; i<nWelch/5; i++){
            area_5_1 += Sw[i];
        }
        area_5_1 *= dw;
        
        output = area_5_1;
    }
    
    free(w);
    free(Sw);
    free(csS);
    free(f);
    free(S);
    
    return output;
    
    
}


double SP_Summaries_welch_rect_area_5_1(const double y[], const int size)
{
    return SP_Summaries_welch_rect(y, size, "area_5_1");
}
double SP_Summaries_welch_rect_centroid(const double y[], const int size)
{
    return SP_Summaries_welch_rect(y, size, "centroid");
    
}
