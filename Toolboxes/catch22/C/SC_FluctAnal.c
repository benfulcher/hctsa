#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include "stats.h"
#include "CO_AutoCorr.h"

double SC_FluctAnal_2_50_1_logi_prop_r1(const double y[], const int size, const int lag, const char how[])
{
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    // generate log spaced tau vector
    double linLow = log(5);
    double linHigh = log(size/2);
    
    int nTauSteps = 50;
    double tauStep = (linHigh - linLow) / (nTauSteps-1);
    
    int tau[50];
    for(int i = 0; i < nTauSteps; i++)
    {
        tau[i] = round(exp(linLow + i*tauStep));
    }
    
    // check for uniqueness, use ascending order
    int nTau = nTauSteps;
    for(int i = 0; i < nTauSteps-1; i++)
    {
        
        while (tau[i] == tau[i+1] && i < nTau-1)
        {
            for(int j = i+1; j < nTauSteps-1; j++)
            {
                tau[j] = tau[j+1];
            }
            // lost one
            nTau -= 1;
        }
    }
    
    // fewer than 12 points -> leave.
    if(nTau < 12){
        return 0;
    }
    
    int sizeCS = size/lag;
    double * yCS = malloc(sizeCS * sizeof(double));
    
    /*
    for(int i = 0; i < 50; i++)
    {
        printf("y[%i]=%1.3f\n", i, y[i]);
    }
     */
    
    // transform input vector to cumsum
    yCS[0] = y[0];
    for(int i = 0; i < sizeCS-1; i++)
    {
        yCS[i+1] = yCS[i] + y[(i+1)*lag];
        
        /*
        if(i<300)
            printf("yCS[%i]=%1.3f\n", i, yCS[i]);
         */
    }
    
    //for each value of tau, cut signal into snippets of length tau, detrend and
    
    // first generate a support for regression (detrending)
    double * xReg = malloc(tau[nTau-1] * sizeof * xReg);
    for(int i = 0; i < tau[nTau-1]; i++)
    {
        xReg[i] = i+1;
    }
    
    // iterate over taus, cut signal, detrend and save amplitude of remaining signal
    double * F = malloc(nTau * sizeof * F);
    for(int i = 0; i < nTau; i++)
    {
        int nBuffer = sizeCS/tau[i];
        double * buffer = malloc(tau[i] * sizeof * buffer);
        double m = 0.0, b = 0.0;
        
        //printf("tau[%i]=%i\n", i, tau[i]);
        
        F[i] = 0;
        for(int j = 0; j < nBuffer; j++)
        {
            
            //printf("%i th buffer\n", j);
            
            linreg(tau[i], xReg, yCS+j*tau[i], &m, &b);
            
            
            for(int k = 0; k < tau[i]; k++)
            {
                buffer[k] = yCS[j*tau[i]+k] - (m * (k+1) + b);
                //printf("buffer[%i]=%1.3f\n", k, buffer[k]);
            }
            
            if (strcmp(how, "rsrangefit") == 0) {
                F[i] += pow(max_(buffer, tau[i]) - min_(buffer, tau[i]), 2);
            }
            else if (strcmp(how, "dfa") == 0) {
                for(int k = 0; k<tau[i]; k++){
                    F[i] += buffer[k]*buffer[k];
                }
            }
            else{
                return 0.0;
            }
        }
        
        if (strcmp(how, "rsrangefit") == 0) {
            F[i] = sqrt(F[i]/nBuffer);
        }
        else if (strcmp(how, "dfa") == 0) {
            F[i] = sqrt(F[i]/(nBuffer*tau[i]));
        }
        //printf("F[%i]=%1.3f\n", i, F[i]);
        
        free(buffer);
        
    }
    
    double * logtt = malloc(nTau * sizeof * logtt);
    double * logFF = malloc(nTau * sizeof * logFF);
    int ntt = nTau;
    
    for (int i = 0; i < nTau; i++)
    {
        logtt[i] = log(tau[i]);
        logFF[i] = log(F[i]);
    }
    
    int minPoints = 6;
    int nsserr = (ntt - 2*minPoints + 1);
    double * sserr = malloc(nsserr * sizeof * sserr);
    double * buffer = malloc((ntt - minPoints + 1) * sizeof * buffer);
    for (int i = minPoints; i < ntt - minPoints + 1; i++)
    {
        // this could be done with less variables of course
        double m1 = 0.0, b1 = 0.0;
        double m2 = 0.0, b2 = 0.0;
        
        sserr[i - minPoints] = 0.0;
        
        linreg(i, logtt, logFF, &m1, &b1);
        linreg(ntt-i+1, logtt+i-1, logFF+i-1, &m2, &b2);
        
        for(int j = 0; j < i; j ++)
        {
            buffer[j] = logtt[j] * m1 + b1 - logFF[j];
        }
        
        sserr[i - minPoints] += norm_(buffer, i);
        
        for(int j = 0; j < ntt-i+1; j++)
        {
            buffer[j] = logtt[j+i-1] * m2 + b2 - logFF[j+i-1];
        }
        
        sserr[i - minPoints] += norm_(buffer, ntt-i+1);
        
    }
    
    double firstMinInd = 0.0;
    double minimum = min_(sserr, nsserr);
    for(int i = 0; i < nsserr; i++)
    {
        if(sserr[i] == minimum)
        {
            firstMinInd = i + minPoints - 1;
            break;
        }
    }
    
    free(yCS); // new
    
    free(xReg);
    free(F);
    free(logtt);
    free(logFF);
    free(sserr);
    free(buffer);
    
    return (firstMinInd+1)/ntt;
    
}

double SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(const double y[], const int size){
    return SC_FluctAnal_2_50_1_logi_prop_r1(y, size, 2, "dfa");
}

double SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(const double y[], const int size){
    return SC_FluctAnal_2_50_1_logi_prop_r1(y, size, 1, "rsrangefit");
}

/*
double SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(double y[], int size)
{
    // generate log spaced tau vector
    double linLow = log(5);
    double linHigh = log(size/2);
    
    int nTauSteps = 50;
    double tauStep = (linHigh - linLow) / (nTauSteps-1);
    
    int tau[50];
    for(int i = 0; i < nTauSteps; i++)
    {
        tau[i] = round(exp(linLow + i*tauStep));
    }
    
   // check for uniqueness, use ascending order
    int nTau = nTauSteps;
    for(int i = 0; i < nTauSteps; i++)
    {
        
        while (tau[i] == tau[i+1] && i < nTau-1)
        {
            for(int j = i+1; j < nTauSteps-1; j++)
            {
                tau[j] = tau[j+1];
            }
            // lost one
            nTau -= 1;
        }
    }
    
    // fewer than 8 points -> leave.
    if(nTau < 8){
        return 0;
    }
    
    // transform input vector to cumsum
    for(int i = 0; i < size-1; i++)
    {
        y[i+1] = y[i] + y[i+1];
    }
    
    //for each value of tau, cut signal into snippets of length tau, detrend and
    
    // first generate a support for regression (detrending)
    double * xReg = malloc(tau[nTau-1] * sizeof * xReg);
    for(int i = 0; i < tau[nTau-1]; i++)
    {
        xReg[i] = i+1;
    }
    
    // iterate over taus, cut signal, detrend and save amplitude of remaining signal
    double * F = malloc(nTau * sizeof * F);
    for(int i = 0; i < nTau; i++)
    {
        int nBuffer = size/tau[i];
        double * buffer = malloc(tau[i] * sizeof * buffer);
        double m = 0.0, b = 0.0;
        
        F[i] = 0;
        for(int j = 0; j < nBuffer; j++)
        {
            
            linreg(tau[i], xReg, y+j*tau[i], &m, &b);
            
            for(int k = 0; k < tau[i]; k++)
            {
                buffer[k] = y[j*tau[i]+k] - (m * (k+1) + b);
            }
            
            F[i] += pow(max(buffer, tau[i]) - min(buffer, tau[i]), 2);
        }
        
        F[i] = sqrt(F[i]/nBuffer);
        
        free(buffer);
        
    }
    
    double * logtt = malloc(nTau * sizeof * logtt);
    double * logFF = malloc(nTau * sizeof * logFF);
    int ntt = nTau;
    
    for (int i = 0; i < nTau; i++)
    {
        logtt[i] = log(tau[i]);
        logFF[i] = log(F[i]);
    }
    
    int minPoints = 6;
    int nsserr = (ntt - 2*minPoints + 1);
    double * sserr = malloc(nsserr * sizeof * sserr);
    double * buffer = malloc((ntt - minPoints + 1) * sizeof * buffer);
    for (int i = minPoints; i < ntt - minPoints + 1; i++)
    {
        // this could be done with less variables of course
        double m1 = 0.0, b1 = 0.0;
        double m2 = 0.0, b2 = 0.0;
        
        sserr[i - minPoints] = 0.0;
        
        linreg(i, logtt, logFF, &m1, &b1);
        linreg(ntt-i+1, logtt+i-1, logFF+i-1, &m2, &b2);
        
        for(int j = 0; j < i; j ++)
        {
            buffer[j] = logtt[j] * m1 + b1 - logFF[j];
        }
        
        sserr[i - minPoints] += norm(buffer, i);
        
        for(int j = 0; j < ntt-i+1; j++)
        {
            buffer[j] = logtt[j+i-1] * m2 + b2 - logFF[j+i-1];
        }
        
        sserr[i - minPoints] += norm(buffer, ntt-i+1);
        
    }
    
    double firstMinInd = 0.0;
    double minimum = min(sserr, nsserr);
    for(int i = 0; i < nsserr; i++)
    {
        if(sserr[i] == minimum)
        {
            firstMinInd = i + minPoints - 1;
            break;
        }
    }
    
    free(xReg);
    free(F);
    free(logtt);
    free(logFF);
    free(sserr);
    free(buffer);
    
    return (firstMinInd+1)/ntt;
    
}
 */
