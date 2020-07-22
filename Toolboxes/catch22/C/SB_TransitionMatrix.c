//
//  SB_TransitionMatrix.c
//  
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#include "SB_TransitionMatrix.h"
#include "butterworth.h"
#include "CO_AutoCorr.h"
#include "SB_CoarseGrain.h"
#include "stats.h"

double SB_TransitionMatrix_3ac_sumdiagcov(const double y[], const int size)
{
    
    // NaN and const check
    int constant = 1;
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
        if(y[i] != y[0]){
            constant = 0;
        }
    }
    if (constant){
        return NAN;
    }
    
    const int numGroups = 3;
    
    int tau = co_firstzero(y, size, size);
    
    double * yFilt = malloc(size * sizeof(double));
    
    // sometimes causes problems in filt!!! needs fixing.
    /*
    if(tau > 1){
        butterworthFilter(y, size, 4, 0.8/tau, yFilt);
    }
    */
    
    for(int i = 0; i < size; i++){
        yFilt[i] = y[i];
    }
    
    /*
    for(int i = 0; i < size; i++){
        printf("yFilt[%i]=%1.4f\n", i, yFilt[i]);
    }
     */
    
    int nDown = (size-1)/tau+1;
    double * yDown = malloc(nDown * sizeof(double));
    
    for(int i = 0; i < nDown; i++){
        yDown[i] = yFilt[i*tau];
    }
    
    /*
    for(int i = 0; i < nDown; i++){
        printf("yDown[%i]=%1.4f\n", i, yDown[i]);
    }
     */
    
    
    // transfer to alphabet
    int * yCG = malloc(nDown * sizeof(double));
    sb_coarsegrain(yDown, nDown, "quantile", numGroups, yCG);
    
    /*
    for(int i = 0; i < nDown; i++){
        printf("yCG[%i]=%i\n", i, yCG[i]);
    }
     */
    
    
    double T[3][3];
    for(int i = 0; i < numGroups; i++){
        for(int j = 0; j < numGroups; j++){
            T[i][j] = 0;
        }
    }
    
    // more efficient way of doing the below 
    for(int j = 0; j < nDown-1; j++){
        T[yCG[j]-1][yCG[j+1]-1] += 1;
    }
    
    /*
    for(int i = 0; i < numGroups; i++){
        for(int j = 0; j < numGroups; j++){
            printf("%1.f, ", T[i][j]);
        }
        printf("\n");
    }
     */
    
    /*
    for(int i = 0; i < numGroups; i++){
        for(int j = 0; j < nDown-1; j++){
            if(yCG[j] == i+1){
                T[i][yCG[j+1]-1] += 1;
            }
        }
    }
     */
    
    for(int i = 0; i < numGroups; i++){
        for(int j = 0; j < numGroups; j++){
            T[i][j] /= (nDown-1);
            // printf("T(%i, %i) = %1.3f\n", i, j, T[i][j]);
            
        }
    }
    
    double column1[3] = {0};
    double column2[3] = {0};
    double column3[3] = {0};
    
    for(int i = 0; i < numGroups; i++){
        column1[i] = T[i][0];
        column2[i] = T[i][1];
        column3[i] = T[i][2];
        // printf("column3(%i) = %1.3f\n", i, column3[i]);
    }
    
    double *columns[3];
    columns[0] = &(column1[0]);
    columns[1] = &(column2[0]);
    columns[2] = &(column3[0]);
    
    
    double COV[3][3];
    double covTemp = 0;
    for(int i = 0; i < numGroups; i++){
        for(int j = i; j < numGroups; j++){
            
            covTemp = cov(columns[i], columns[j], 3);
            
            COV[i][j] = covTemp;
            COV[j][i] = covTemp;
            
            // printf("COV(%i , %i) = %1.3f\n", i, j, COV[i][j]);
        }
    }
    
    double sumdiagcov = 0;
    for(int i = 0; i < numGroups; i++){
        sumdiagcov += COV[i][i];
    }
    
    free(yFilt);
    free(yDown);
    free(yCG);
    
    return sumdiagcov;
    
    
}
