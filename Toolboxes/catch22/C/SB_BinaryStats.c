//
//  SB_BinaryStats.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include "SB_BinaryStats.h"
#include "stats.h"

double SB_BinaryStats_diff_longstretch0(const double y[], const int size){
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    // binarize
    int * yBin = malloc((size-1) * sizeof(int));
    for(int i = 0; i < size-1; i++){
        
        double diffTemp = y[i+1] - y[i];
        yBin[i] = diffTemp < 0 ? 0 : 1;
        
        /*
        if( i < 300)
            printf("%i, y[i+1]=%1.3f, y[i]=%1.3f, yBin[i]=%i\n", i, y[i+1], y[i], yBin[i]);
         */
        
    }
    
    int maxstretch0 = 0;
    int last1 = 0;
    for(int i = 0; i < size-1; i++){
        if(yBin[i] == 1 | i == size-2){
            double stretch0 = i - last1;
            if(stretch0 > maxstretch0){
                maxstretch0 = stretch0;
            }
            last1 = i;
        }
    }
    
    free(yBin);
    
    return maxstretch0;
}

double SB_BinaryStats_mean_longstretch1(const double y[], const int size){
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    // binarize
    int * yBin = malloc((size-1) * sizeof(int));
    double yMean = mean(y, size);
    for(int i = 0; i < size-1; i++){
        
        yBin[i] = (y[i] - yMean <= 0) ? 0 : 1;
        //printf("yBin[%i]=%i\n", i, yBin[i]);
        
    }
    
    int maxstretch1 = 0;
    int last1 = 0;
    for(int i = 0; i < size-1; i++){
        if(yBin[i] == 0 | i == size-2){
            double stretch1 = i - last1;
            if(stretch1 > maxstretch1){
                maxstretch1 = stretch1;
            }
            last1 = i;
        }
        
    }
    
    free(yBin);
    
    return maxstretch1;
}
