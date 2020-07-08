//
//  histcounts.c
//  C_polished
//
//  Created by Carl Henning Lubba on 19/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include <stdlib.h>
#include <float.h>

#include "stats.h"
#include "histcounts.h"

int num_bins_auto(const double y[], const int size){
    
    double maxVal = max_(y, size);
    double minVal = min_(y, size);
    
    if (stddev(y, size) < 0.001){
        return 0;
    }
    
    return ceil((maxVal-minVal)/(3.5*stddev(y, size)/pow(size, 1/3.)));
    
}

int histcounts_preallocated(const double y[], const int size, int nBins, int * binCounts, double * binEdges)
{
    
    int i = 0;
    
    // check min and max of input array
    double minVal = DBL_MAX, maxVal=-DBL_MAX;
    for(int i = 0; i < size; i++)
    {
        // printf("histcountInput %i: %1.3f\n", i, y[i]);
        
        if (y[i] < minVal)
        {
            minVal = y[i];
        }
        if (y[i] > maxVal)
        {
            maxVal = y[i];
        }
    }
    
    // and derive bin width from it
    double binStep = (maxVal - minVal)/nBins;
    
    // variable to store counted occurances in
    for(i = 0; i < nBins; i++)
    {
        binCounts[i] = 0;
    }
    
    for(i = 0; i < size; i++)
    {
        
        int binInd = (y[i]-minVal)/binStep;
        if(binInd < 0)
            binInd = 0;
        if(binInd >= nBins)
            binInd = nBins-1;
        //printf("histcounts, i=%i, binInd=%i, nBins=%i\n", i, binInd, nBins);
        binCounts[binInd] += 1;
        
    }
    
    for(i = 0; i < nBins+1; i++)
    {
        binEdges[i] = i * binStep + minVal;
    }
    
    /*
     // debug
     for(i=0;i<nBins;i++)
     {
     printf("%i: count %i, edge %1.3f\n", i, binCounts[i], binEdges[i]);
     }
     */
    
    return 0;
    
}

int histcounts(const double y[], const int size, int nBins, int ** binCounts, double ** binEdges)
{

    int i = 0;
    
    // check min and max of input array
    double minVal = DBL_MAX, maxVal=-DBL_MAX;
    for(int i = 0; i < size; i++)
    {
        // printf("histcountInput %i: %1.3f\n", i, y[i]);
        
        if (y[i] < minVal)
        {
            minVal = y[i];
        }
        if (y[i] > maxVal)
        {
            maxVal = y[i];
        }
    }
    
    // if no number of bins given, choose spaces automatically
    if (nBins <= 0){
        nBins = ceil((maxVal-minVal)/(3.5*stddev(y, size)/pow(size, 1/3.)));
    }
    
    // and derive bin width from it
    double binStep = (maxVal - minVal)/nBins;
    
    // variable to store counted occurances in
    *binCounts = malloc(nBins * sizeof(int));
    for(i = 0; i < nBins; i++)
    {
        (*binCounts)[i] = 0;
    }
    
    for(i = 0; i < size; i++)
    {
        
        int binInd = (y[i]-minVal)/binStep;
        if(binInd < 0)
            binInd = 0;
        if(binInd >= nBins)
            binInd = nBins-1;
        (*binCounts)[binInd] += 1;
        
    }
    
    *binEdges = malloc((nBins+1) * sizeof(double));
    for(i = 0; i < nBins+1; i++)
    {
        (*binEdges)[i] = i * binStep + minVal;
    }
   
    /*
    // debug
    for(i=0;i<nBins;i++)
    {
        printf("%i: count %i, edge %1.3f\n", i, binCounts[i], binEdges[i]);
    }
    */
    
    return nBins;
    
}

int * histbinassign(const double y[], const int size, const double binEdges[], const int nEdges)
{
    
    
    // variable to store counted occurances in
    int * binIdentity = malloc(size * sizeof(int));
    for(int i = 0; i < size; i++)
    {
        // if not in any bin -> 0
        binIdentity[i] = 0;
        
        // go through bin edges
        for(int j = 0; j < nEdges; j++){
            if(y[i] < binEdges[j]){
                binIdentity[i] = j;
                break;
            }
        }
    }
    
    return binIdentity;
    
}

int * histcount_edges(const double y[], const int size, const double binEdges[], const int nEdges)
{
    
    
    int * histcounts = malloc(nEdges * sizeof(int));
    for(int i = 0; i < nEdges; i++){
        histcounts[i] = 0;
    }
    
    for(int i = 0; i < size; i++)
    {
        // go through bin edges
        for(int j = 0; j < nEdges; j++){
            if(y[i] <= binEdges[j]){
                histcounts[j] += 1;
                break;
            }
        }
    }
    
    return histcounts;
    
}
