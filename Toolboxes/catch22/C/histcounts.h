//
//  histcounts.h
//  C_polished
//
//  Created by Carl Henning Lubba on 19/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#ifndef histcounts_h
#define histcounts_h

#include <stdlib.h>
#include <float.h>
#include <stdio.h>

extern int num_bins_auto(const double y[], const int size);
extern int histcounts(const double y[], const int size, int nBins, int ** binCounts, double ** binEdges);
extern int histcounts_preallocated(const double y[], const int size, int nBins, int * binCounts, double * binEdges);
extern int * histcount_edges(const double y[], const int size, const double binEdges[], const int nEdges);
extern int * histbinassign(const double y[], const int size, const double binEdges[], const int nEdges);

#endif /* histcounts_h */
