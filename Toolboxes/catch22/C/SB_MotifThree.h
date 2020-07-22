#ifndef SB_MOTIFTHREE_H
#define SB_MOTIFTHREE_H
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "SB_CoarseGrain.h"
#include "helper_functions.h"

extern double SB_MotifThree_quantile_hh(const double y[], const int size);
extern double * sb_motifthree(const double y[], int size, const char how[]);

#endif
