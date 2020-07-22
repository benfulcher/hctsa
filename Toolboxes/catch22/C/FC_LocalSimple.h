#ifndef FC_LOCALSIMPLE_H
#define FC_LOCALSIMPLE_H
#include <math.h>
#include <string.h>
#include "stats.h"
#include "CO_AutoCorr.h"

extern double fc_local_simple(const double y[], const int size, const int train_length);
extern double FC_LocalSimple_mean_taures(const double y[], const int size, const int train_length);
extern double FC_LocalSimple_lfit_taures(const double y[], const int size);
extern double FC_LocalSimple_mean_tauresrat(const double y[], const int size, const int train_length);
extern double FC_LocalSimple_mean1_tauresrat(const double y[], const int size);
extern double FC_LocalSimple_mean_stderr(const double y[], const int size, const int train_length);
extern double FC_LocalSimple_mean3_stderr(const double y[], const int size);

#endif
