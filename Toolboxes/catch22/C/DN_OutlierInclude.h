#ifndef DN_OUTLIERINCLUDE_ABS_001
#define DN_OUTLIERINCLUDE_ABS_001
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "stats.h"

extern double DN_OutlierInclude_abs_001(const double y[], const int size);
extern double DN_OutlierInclude_np_001_mdrmd(const double y[], const int size, const int sign);
extern double DN_OutlierInclude_p_001_mdrmd(const double y[], const int size);
extern double DN_OutlierInclude_n_001_mdrmd(const double y[], const int size);

#endif
