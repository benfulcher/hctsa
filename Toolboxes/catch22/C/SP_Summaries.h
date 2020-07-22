//
//  SP_Summaries.h
//  
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#ifndef SP_Summaries_h
#define SP_Summaries_h

#include <stdio.h>

extern double SP_Summaries_welch_rect(const double y[], const int size, const char what[]);
extern double SP_Summaries_welch_rect_area_5_1(const double y[], const int size);
extern double SP_Summaries_welch_rect_centroid(const double y[], const int size);

#endif /* SP_Summaries_h */
