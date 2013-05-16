/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Rainer Hegger. Last modified: Feb 12, 2006 */
/* Changes:
   Sep 5, 2004 : add extern definition of check_alloc
   Feb 12, 2006: add was_set to avoid multiple initialisations
*/

#define __RANDOM

#ifndef	_STDLIB_H
#include <stdlib.h>
#endif

#ifndef _LIMITS_H
#include <limits.h>
#endif

#ifndef _MATH_H
#include <math.h>
#endif


#ifndef M_PI
#define M_PI 3.1415926535897932385E0
#endif

extern void check_alloc(void*);

static unsigned long *rnd_array,rnd69,*rnd1279,factor;
static unsigned long *nexti,rndtime,rndtime1,rndtime2,rndtime3,*next1279;
static unsigned long t1279,t1279_1,t1279_2,t1279_3;
static double lo_limit;
static char rnd_init_was_set=0;

void rnd_init(unsigned long iseed)
{
  int i;
  unsigned long z,index;

  if (rnd_init_was_set == 1)
    return ;

  rnd_init_was_set=1;

  if (sizeof(long) == 8) {
    factor=13*13*13*13;
    factor=factor*factor*factor*13;
  }
  else
    factor=69069;
  lo_limit=(double)ULONG_MAX;

  check_alloc(rnd_array=(unsigned long *)malloc(9689*sizeof(unsigned long)));
  check_alloc(nexti=(unsigned long *)malloc(9689*sizeof(long)));
  check_alloc(rnd1279=(unsigned long *)malloc(1279*sizeof(unsigned long)));
  check_alloc(next1279=(unsigned long *)malloc(1279*sizeof(long)));

  rnd_array[0]=rnd1279[0]=iseed;
  rnd69=iseed;
  index=iseed;
  nexti[0]=next1279[0]=1;

  for (i=1;i<9689;i++) {
    rnd_array[i]=factor*rnd_array[i-1]+1;
    nexti[i]=i+1;
  }

  for (i=1;i<1279;i++) {
    rnd1279[i]=factor*rnd1279[i-1]+1;
    next1279[i]=i+1;
  }
  nexti[9688]=next1279[1278]=0;

  for (i=1;i<2000;i++) {
    index=factor*index+1;
    z=rnd1279[((index>>10)%1279)];
    z=(z<<10)+(z>>10);
    index=factor*index+1;
    rnd1279[((index>>10)%1279)] += z;
  }

  nexti[9688]=next1279[1278]=0;
  rndtime=9688;
  rndtime1=9688-157;
  rndtime2=9688-314;
  rndtime3=9688-471;
  t1279=1278;
  t1279_1=1278-216;
  t1279_2=1278-299;
  t1279_3=1278-598;
}

unsigned long rnd_long(void)
{
  rndtime=nexti[rndtime];
  rndtime1=nexti[rndtime1];
  rndtime2=nexti[rndtime2];
  rndtime3=nexti[rndtime3];
  rnd_array[rndtime] ^= rnd_array[rndtime1]
    ^rnd_array[rndtime2]^rnd_array[rndtime3];

  return rnd_array[rndtime];
}

unsigned long rnd_1279(void)
{
  t1279=next1279[t1279];
  t1279_1=next1279[t1279_1];
  t1279_2=next1279[t1279_2];
  t1279_3=next1279[t1279_3];

  rnd1279[t1279] += (rnd1279[t1279_1] + rnd1279[t1279_2] 
		     + rnd1279[t1279_3]);
  return rnd1279[t1279];
}

unsigned long rnd69069(void)
{
  return (rnd69=rnd69*factor+1);
}

double gaussian(double sigma)
{
  static unsigned long gausscount=0;
  double x,r,u,phi;
  static double y;
  
  if (!(gausscount++ & 0x1)) {
    phi=2.0*M_PI*(double)rnd_1279()/lo_limit;
    u=(double)rnd_1279()/lo_limit;
    r=sqrt(-2.0*sigma*sigma*log(u));
    x=r*cos(phi);
    y=r*sin(phi);

    return x;
  }
  else
    return y;
}
