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
/*Author: Rainer Hegger Last modified: Feb 11, 2006 */
/*Changes:
  Feb 11, 2006: First version
*/
/*Comment:
  Creates a sequence of random numbers with arbitrary distribution
  Input Paramters: x=original data defining the distribution
                   nx= length of x
                   nc= number of random numbers to create
                   nb= number of bins for the distribution
		   iseed = seed for the random number generator
  Return: rand = array containing the nc random numbers
*/

#ifndef _STDLIB_H
#include <stdlib.h>
#endif

#ifndef _LIMITS_H
#include <limits.h>
#endif

#ifndef _TIME_H
#include <time.h>
#endif

extern void rescale_data(double*,unsigned long,double*,double*);
extern void check_alloc(void*);
extern unsigned long rnd_long(void);
extern void  rnd_init(unsigned long);

double *rand_arb_dist(double *x,unsigned long nx,unsigned long nc,
		      unsigned int nb,unsigned long iseed)
{
  double h,min,inter,*randarb,drnd,epsinv=1.0/(double)nb;
  unsigned long i,j,*box,hrnd,nall=nx+nb;

  rescale_data(x,nx,&min,&inter);

  check_alloc(box=(unsigned long*)malloc(sizeof(unsigned long)*nb));
  for (i=0;i<nb;i++)
    box[i]=1;

  for (i=0;i<nx;i++) {
    h=x[i];
    if (h >= 1.0) 
      h -= epsinv/2.0;
    j=(unsigned int)(h*nb);
    box[j]++;
  }
  for (i=1;i<nb;i++)
    box[i] += box[i-1];

  check_alloc(randarb=(double*)malloc(sizeof(double)*nc));

  if (iseed == 0)
    iseed=(unsigned long)time((time_t*)&iseed);

  rnd_init(iseed);
  for (i=0;i<1000;i++)
    rnd_long();

  for (i=0;i<nc;i++) {
    hrnd=rnd_long()%nall;
    for (j=0;j<nb;j++)
      if (box[j] >= hrnd)
	break;
    drnd=(double)rnd_long()/(double)ULONG_MAX*epsinv;
    randarb[i]=min+((double)j*epsinv+drnd)*inter;
  }

  free(box);

  return randarb;
}
