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
/* Author: Rainer Hegger */
/* Changes:
   10/12/05: First version
*/
/* Comments: Parameters are no. of components of the ts, embedding
             dimension and optionally the delay
             return: [0][i] components, [1][i] delay
*/

#include <stdlib.h>

extern void check_alloc(void *);

unsigned int **make_multi_index(unsigned int comps,unsigned int emb,
				unsigned int del)
{
  unsigned long i,alldim;
  unsigned int **mmi;

  alldim=comps*emb;
  check_alloc(mmi=(unsigned int**)malloc(sizeof(unsigned int*)*2));
  for (i=0;i<2;i++)
    check_alloc(mmi[i]=(unsigned int*)malloc(sizeof(unsigned int)*alldim));

  for (i=0;i<alldim;i++) {
    mmi[0][i]=i%comps;
    mmi[1][i]=(i/comps)*del;
  }

  return mmi;
}
