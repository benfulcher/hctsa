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
/*Author: Rainer Hegger Last modified: March 1st, 1998 */
void make_box(double *ser,long **box,long *list,unsigned long l,
		 unsigned int bs,unsigned int dim,unsigned int del,double eps)
{
  int i,x,y;
  int ib=bs-1;

  for (x=0;x<bs;x++)
    for (y=0;y<bs;y++)
      box[x][y] = -1;
  
  for (i=(dim-1)*del;i<l;i++) {
    x=(int)(ser[i-(dim-1)*del]/eps)&ib;
    y=(int)(ser[i]/eps)&ib;
    list[i]=box[x][y];
    box[x][y]=i;
  }
}

