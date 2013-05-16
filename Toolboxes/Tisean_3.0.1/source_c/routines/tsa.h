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
/*Author: Rainer Hegger Last modified: Sep 3, 1999 */

#ifndef _TSA_ROUTINES_H
#define _TSA_ROUTINES_H

#ifndef _TISEAN_CEC_H
#include "tisean_cec.h"
#endif

/* size of the string which reads the input data
   if your lines are longer than some 500 reals, increase the value
   */
#define INPUT_SIZE 1024

/* The possible names of the verbosity levels */
#define VER_INPUT 0x1
#define VER_USR1 0x2
#define VER_USR2 0x4
#define VER_USR3 0x8
#define VER_USR4 0x10
#define VER_USR5 0x20
#define VER_USR6 0x40
#define VER_FIRST_LINE 0x80

/* Uncomment the variable to get rid of the initial Version message */
/*#define OMIT_WHAT_I_DO*/

#define sqr(x) ((x)*(x))

#ifdef __cplusplus
extern "C" {
#endif

extern int scan_help(int,char**);
extern double *get_series(char *,unsigned long *,unsigned long,
		       unsigned int,unsigned int);
extern double **get_multi_series(char *,unsigned long *,unsigned long,
				 unsigned int *,char *,char,unsigned int);
extern void rescale_data(double *,unsigned long,double *,double *);
extern void variance(double *,unsigned long,double *,double *);
extern void make_box(double *,long **,long *,unsigned long,
			unsigned int,unsigned int,unsigned int,double);
extern unsigned long find_neighbors(double *,long **,long *,double *,
				    unsigned long,unsigned int,unsigned int,
				    unsigned int,double,unsigned long *);
extern char* search_datafile(int, char**,unsigned int*,unsigned int);
extern char* check_option(char**,int,int,int);
extern void  solvele(double**,double *,unsigned int);
extern void test_outfile(char*);
extern double** invert_matrix(double**,unsigned int);
extern unsigned long exclude_interval(unsigned long,long,long,
				      unsigned long*,unsigned long*);
extern void make_multi_box(double **,long **,long *,unsigned long,
			   unsigned int,unsigned int,unsigned int,
			   unsigned int,double);
  /*only used for nrlazy. Will be removed with nrlazy */
extern void make_multi_box2(double **,long **,long *,unsigned long,
			   unsigned int,unsigned int,unsigned int,
			   unsigned int,double);
extern unsigned long find_multi_neighbors(double **,long **,long *,double **,
					  unsigned long,unsigned int,
					  unsigned int,unsigned int,
					  unsigned int,double,unsigned long *);
extern unsigned int** make_multi_index(unsigned int,unsigned int,unsigned int);

extern void check_alloc(void *);
extern char* myfgets(char *,int *,FILE *,unsigned int);
extern void what_i_do(char *, char *);
extern double* rand_arb_dist(double *,unsigned long,unsigned long,
			     unsigned int,unsigned long);

/* routines from rand.c */
extern void rnd_init(unsigned long);
extern unsigned long rnd_long();
extern unsigned long rnd_1279();
extern unsigned long rnd69069();
extern double gaussian(double);

/* routines from eigen.c */
extern void eigen(double**,unsigned long,double*);

#ifdef __cplusplus
}
#endif

#endif
