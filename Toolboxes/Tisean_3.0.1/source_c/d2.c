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
/*Author: Rainer Hegger. Last modified May 10, 2000 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the correlation sum, -dimension and -entropy"

/* output is written every WHEN seconds */
#define WHEN 120
/* Size of the field for box assisted neighbour searching 
   (has to be a power of 2)*/
#define NMAX 256
/* Size of the box for the scramble routine */
#define SCBOX 4096

double **series;
long *scr;
char dimset=0,rescale_set=0,eps_min_set=0,eps_max_set=0;
char *FOUT=NULL;
double epsfactor,epsinv,lneps,lnfac;
double EPSMAX=1.0,EPSMIN=1.e-3;
double min,interval;
int imax=NMAX-1,howoften1,imin;
long box[NMAX][NMAX],*list,boxc1[NMAX],*listc1;
unsigned long nmax;
double **found,*norm;
unsigned long MINDIST=0,MAXFOUND=1000;
unsigned long length=ULONG_MAX,exclude=0;
unsigned int DIM=1,EMBED=10,HOWOFTEN=100,DELAY=1;
unsigned int verbosity=0x1;
char *column=NULL;
char *infile=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"  Usage: %s [options]\n",progname);
  fprintf(stderr,"  Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l datapoints [default is whole file]\n");
  fprintf(stderr,"\t-x exclude # points [default 0]\n");
  fprintf(stderr,"\t-d delay  [default 1]\n");
  fprintf(stderr,"\t-M # of components, max. embedding dim. [default 1,10]\n");
  fprintf(stderr,"\t-c columns [default 1,...,# of components]\n");
  fprintf(stderr,"\t-t theiler-window [default 0]\n");
  fprintf(stderr,"\t-R max-epsilon "
	  "[default: max data interval]\n");
  fprintf(stderr,"\t-r min-epsilon [default: (max data interval)/1000]\n");
  fprintf(stderr,"\t-# #-of-epsilons [default 100]\n");
  fprintf(stderr,"\t-N max-#-of-pairs (0 means all) [default 1000]\n");
  fprintf(stderr,"\t-E use rescaled data [default: not rescaled]\n");
  fprintf(stderr," \t-o outfiles"
	  " [without exts.! default datafile[.d2][.h2][.stat][.c2]]\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
	  "2='+ output message each time output is done\n");
  
  fprintf(stderr,"\t-h show these options\n");
  fprintf(stderr,"\n");
  exit(0);
}

void scan_options(int n,char **argv)
{
  char *out;
  
  if ((out=check_option(argv,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(argv,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(argv,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(argv,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(argv,n,'M','2')) != NULL) {
    sscanf(out,"%u,%u",&DIM,&EMBED);
    dimset=1;
  }
  if ((out=check_option(argv,n,'t','u')) != NULL)
    sscanf(out,"%lu",&MINDIST);
  if ((out=check_option(argv,n,'R','f')) != NULL) {
    sscanf(out,"%lf",&EPSMAX);
    eps_max_set=1;
  }
  if ((out=check_option(argv,n,'r','f')) != NULL) {
    sscanf(out,"%lf",&EPSMIN);
    eps_min_set=1;
  }
  if ((out=check_option(argv,n,'#','u')) != NULL)
    sscanf(out,"%u",&HOWOFTEN);
  if ((out=check_option(argv,n,'N','u')) != NULL) {
    sscanf(out,"%lu",&MAXFOUND);
    if (MAXFOUND == 0)
      MAXFOUND=ULONG_MAX;
  }
  if ((out=check_option(argv,n,'E','n')) != NULL)
    rescale_set=1;
  if ((out=check_option(argv,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      FOUT=out;
}
      
void scramble(void)
{
  long i,j,k,m;
  unsigned long rnd,rndf,hlength,allscr=0;
  long *scfound,*scnhelp,scnfound;
  long scbox[SCBOX],lswap,element,scbox1=SCBOX-1;
  double *rz,*schelp,sceps=(double)SCBOX-0.001,swap;
  
  hlength=length-(EMBED-1)*DELAY;

  if (sizeof(long) == 8) {
    rndf=13*13*13*13;
    rndf=rndf*rndf*rndf*13;
    rnd=0x849178L;
  }
  else {
    rndf=69069;
    rnd=0x234571L;
  }
  for (i=0;i<1000;i++)
    rnd=rnd*rndf+1;

  check_alloc(rz=(double*)malloc(sizeof(double)*hlength));
  check_alloc(scfound=(long*)malloc(sizeof(long)*hlength));
  check_alloc(scnhelp=(long*)malloc(sizeof(long)*hlength));
  check_alloc(schelp=(double*)malloc(sizeof(double)*hlength));

  for (i=0;i<hlength;i++)
    rz[i]=(double)(rnd=rnd*rndf+1)/ULONG_MAX;
  
  for (i=0;i<SCBOX;i++)
    scbox[i]= -1;
  for (i=0;i<hlength;i++) {
    m=(int)(rz[i]*sceps)&scbox1;
    scfound[i]=scbox[m];
    scbox[m]=i;
  }
  for (i=0;i<SCBOX;i++) {
    scnfound=0;
    element=scbox[i];
    while(element != -1) {
      scnhelp[scnfound]=element;
      schelp[scnfound++]=rz[element];
      element=scfound[element];
    }
    
    for (j=0;j<scnfound-1;j++)
      for (k=j+1;k<scnfound;k++)
	if (schelp[k] < schelp[j]) {
	  swap=schelp[k];
	  schelp[k]=schelp[j];
	  schelp[j]=swap;
	  lswap=scnhelp[k];
	  scnhelp[k]=scnhelp[j];
	  scnhelp[j]=lswap;
	}
    for (j=0;j<scnfound;j++)
      scr[allscr+j]=scnhelp[j];
    allscr += scnfound;
  }

  free(rz);
  free(scfound);
  free(schelp);
}

void make_c2_dim(int n)
{
  char small;
  long i,j,k,x,y,i1,i2,j1,element,n1,maxi,count,hi;
  double *hs,max,dx;
  
  check_alloc(hs=(double*)malloc(sizeof(double)*EMBED*DIM));
  n1=scr[n];

  count=0;
  for (i1=0;i1<EMBED;i1++) {
    i2=i1*DELAY;
    for (j=0;j<DIM;j++)
      hs[count++]=series[j][n1+i2];
  }

  x=(int)(hs[0]*epsinv)&imax;
  y=(int)(hs[1]*epsinv)&imax;
  
  for (i1=x-1;i1<=x+1;i1++) {
    i2=i1&imax;
    for (j1=y-1;j1<=y+1;j1++) {
      element=box[i2][j1&imax];
      while (element != -1) {
	if (labs((long)(element-n1)) > MINDIST) {
	  count=0;
	  max=0.0;
	  maxi=howoften1;
	  small=0;
	  for (i=0;i<EMBED;i++) {
	    hi=i*DELAY;
	    for (j=0;j<DIM;j++) {
	      dx=fabs(hs[count]-series[j][element+hi]);
	      if (dx <= EPSMAX) {
		if (dx > max) {
		  max=dx;
		  if (max < EPSMIN) {
		    maxi=howoften1;
		  }
		  else {
		    maxi=(lneps-log(max))/lnfac;
		  }
		}
		if (count > 0)
		  for (k=imin;k<=maxi;k++)
		    found[count][k] += 1.0;
	      }
	      else {
		small=1;
		break;
	      }
	      count++;
	    }
	    if (small)
	      break;
	  }
	}
	element=list[element];
      }
    }
  }

  free(hs);
}

void make_c2_1(int n)
{
  int i,x,i1,maxi;
  long element,n1;
  double hs,max;
  
  n1=scr[n];
  hs=series[0][n1];
  
  x=(int)(hs*epsinv)&imax;
  
  for (i1=x-1;i1<=x+1;i1++) {
    element=boxc1[i1&imax];
    while (element != -1) {
      if (labs(element-n1) > MINDIST) {
	max=fabs(hs-series[0][element]);
	if (max <= EPSMAX) {
	  if (max < EPSMIN)
	    maxi=howoften1;
	  else
	    maxi=(lneps-log(max))/lnfac;
	  for (i=imin;i<=maxi;i++)
	    found[0][i] += 1.0;
	}
      }
      element=listc1[element];
    }
  }
}

int main(int argc,char **argv)
{
  char smaller,stdi=0;
  FILE *fout,*fstat;
  char *outd1,*outc1,*outh1,*outstat;
  int maxembed;
  long i1,j1,x,y,sn,n,i,j,n1,n2;
  long *oscr;
  long lnorm;
  double eps,*epsm,EPSMAX1,maxinterval;
  time_t mytime,lasttime;
  
  if (scan_help(argc,argv)) 
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;
  
  if (FOUT == NULL) {
    if (!stdi) {
      check_alloc(FOUT=calloc(strlen(infile)+1,(size_t)1));
      strcpy(FOUT,infile);
    }
    else {
      check_alloc(FOUT=calloc((size_t)6,(size_t)1));
      strcpy(FOUT,"stdin");
    }
  }
  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&DIM,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&DIM,column,
				      dimset,verbosity);

  if (rescale_set) {
    for (i=0;i<DIM;i++)
      rescale_data(series[i],length,&min,&interval);
    maxinterval=1.0;
  }
  else {
    maxinterval=0.0;
    for (i=0;i<DIM;i++) {
      min=interval=series[i][0];
      for (j=1;j<length;j++) {
	if (min > series[i][j])
	  min=series[i][j];
	if (interval < series[i][j])
	  interval=series[i][j];
      }
      interval -= min;
      if (interval > maxinterval)
	maxinterval=interval;
    }
  }
  if (!eps_max_set)
    EPSMAX *= maxinterval;
  if (!eps_min_set)
    EPSMIN *= maxinterval;
  EPSMAX=(fabs(EPSMAX)<maxinterval) ? fabs(EPSMAX) : maxinterval;
  EPSMIN=(fabs(EPSMIN)<EPSMAX) ? fabs(EPSMIN) : EPSMAX/2.;
  EPSMAX1=EPSMAX;

  howoften1=HOWOFTEN-1;
  maxembed=DIM*EMBED-1;

  check_alloc(outd1=(char*)calloc(strlen(FOUT)+4,(size_t)1));
  check_alloc(outc1=(char*)calloc(strlen(FOUT)+4,(size_t)1));
  check_alloc(outh1=(char*)calloc(strlen(FOUT)+4,(size_t)1));
  check_alloc(outstat=(char*)calloc(strlen(FOUT)+6,(size_t)1));
  strcpy(outd1,FOUT);
  strcpy(outc1,FOUT);
  strcpy(outh1,FOUT);
  strcpy(outstat,FOUT);
  strcat(outd1,".d2");
  strcat(outc1,".c2");
  strcat(outh1,".h2");
  strcat(outstat,".stat");
  test_outfile(outd1);
  test_outfile(outc1);
  test_outfile(outh1);
  test_outfile(outstat);

  check_alloc(list=(long*)malloc(length*sizeof(long)));
  check_alloc(listc1=(long*)malloc(length*sizeof(long)));
  if ((long)(length-(EMBED-1)*DELAY) <= 0) {
    fprintf(stderr,"Embedding dimension and delay are too large.\n"
	    "The delay vector would be longer than the whole series."
	    " Exiting\n");
    exit(VECTOR_TOO_LARGE_FOR_LENGTH);
  }
  check_alloc(scr=(long*)malloc(sizeof(long)*(length-(EMBED-1)*DELAY)));
  check_alloc(oscr=(long*)malloc(sizeof(long)*(length-(EMBED-1)*DELAY)));
  check_alloc(found=(double**)malloc(DIM*EMBED*sizeof(double*)));
  for (i=0;i<EMBED*DIM;i++)
    check_alloc(found[i]=(double*)malloc(HOWOFTEN*sizeof(double)));
  check_alloc(norm=(double*)malloc(HOWOFTEN*sizeof(double)));
  check_alloc(epsm=(double*)malloc(HOWOFTEN*sizeof(double)));
  
  epsinv=1.0/EPSMAX;
  epsfactor=pow(EPSMAX/EPSMIN,1.0/(double)howoften1);
  lneps=log(EPSMAX);
  lnfac=log(epsfactor);

  epsm[0]=EPSMAX;
  norm[0]=0.0;
  for (i=1;i<HOWOFTEN;i++) {
    norm[i]=0.0;
    epsm[i]=epsm[i-1]/epsfactor;
  }
  imin=0;

  scramble();
  for (i=0;i<(length-(EMBED-1)*DELAY);i++)
    oscr[scr[i]]=i;

  for (i=0;i<DIM*EMBED;i++)
    for (j=0;j<HOWOFTEN;j++)
      found[i][j]=0.0;
  
  nmax=length-DELAY*(EMBED-1);

  for (i1=0;i1<NMAX;i1++) {
    boxc1[i1]= -1;
    for (j1=0;j1<NMAX;j1++)
      box[i1][j1]= -1;
  }
  time(&lasttime);
  lnorm=0;
  
  for (n=1;n<nmax;n++) {
    smaller=0;
    sn=scr[n-1];
    if (DIM > 1) {
      x=(long)(series[0][sn]*epsinv)&imax;
      y=(long)(series[1][sn]*epsinv)&imax;
    }
    else {
      x=(long)(series[0][sn]*epsinv)&imax;
      y=(long)(series[0][sn+DELAY]*epsinv)&imax;
    }
    list[sn]=box[x][y];
    box[x][y]=sn;
    listc1[sn]=boxc1[x];
    boxc1[x]=sn;

    i=imin;
    while (found[maxembed][i] >= MAXFOUND) {
      smaller=1;
      if (++i > howoften1)
	break;
    }
    if (smaller) {
      imin=i;
      if (imin <= howoften1) {
	EPSMAX=epsm[imin];
	epsinv=1.0/EPSMAX;
	for (i1=0;i1<NMAX;i1++) {
	  boxc1[i1]= -1;
	  for (j1=0;j1<NMAX;j1++)
	    box[i1][j1]= -1;
	}
	for (i1=0;i1<n;i1++) {
	  sn=scr[i1];
	  if (DIM > 1) {
	    x=(long)(series[0][sn]*epsinv)&imax;
	    y=(long)(series[1][sn]*epsinv)&imax;
	  }
	  else {
	    x=(long)(series[0][sn]*epsinv)&imax;
	    y=(long)(series[0][sn+DELAY]*epsinv)&imax;
	  }
	  list[sn]=box[x][y];
	  box[x][y]=sn;
	  listc1[sn]=boxc1[x];
	  boxc1[x]=sn;
	}
      }
    }

    if (imin <= howoften1) {
      lnorm=n;
      if (MINDIST > 0) {
	sn=scr[n];
	n1=(sn-(long)MINDIST>=0)?sn-(long)MINDIST:0;
	n2=(sn+MINDIST<length-(EMBED-1)*DELAY)?sn+MINDIST:
	  length-(EMBED-1)*DELAY-1;
	for (i1=n1;i1<=n2;i1++)
	  if ((oscr[i1] < n))
	    lnorm--;
      }
      
      if (EMBED*DIM > 1)
	make_c2_dim(n);
      make_c2_1(n);
      for (i=imin;i<HOWOFTEN;i++)
	norm[i] += (double)(lnorm);
    }
    
    if (((time(&mytime)-lasttime) > WHEN) || (n == (nmax-1)) || 
	(imin > howoften1)) {
      time(&lasttime);
      fstat=fopen(outstat,"w");
      if (verbosity&VER_USR1)
	fprintf(stderr,"Opened %s for writing\n",outstat);
      fprintf(fstat,"Center points treated so far= %ld\n",n);
      fprintf(fstat,"Maximal epsilon in the moment= %e\n",epsm[imin]);
      fclose(fstat);
      fout=fopen(outc1,"w");
      if (verbosity&VER_USR1)
	fprintf(stderr,"Opened %s for writing\n",outc1);
      fprintf(fout,"#center= %ld\n",n);
      for (i=0;i<EMBED*DIM;i++) {
	fprintf(fout,"#dim= %ld\n",i+1);
	eps=EPSMAX1*epsfactor;
	for (j=0;j<HOWOFTEN;j++) {
	  eps /= epsfactor;
	  if (norm[j] > 0.0)
	    fprintf(fout,"%e %e\n",eps,found[i][j]/norm[j]);
	}
	fprintf(fout,"\n\n");
      }
      fclose(fout);
      fout=fopen(outh1,"w");
      if (verbosity&VER_USR1)
	fprintf(stderr,"Opened %s for writing\n",outh1);
      fprintf(fout,"#center= %ld\n",n);
      fprintf(fout,"#dim= 1\n");
      eps=EPSMAX1*epsfactor;
      for (j=0;j<HOWOFTEN;j++) {
	eps /= epsfactor;
	if (found[0][j] > 0.0)
	  fprintf(fout,"%e %e\n",eps,-log(found[0][j]/norm[j]));
      }
      fprintf(fout,"\n\n");
      for (i=1;i<DIM*EMBED;i++) {
	fprintf(fout,"#dim= %ld\n",i+1);
	eps=EPSMAX1*epsfactor;
	for (j=0;j<HOWOFTEN;j++) {
	  eps /= epsfactor;
	  if ((found[i-1][j] > 0.0) && (found[i][j] > 0.0))
	    fprintf(fout,"%e %e\n",eps,log(found[i-1][j]/found[i][j]));
	}
	fprintf(fout,"\n\n");
      }
      fclose(fout);
      fout=fopen(outd1,"w");
      if (verbosity&VER_USR1)
	fprintf(stderr,"Opened %s for writing\n",outd1);
      fprintf(fout,"#center= %ld\n",n);
      for (i=0;i<DIM*EMBED;i++) {
	fprintf(fout,"#dim= %ld\n",i+1);
	eps=EPSMAX1;
	for (j=1;j<HOWOFTEN;j++) {
	  eps /= epsfactor;
	  if ((found[i][j] > 0.0) && (found[i][j-1] > 0.0))
	    fprintf(fout,"%e %e\n",eps,log(found[i][j-1]/found[i][j]
					   /norm[j-1]*norm[j])/lnfac);
	}
	fprintf(fout,"\n\n");
      }
      fclose(fout);
      if (imin > howoften1)
	exit(0);
    }
  }

  if (infile != NULL)
    free(infile);
  free(outd1);
  free(outh1);
  free(outc1);
  free(outstat);
  free(list);
  free(listc1);
  free(scr);
  free(oscr);
  free(norm);
  free(epsm);
  for (i=0;i<EMBED*DIM;i++)
    free(found[i]);
  free(found);
  for (i=0;i<DIM;i++)
    free(series[i]);
  free(series);

  return 0;
}
