/*Author: Rainer Hegger. Last modified Mar 20, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "tsa.h"

/* output is written every WHEN seconds */
#define WHEN 120
/* Size of the field for box assisted neighbour searching 
   (has to be a power of 2)*/
#define NMAX 256
/* Size of the box for the scramble routine */
#define SCBOX 4096

double *series;
long *scr;
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
unsigned int MAXEMBED=10,HOWOFTEN=100,DELAY=1;
unsigned int column=1;
char *infile=NULL;

void show_options(char *progname)
{
  fprintf(stderr,"Estimates the correlation sum, -dimension and -entropy\n\n");
  fprintf(stderr,"  Usage: %s [options]\n",progname);
  fprintf(stderr,"  Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l datapoints [default is whole file]\n");
  fprintf(stderr,"\t-x exclude # points [default 0]\n");
  fprintf(stderr,"\t-c column [default 1]\n");
  fprintf(stderr,"\t-d delay  [default 1]\n");
  fprintf(stderr,"\t-M maxdim [default 10]\n");
  fprintf(stderr,"\t-t theiler-window [default 0]\n");
  fprintf(stderr,"\t-E max-epsilon "
	  "[default 1 (data are normalized to [0:1])]\n");
  fprintf(stderr,"\t-e min-epsilon [default 1e-3]\n");
  fprintf(stderr,"\t-# #-of-epsilons [default 100]\n");
  fprintf(stderr,"\t-N max-#-of-pairs (0 means all) [default 1000]\n");
  fprintf(stderr," \t-o outfiles  "
	  " [without extensions! default datafile[.d2][.h2][.stat][.c2]]\n");
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
  if ((out=check_option(argv,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(argv,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(argv,n,'M','u')) != NULL)
    sscanf(out,"%u",&MAXEMBED);
  if ((out=check_option(argv,n,'t','u')) != NULL)
    sscanf(out,"%lu",&MINDIST);
  if ((out=check_option(argv,n,'E','f')) != NULL)
    sscanf(out,"%lf",&EPSMAX);
  if ((out=check_option(argv,n,'e','f')) != NULL)
    sscanf(out,"%lf",&EPSMIN);
  if ((out=check_option(argv,n,'#','u')) != NULL)
    sscanf(out,"%u",&HOWOFTEN);
  if ((out=check_option(argv,n,'N','u')) != NULL) {
    sscanf(out,"%lu",&MAXFOUND);
    if (MAXFOUND == 0)
      MAXFOUND=ULONG_MAX;
  }
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
  
  hlength=length-(MAXEMBED-1)*DELAY;

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
  long i,j,x,y,i1,i2,j1,element,n1,maxi;
  double *hs,max,dx;
  
  check_alloc(hs=(double*)malloc(sizeof(double)*MAXEMBED));
  n1=scr[n];
  for (i1=0;i1<MAXEMBED;i1++)
    hs[i1]=series[n1+i1*DELAY];
  
  x=(int)(hs[0]*epsinv)&imax;
  y=(int)(hs[1]*epsinv)&imax;
  
  for (i1=x-1;i1<=x+1;i1++) {
    i2=i1&imax;
    for (j1=y-1;j1<=y+1;j1++) {
      element=box[i2][j1&imax];
      while (element != -1) {
	if (labs((long)(element-n1)) > MINDIST) {
	  max=fabs(hs[0]-series[element]);
	  if (max <= EPSMAX) {
	    if (max < EPSMIN)
	      maxi=howoften1;
	    else
	      maxi=(lneps-log(max))/lnfac;
	    for (i=1;i<MAXEMBED;i++) {
	      dx=fabs(hs[i]-series[element+i*DELAY]);
	      if (dx <= EPSMAX) {
		if (dx > max) {
		  max=dx;
		  if (max < EPSMIN)
		    maxi=howoften1;
		  else
		    maxi=(lneps-log(max))/lnfac;
		}
		for (j=imin;j<=maxi;j++)
		  found[i][j] += 1.0;
	      }
	      else
		break;
	    }
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
  hs=series[n1];
  
  x=(int)(hs*epsinv)&imax;
  
  for (i1=x-1;i1<=x+1;i1++) {
    element=boxc1[i1&imax];
    while (element != -1) {
      if (labs(element-n1) > MINDIST) {
	max=fabs(hs-series[element]);
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
  double eps,*epsm,EPSMAX1;
  time_t mytime,lasttime;

  if (scan_help(argc,argv)) 
    show_options(argv[0]);
  
  scan_options(argc,argv);
  infile=search_datafile(argc,argv,&column);
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
  
  series=(double*)get_series(infile,&length,exclude,column);
  rescale_data(series,length,&min,&interval);
  
  EPSMAX=(fabs(EPSMAX)<1.0) ? fabs(EPSMAX) : 1.0;
  EPSMIN=(fabs(EPSMIN)<EPSMAX) ? fabs(EPSMIN) : EPSMAX/2.;
  EPSMAX1=EPSMAX;
  howoften1=HOWOFTEN-1;
  maxembed=MAXEMBED-1;

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
  check_alloc(scr=(long*)malloc(sizeof(long)*(length-(MAXEMBED-1)*DELAY)));
  check_alloc(oscr=(long*)malloc(sizeof(long)*(length-(MAXEMBED-1)*DELAY)));
  check_alloc(found=(double**)malloc(MAXEMBED*sizeof(double*)));
  for (i=0;i<MAXEMBED;i++)
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
  for (i=0;i<(length-(MAXEMBED-1)*DELAY);i++)
    oscr[scr[i]]=i;

  for (i=0;i<MAXEMBED;i++)
    for (j=0;j<HOWOFTEN;j++)
      found[i][j]=0.0;
  
  nmax=length-DELAY*(MAXEMBED-1);

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
    x=(long)(series[sn]*epsinv)&imax;
    y=(long)(series[sn+DELAY]*epsinv)&imax;
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
	  x=(long)(series[sn]*epsinv)&imax;
	  y=(long)(series[sn+DELAY]*epsinv)&imax;
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
	n2=(sn+MINDIST<length-(MAXEMBED-1)*DELAY)?sn+MINDIST:
	  length-(MAXEMBED-1)*DELAY-1;
	for (i1=n1;i1<=n2;i1++)
	  if ((oscr[i1] < n))
	    lnorm--;
      }

      make_c2_dim(n);
      make_c2_1(n);
      for (i=imin;i<HOWOFTEN;i++)
	norm[i] += (double)(lnorm);
    }
    
    if (((time(&mytime)-lasttime) > WHEN) || (n == (nmax-1)) || 
	(imin > howoften1)) {
      time(&lasttime);
      fstat=fopen(outstat,"w");
      fprintf(fstat,"Center points treated so far= %ld\n",n);
      fprintf(fstat,"Maximal epsilon in the moment= %e\n",epsm[imin]);
      fclose(fstat);
      fout=fopen(outc1,"w");
      fprintf(fout,"#center= %ld\n",n);
      for (i=0;i<MAXEMBED;i++) {
	fprintf(fout,"#dim= %ld\n",i+1);
	eps=EPSMAX1*epsfactor;
	for (j=0;j<HOWOFTEN;j++) {
	  eps /= epsfactor;
	  if (norm[j] > 0.0)
	    fprintf(fout,"%e %e %e\n",eps*interval,found[i][j]/norm[j],
		    eps);
	}
	fprintf(fout,"\n\n");
      }
      fclose(fout);
      fout=fopen(outh1,"w");
      fprintf(fout,"#center= %ld\n",n);
      for (i=0+1;i<MAXEMBED;i++) {
	fprintf(fout,"#dim= %ld\n",i+1);
	eps=EPSMAX1*epsfactor;
	for (j=0;j<HOWOFTEN;j++) {
	  eps /= epsfactor;
	  if ((found[i-1][j] > 0.0) && (found[i][j] > 0.0))
	    fprintf(fout,"%e %e %e\n",eps*interval,
		    log(found[i-1][j]/found[i][j]),eps);
	}
	fprintf(fout,"\n\n");
      }
      fclose(fout);
      fout=fopen(outd1,"w");
      fprintf(fout,"#center= %ld\n",n);
      for (i=0;i<MAXEMBED;i++) {
	fprintf(fout,"#dim= %ld\n",i+1);
	eps=EPSMAX1;
	for (j=1;j<HOWOFTEN;j++) {
	  eps /= epsfactor;
	  if ((found[i][j] > 0.0) && (found[i][j-1] > 0.0))
	    fprintf(fout,"%e %e %e\n",eps*interval,
		    log(found[i][j-1]/found[i][j]
			/norm[j-1]*norm[j])/lnfac,eps);
	}
	fprintf(fout,"\n\n");
      }
      fclose(fout);
      if (imin > howoften1)
	exit(0);
    }
  }

  return 0;
}
