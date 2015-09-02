/* file: mse.c			M. Costa		1 August 2004
				Last revised:		4 August 2004 (GBM)
-------------------------------------------------------------------------------
mse: calculates multiscale entropy (MSE) of one or multiple data sets
Copyright (C) 2004 Madalena Costa

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the author by e-mail (mcosta@fsa.harvard.edu).  For updates to
this software, please visit PhysioNet (http://www.physionet.org/).
_______________________________________________________________________________

Compile this program by
    gcc -o mse -O mse.c -lm

There are two major steps in the calculations performed by mse:
1. Time series are coarse-grained.
2. Sample entropy (SampEn) is calculated for each coarse-grained time series.

Output file:
1st line: shows the r value.
2nd line: shows the m values. 

After the 2nd line there are several columns: the first column (of integers)
is the scale factor. The following columns are SampEn values for coarse-grained
time series calculated for the values of r and m specified. If the option for
calculating MSE for several r values is chosen a new line containing the new r
value and new columns with the corresponding results are written.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXSTR 1250     /* maximum string size */
#define DATA_MAX 200000 /* maximum number of data points */  
#define M_MAX 10        /* maximum value of the parameter m */
#define SCALE_MAX 40    /* maximum scale value */ 
#define R_STEP_MAX 10   /* maximum number of different r values */  
#define FILEMAX 100     /* maximum number of data files in the list file*/

/* Global variables */
char *prog, line[MAXSTR];
double SE[FILEMAX][R_STEP_MAX][SCALE_MAX][M_MAX];
double *u, *y, r_min, r_max, r_step;
int c, nlin, m_min, m_max, m_step, scale_max, scale_step, i_min, i_max;

static char file[500];

FILE *fl;
FILE *pin;

main(argc, argv)
int argc;
char *argv[];
{

    int i, j, k, l, m, nfile, flag; 
    double r, sd, tmp;
 
    /* Function prototypes */
    void ReadData(void);
    double *array(int n);
    double StandardDeviation(void);
    void CoarseGraining (int j);
    void SampleEntropy (int ll, double r, double sd, int j);
    void PrintResults (int nfile);
    void PrintAverageResults (int nfile);

    /* Initialize variables. */
    prog = argv[0];
    scale_max = 20;
    scale_step = 1;
    m_min = 2;
    m_max = 2;
    m_step = 1;
    i_min = 0;
    i_max = 39999;
    r_min = 0.15;
    r_max = 0.15;
    r_step = 0.05;
    c = 0;
    nfile = 0;
    flag = 0;
  
    /* Read and interpret the command line. */
    i = 0;
    while (++i < argc && *argv[i] == '-') {
	switch(argv[i][1]) {

	  case 'F':
	    if ((fl = fopen(argv[++i], "r")) == NULL) {
		fprintf(stderr, "%s [-F]: can't open input file %s\n",
			prog, argv[i]);
		exit(1);
	    }
	    flag=1;
	    break;

	  case 'n':
	    if ((scale_max=atoi(argv[++i])) <= 0 || scale_max > SCALE_MAX) {
		fprintf(stderr,
	     "%s [-n]: maximum scale must be between 1 and %d (default: 20)\n",
			prog, SCALE_MAX);
		exit(1);
	    }
	    break;

	  case 'a':
	    if ((scale_step=atoi(argv[++i])) <= 0 || scale_step >= SCALE_MAX) {
		fprintf(stderr,
	    "%s [-a]: scale increment must be between 1 and %d (default: 1)\n",
			prog, SCALE_MAX);
		exit(1);
	    }
	    break;

	  case 'm':
	    if ((m_min = atoi(argv[++i])) >= M_MAX || m_min <= 0 ) {
		fprintf(stderr,
	    "%s [-m]: minimum m value must be between 1 and %d (default: 2)\n",
			prog, M_MAX);
		exit(1);
	    }
	    break;

	  case 'M':
	    if ((m_max = atoi(argv[++i])) >= M_MAX || m_max <= 0) {
		fprintf(stderr,
	    "%s [-M]: maximum m value must be between 1 and %d (default: 2)\n",
			prog, M_MAX);
		exit(1);
	    }
	    break;

	  case 'b':
	    if ((m_step = atoi(argv[++i])) <= 0 || m_step > M_MAX) {
		fprintf(stderr,
	        "%s [-b]: m increment must be between 1 and %d (default: 1)\n",
			prog, M_MAX);
		exit(1);
	    }
	    break;

	  case 'r':
	    if ((r_min = atof(argv[++i])) <= 0) {
		fprintf(stderr,
	         "%s [-r]: minimum r must be greater than 0 (default: 0.15)\n",
			prog);
		exit(1);
	    } 
	    break;

	  case 'R':
	    if ((r_max = atof(argv[++i])) <= 0) {
		fprintf(stderr,
		 "%s [-R]: maximum r must be greater than 0 (default: 0.15)\n",
			prog);
		exit(1);
	    } 
	    break;

	  case 'c':
	    if ((r_step=atof(argv[++i]))<=0||r_step<(r_max-r_min)/R_STEP_MAX) {
		fprintf(stderr,
	      "%s [-c]: r increment must be greater than %g (default: 0.05)\n",
			prog, (r_max-r_min)/R_STEP_MAX);
	    }
	    exit(1);
	    break;
    
	  case 'i':
	    if ((i_min = atoi(argv[++i])) < 0) {
		fprintf(stderr,
	           "%s [-i]: minimum i must not be less than 0 (default: 0)\n",
			prog);
		exit(1);
	    }
	    break;

	  case 'I':
	    if ((i_max = atoi(argv[++i])) <= 0 || i_max <= i_min) {
		fprintf(stderr,
			"%s [-I]: maximum i must be greater than %d "
			"(default: number of data points)\n",
			prog, i_min);
		exit(1);
	    }
	    break;

	  default:
	    usage();
	    exit(1);
	}
    }

    if (m_max < m_min) {
	tmp = m_max;
	m_max = m_min;
	m_min = tmp;
    }
      
    if (r_max < r_min){
	tmp = r_max;
	r_max = r_min;
	r_min = tmp;
    }
   
    /* Memory allocation. */
    u = array(DATA_MAX);
    y = array(DATA_MAX);


    /* Process a single data file. */
    if (flag == 0) {  
	l = 0;
	nfile = 1;
    
	/* Read data from stdin. */
	pin = stdin;
	ReadData();
    
	/* Calculate standard deviation. */
	sd = StandardDeviation();
    
	/* Perform the coarse-graining procedure. */
	for (j = 1; j <= scale_max; j += scale_step){      
	    CoarseGraining(j);

	    /* Calculate SampEn for each scale and each r value. */
	    c = 0;
	    for (r = r_min; r <= (r_max*1.0000000001); r += r_step){
		SampleEntropy(l, r, sd, j);
		c++;
	    }   
	}
    
	/* Print results. */
	PrintResults(nfile);
    }
  

    /* Process multiple data files. */
    if (flag == 1) {

	/* Read the list of data files. */  
	for (l = 0; fscanf(fl, "%s", file) == 1; l++) {  
	    nfile++;   /*count the number of data files*/
      
	    if ((pin = fopen(file, "r")) == NULL) {   /* open each data file */
		fprintf(stderr, "%s : Cannot open file %s\n", prog, file);
		exit(1);
	    }
	    
	    /* Read the data. */
	    ReadData();
           
	    /* Calculate the standard deviation. */
	    sd = StandardDeviation ();
      
	    /* Perform the coarse-graining procedure. */
	    for (j = 1; j <= scale_max; j += scale_step) {
		CoarseGraining(j);
		c = 0;

		/* Calculate SampEn for each scale and each r value. */
		for (r = r_min; r <= (r_max*1.0000000001) ; r += r_step) {
		    SampleEntropy (l, r, sd, j);
		    c++;
		}
	    }      
	}

	/* Print results. */
	PrintResults(nfile);

	/* Print average results. */
	if (nfile > 1)
	    PrintAverageResults(nfile);
      
	fclose(pin);
	fclose(fl);
    }
  
    free(u);
    free(y);
    exit(0);
}

double *array(int n)
{
    int i;
    double *a;
  
    if ((a = calloc (n, sizeof (double))) == NULL) {
	fprintf(stderr, "%s : insufficient memory\n", prog);
	exit(2);
    }
    return (a);
}

void ReadData(void)
{
    int j = -1;

    while (fgets(line, MAXSTR, pin) != NULL) {
	j++;
	if (j >= i_min && j <= i_max) {
	    sscanf(line, "%lf", &u[j-i_min]);
	    nlin=j-i_min+1;
	}
    }
}

double StandardDeviation(void)
{
    double sum=0.0, sum2=0.0, sd;
    int j;
  
    for (j = 0; j < nlin; j++) {
	sum += u[j];
	sum2 += u[j] * u[j];
    }
    sd = sqrt((sum2 - sum*sum/nlin)/(nlin - 1));
    return (sd);
}


void CoarseGraining(int j)
{
    int i, k;

    for (i = 0; i < nlin/j; i++) {
	y[i] = 0;
	for (k = 0; k < j; k++)
	    y[i] += u[i*j+k];
	y[i] /= j; 
    }
}


void SampleEntropy(int ll, double r, double sd, int j)
{
    int i, k, l, nlin_j; 
    int cont[M_MAX+1];
    double r_new;

    nlin_j = (nlin/j) - m_max; 
    r_new = r*sd;              

    for (i = 0; i < M_MAX; i++)
	cont[i]=0;

    for (i = 0; i < nlin_j; ++i) {
	for (l = i+1; l < nlin_j; ++l) { /*self-matches are not counted*/
	    k = 0;
	    while (k < m_max && fabs(y[i+k] - y[l+k]) <= r_new)
		cont[++k]++;
	    if (k == m_max && fabs(y[i+m_max] - y[l+m_max]) <= r_new)
		cont[m_max+1]++;
	} 
    }     

    for (i = 1; i <= m_max; i++)
	if (cont[i+1] == 0 || cont[i] == 0)
	    SE[ll][c][j][i] = -log((double)1/((nlin_j)*(nlin_j-1)));
	else
	    SE[ll][c][j][i] = -log((double)cont[i+1]/cont[i]);
}


void PrintResults(int nfile)
{
    int j, m, k, l, i;

    printf ("\n");
    for (m = m_min; m <= m_max; m += m_step)
	for (k = 0; k < c; k++) {
	    printf ("\nm = %d,  r = %.3f\n\n", m, r_min+k*r_step);
	    if (nfile > 1) {
		fseek(fl, 0, SEEK_SET);
		for(l = 0; fscanf(fl, "%s", file) == 1; l++)  
		    printf("\t%.6s", file);
		printf("\n");
	    }
	    for (j = 1; j <= scale_max; j += scale_step) {
		printf("%d\t", j);
		for (l=0; l<nfile; l++)
		    printf("%.3lf\t", SE[l][k][j][m]);
		printf("\n");   
	    }  
	}  
}

void PrintAverageResults(int nfile)
{
    int k, m, j, i, l; 
    double av, av2, sd1;
 
    printf("\n**************************\n");
    printf("Mean and SD over all files\n");
    printf("**************************\n");
      
    for (k = 0; k < c; k++) {
	printf("\n");
	
	for (m = m_min; m <= m_max; m += m_step)
	    printf("\tm=%d, r=%5.3lf", m, r_min+k*r_step);
	printf("\n");
	for (m = m_min; m <= m_max; m += m_step)
	    printf("\tmean\tsd");
	printf("\n");
	
	for (j = 1; j <= scale_max; j += scale_step) {
	    printf("%d\t", j);
	    for (i = m_min; i <= m_max; i += m_step) {
		av = 0.0;
		av2 = 0.0;
		/* Calculate entropy mean values and SD over all files. */ 
		for (l = 0; l < nfile; l++) {
		    av += SE[l][k][j][i];
		    av2 += (SE[l][k][j][i]*SE[l][k][j][i]);
		}
		sd1 = sqrt((av2-av*av/nfile) / (nfile-1));
		av /= nfile;
		printf("%.3lf\t%.3lf\t", av, sd1);
	    }
	    printf("\n"); 
	}
    }
}

usage()
{
    fprintf(stderr, "usage: %s [options]\n", prog);
    fprintf(stderr, "\nTo calculate MSE for a single data file:\n"
	    "    %s <datafile >outputfile\n"
	    "To calculate MSE for multiple data files:\n"
	    "    %s -F listfile >outputfile\n"
	    "(where listfile contains a list of data files).\n\n", prog, prog);
    fprintf(stderr, "Data files should contain a single column of numbers\n");
    fprintf(stderr, "Options may include:\n");
    fprintf(stderr, "  -a N   set scale increment to N [1-%d; default: 1]\n",
	    SCALE_MAX);
    fprintf(stderr, "  -b N   set m increment to N [1-%d; default: 1]\n",
	    M_MAX);
    fprintf(stderr, "  -c X   set r increment to X [>%g; default: 0.05]\n",
	    (r_max-r_min)/R_STEP_MAX);
    fprintf(stderr, "  -i N   set minimum i to N [0-39998; default: 0]\n");
    fprintf(stderr, "  -I N   set maximum i to N [1-39999: default: 39999]\n");
    fprintf(stderr, "  -m N   set minimum m to N [1-%d; default: 2]\n", M_MAX);
    fprintf(stderr, "  -M N   set maximum m to N [1-%d; default: 2]\n", M_MAX);
    fprintf(stderr, "  -n N   set maximum scale to N [1-%d; default: 20]\n",
	    SCALE_MAX);
    fprintf(stderr, "  -r X   set minimum r to X [>0; default: 0.15]\n");
    fprintf(stderr, "  -R X   set maximum r to X [>0; default: 0.15]\n");
    fprintf(stderr,
 "Option arguments indicated as N are integers; those shown as X may be given"
 "\nin any floating point format.  For further information, visit:\n"
 "    http://physionet.org/physiotools/mse/\n");
}



