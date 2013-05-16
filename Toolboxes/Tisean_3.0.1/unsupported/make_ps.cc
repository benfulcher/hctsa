/*Author: Rainer Hegger, last modified: Mar 20, 1999
  I'm still working on it. Especially the tics and the title are not running
  properly yet
  */
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "tsa.h"

unsigned long length=ULONG_MAX;
unsigned int delay=1;
unsigned int column=1;
unsigned int maxcolor=32;
unsigned int linewidth=1;
unsigned int textsize=15;
unsigned long exclude=0;
double size=16.0;
char pline=1;
char tics=0;
char onecolumn=1;
char *outfile=NULL;
char *infile=NULL;
char *done;
char *cformat;

void error(char *str1,char *str2="")
{
  cerr << str1 << " " << str2 << '\n';
  exit(127);
}

void show_options(char *progname)
{
  cerr << " Usage: " << progname << " [Options]\n";
  cerr << " Options:\n";
  cerr << "Everything not being a valid option will be interpreted as a"
    " possible datafile.\nIf no datafile is given stdin is read. Just - also"
    " means stdin\n";
  cerr << "  -l # of data [default whole file]\n";
  cerr << "  -x # of rows to ignore [default: 0]\n";
  cerr << "  -c column to read (only if -3 is not set) [default: 1]\n";
  cerr << "  -d delay [default 1]\n";
  cerr << "  -C # of colors [default 32]\n";
  cerr << "  -s size in cm [default is 16 cm]\n";
  cerr << "  -w line[point]width [default 1]\n";
  cerr << "  -t textwidth [default 15]\n";
  cerr << "  -p plot with points [default with lines]\n";
  cerr << "  -a make tics [default is notics]\n";
  cerr << "  -3 read 3 columns [default read 1 col. and build delay vect.]\n";
  cerr << "  -F which columns (separated by commas) [default: 1,2,3]\n";
  cerr << "  -o outfile [default datafile.ps]\n";
  cerr << "  -h show this options\n\n";
  exit(127);
}

void scan_options(int n,char **str)
{
  char *out;
  
  if ((out=check_option(str,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(str,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(str,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(str,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(str,n,'C','u')) != NULL)
    sscanf(out,"%u",&maxcolor);
  if ((out=check_option(str,n,'s','f')) != NULL)
    sscanf(out,"%lf",&size);
  if ((out=check_option(str,n,'w','u')) != NULL)
    sscanf(out,"%u",&linewidth);
  if ((out=check_option(str,n,'t','u')) != NULL)
    sscanf(out,"%u",&textsize);
  if ((out=check_option(str,n,'3','n')) != NULL)
    onecolumn=0;
  if ((out=check_option(str,n,'p','n')) != NULL)
    pline=0;
  if ((out=check_option(str,n,'a','n')) != NULL)
    tics=1;
  if ((out=check_option(str,n,'F','s')) != NULL)
    cformat=out;
  if ((out=check_option(str,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      outfile=out;
}

int main(int argc,char **argv)
{
  char stdi=0;
  unsigned long i,j;
  double min,max;
  double rot,gruen,blau,rgbstep;
  double dxmax,dxmin,dymin,dymax,dzmin,dzmax;
  double **s3,*s1;
  double *sx,*sy,*sz;
  double zvar,zav;
  
  if (scan_help(argc,argv))
    show_options(argv[0]);

  cformat=(char*)calloc(1,1);
  cformat[0]='\0';

  scan_options(argc,argv);
  infile=search_datafile(argc,argv,&column);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL)
    if (!stdi) {
      outfile=(char*)calloc(strlen(infile)+4,1);
      outfile=strcpy(outfile,infile);
      outfile=strcat(outfile,".ps");
    }
    else {
      outfile=(char*)calloc(9,1);
      strcpy(outfile,"stdin.ps");
    }
  test_outfile(outfile);

  if (onecolumn) {
    s1=(double*)get_series(infile,&length,exclude,column);
    sx=(double*)malloc(sizeof(double)*length);
    sy=(double*)malloc(sizeof(double)*length);
    sz=(double*)malloc(sizeof(double)*length);
    for (i=2*delay;i<length;i++) {
      sx[i-2*delay]=s1[i-2*delay];
      sy[i-2*delay]=s1[i-delay];
      sz[i-2*delay]=s1[i];
    }
    length -= 2*delay;
  }
  else {
    s3=(double**)get_multi_series(infile,&length,exclude,3,cformat);
    sx=s3[0];
    sy=s3[1];
    sz=s3[2];
  }
  rescale_data(sx,length,&dxmin,&dxmax);
  rescale_data(sy,length,&dymin,&dymax);
  rescale_data(sz,length,&dzmin,&dzmax);

  ofstream fout(outfile);

  if (fout == NULL)
    error("Could not create:",outfile);

  done=(char*)calloc(length,1);
  for (i=0;i<length;i++)
    done[i]=0;

  fout << "%!PS-Adobe-3.0\n";
  fout << "%%Creator: make_ps Version 1.0\n";
  //fout << "%%BoundingBox: 0 0 510 510\n";
  fout << "%%EndComments\n";
  fout << "/size " << size << " def\n";
  fout << "/linewidth " << linewidth << " def\n";
  fout << "/outlinewidth 4 def\n";
  fout << "/ticlength 0.2 def\n";
  fout << "/textsize " << textsize << " def\n";
  fout << "/factor 10 def\n";
  fout << "/sfactor {1 factor div} def\n";
  fout << "/cm {28.35 mul} def\n";
  fout << "/ccm {cm factor mul} def\n";
  fout << "/xpage 21 cm def\n";
  fout << "/ypage 29.7 cm def\n";
  fout << "/xoffset {xpage size cm sub 2 div} def\n";
  fout << "/yoffset {ypage size cm sub 2 div} def\n";
  fout << "/soffset {linewidth outlinewidth add 2 div} def\n";
  fout << "/ssoffset {soffset factor mul} def\n";
  fout << "/box [soffset neg soffset neg\n"
    " size cm soffset 2 mul add size cm soffset 2 mul add] def\n";
  fout << "/xmaketic {size mul ccm ssoffset neg moveto\n"
    "0 ticlength ccm rlineto stroke} def\n";
  fout << "/ymaketic {ssoffset neg exch size ccm mul moveto\n"
    "ticlength ccm 0 rlineto stroke} def\n";
  fout << "/xtextpos {size mul ccm 1.05 textsize mul\n"
    " factor mul ssoffset add neg moveto} def\n";
  fout << "/ytextpos {ssoffset neg exch size mul ccm moveto\n"
    "0 textsize factor mul 3 div neg rmoveto} def\n";
  fout << "/pl {size ccm mul 4 1 roll size ccm mul 4 1 roll\n"
    " size ccm mul 4 1 roll size ccm mul 4 3 roll\n" 
    "1 setlinecap moveto lineto stroke} def\n";
  fout << "/pp {size ccm mul 2 1 roll size ccm mul 2 1 roll\n"
    "linewidth factor mul 2 div pop moveto linewidth 0 rlineto 0 linewidth"
    " rlineto -1 linewidth mul 0 rlineto closepath fill stroke} def\n";
  //moveto 0 0 rlineto stroke} def\n";
  fout << "/Optima-Bold findfont\ntextsize factor mul"
    " scalefont\nsetfont\n";
  fout << "%%EndProlog\n";
  fout << "xoffset yoffset translate\n";
  fout << "currentgray\n0.85 setgray\nbox rectfill\n";
  fout << "outlinewidth setlinewidth\n";
  fout << "setgray\n";
  fout << "box rectstroke\n";
  fout << "sfactor sfactor scale\n";
  fout << "linewidth factor mul setlinewidth\n";
  fout << "newpath\n";
  rgbstep=maxcolor/3.0;
  rot=0.0;
  gruen=0.0;
  blau=1.0;

  for (i=0;i<=maxcolor;i++) {
    min=(double)i/maxcolor;
    max=(double)(i+1)/maxcolor;
    
    //rot=(double)i/maxcolor;
    //blau=1.0-(double)i/maxcolor;
    if (i < rgbstep) {
      blau=1.-(double)i/rgbstep;
      rot=gruen=0.0;
    }
    else {
      if (i < 2*rgbstep) {
	//gruen=1.0-(i-rgbstep)/(double)rgbstep;
	blau=0.0;
	gruen=0.0;
	rot=(i-rgbstep)/(double)rgbstep;
      }
      else {
	blau=0.;
	rot=1.;
	gruen=(i-2*rgbstep)/(double)rgbstep;
      }
    }

    fout << rot << " " << gruen << " " << blau << " setrgbcolor\n";
    for (j=0;j<length-1;j++)
      if (!done[j])
	if ((sz[j] >= min) && (sz[j] < max)) {
	  if (pline)
	    fout << sx[j] << " " << sy[j] << " " << sx[j+1] 
		 << " " << sy[j+1] << " pl\n";
	  else
	    fout << sx[j] << " " << sy[j] << " pp\n";
	  done[j]=1;
	}
  }
  if (tics) {
    fout << "outlinewidth 2 div factor mul setlinewidth\n";
    fout << "0.0 setgray\n";

    int tx1=(int)(0.5*dxmax+dxmin+0.5);
    double x1=(tx1-dxmin)/dxmax;
    int tx2=(int)(0.25*dxmax+dxmin+0.5);
    double x2=(tx2-dxmin)/dxmax;
    int tx3=2*tx1-tx2;
    double x3=(tx3-dxmin)/dxmax;

    int ty1=(int)(0.5*dymax+dymin+0.5);
    double y1=(tx1-dymin)/dymax;
    int ty2=(int)(0.25*dymax+dymin+0.5);
    double y2=(ty2-dymin)/dymax;
    int ty3=2*ty1-ty2;
    double y3=(ty3-dymin)/dymax;

    if ((x1>= 0.0) && (x1 <= 1.0)) {
      fout << x1 << " dup xmaketic xtextpos\n";
      fout << "(";
      fout << tx1 << ")\n";
      fout << "dup stringwidth pop 2 div neg 0 rmoveto show\n";
    }
    else 
      cerr << "X-Tic 1 out of range. Won't plot it\n";
    if ((y1 >= 0.0) && (y1 <= 1.0)) {
      fout << y1 << " dup ymaketic ytextpos\n";
      fout << "(";
      fout << ty1 << ")\n";
      fout << "dup stringwidth pop neg textsize 0.4 mul factor mul"
	" sub 0 rmoveto show\n";
    }
    else
      cerr << "Y-Tic 1 out of range. Won't plot it\n";

    if ((x2 >= 0.0) && (x2 <= 1.0)) {
      fout << x2 << " dup xmaketic xtextpos\n";
      fout << "(";
      fout << tx2 << ")\n";
      fout << "dup stringwidth pop 2 div neg 0 rmoveto show\n";
    }
    else 
      cerr << "X-Tic 2 out of range. Won't plot it\n";
    if ((y2 >= 0.0) && (y2 <= 1.0)) {
      fout << y2 << " dup ymaketic ytextpos\n";
      fout << "(";
      fout << ty2 << ")\n";
      fout << "dup stringwidth pop neg textsize 0.4 mul factor mul"
	" sub 0 rmoveto show\n";
    }
    else
      cerr << "Y-Tic 2 out of range. Won't plot it\n";

    if ((x3 >= 0.0) && (x3 <= 1.0)) {
      fout << x3 << " dup xmaketic xtextpos\n";
      fout << "(";
      fout << tx3 << ")\n";
      fout << "dup stringwidth pop 2 div neg 0 rmoveto show\n";
    }
    else 
      cerr << "X-Tic 3 out of range. Won't plot it\n";
    if ((y3 >= 0.0) && (y3 <= 1.0)) {
      fout << y3 << " dup ymaketic ytextpos\n";
      fout << "(";
      fout << ty3 << ")\n";
      fout << "dup stringwidth pop neg textsize 0.4 mul factor mul"
	" sub 0 rmoveto show\n";
    }
    else
      cerr << "Y-Tic 3 out of range. Won't plot it\n";
  }
  fout << "showpage\n";
  fout << "%%EOF\n";
  fout.close();

  return 0;
}
