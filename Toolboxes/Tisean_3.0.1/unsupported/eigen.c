/* Copyright (C) Marcus Richter and Thomas Schreiber (1997)         */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void esort(double *d, double **v, long n, long nq)
{
   long k,j,i;
   double p;

   for (i=0; i<nq; i++) {
      p=d[k=i];
      for (j=i+1;j<n;j++)
         if (d[j] >= p) p=d[k=j];
      if (k != i) {
         d[k]=d[i];
         d[i]=p;
         for (j=0; j<n; j++) {
            p=v[j][i];
            v[j][i]=v[j][k];
            v[j][k]=p;
         }
      }
   }
}

double pythag(double a, double b)
{
   double absa,absb;
   absa=fabs(a);
   absb=fabs(b);
   if (absa > absb) 
      return absa*sqrt(1.0+(absb/absa)*(absb/absa));
   else 
      return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

void eig1(double **a, long n, double *d, double *e)
{
   long l,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>0;i--) {
      l=i-1;
      h=scale=0.0;
      if (l > 0) {
         for (k=0; k<=l; k++)
            scale += fabs(a[i][k]);
         if (scale == 0.0)
            e[i]=a[i][l];
         else {
            for (k=0; k<=l; k++) {
               a[i][k] /= scale;
               h += a[i][k]*a[i][k];
            }
            f=a[i][l];
            g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i][l]=f-g;
            f=0.0;
            for (j=0; j<=l; j++) {
               a[j][i]=a[i][j]/h;
               g=0.0;
               for (k=0; k<=j; k++)
                  g += a[j][k]*a[i][k];
               for (k=j+1; k<=l; k++)
                  g += a[k][j]*a[i][k];
               e[j]=g/h;
               f += e[j]*a[i][j];
            }
            hh=f/(h+h);
            for (j=0; j<=l; j++) {
               f=a[i][j];
               e[j]=g=e[j]-hh*f;
               for (k=0; k<=j; k++)
                  a[j][k] -= (f*e[k]+g*a[i][k]);
            }
         }
      } else
         e[i]=a[i][l];
      d[i]=h;
   }
   d[0]=0.0;
   e[0]=0.0;
   for (i=0; i<n; i++) {
      l=i-1;
      if (d[i]) {
         for (j=0; j<=l; j++) {
            g=0.0;
            for (k=0; k<=l; k++)
               g += a[i][k]*a[k][j];
            for (k=0; k<=l; k++)
               a[k][j] -= g*a[k][i];
         }
      }
      d[i]=a[i][i];
      a[i][i]=1.0;
      for (j=0; j<=l; j++) a[j][i]=a[i][j]=0.0;
   }
}

void eig2(double d[], double e[], long n, double **z)
{
   double pythag(double a, double b);
   long m,l,iter,i,k;
   double s,r,p,g,f,dd,c,b;

   for (i=1; i<n; i++) e[i-1]=e[i];
   e[n-1]=0.0;
   for (l=0; l<n; l++) {
      iter=0;
      do {
         for (m=l;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if ((double)(fabs(e[m])+dd) == dd) break;
         }
         if (m != l) {
	    if (iter++ == 30) {
               fprintf(stderr,"Iteration count exceeded in eig2");
	       exit(5);
	    }
            g=(d[l+1]-d[l])/(2.0*e[l]);
            r=pythag(g,1.0);
            g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s=c=1.0;
            p=0.0;
            for (i=m-1;i>=l;i--) {
               f=s*e[i];
               b=c*e[i];
               e[i+1]=(r=pythag(f,g));
               if (r == 0.0) {
                  d[i+1] -= p;
                  e[m]=0.0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2.0*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0; k<n; k++) {
                  f=z[k][i+1];
                  z[k][i+1]=s*z[k][i]+c*f;
                  z[k][i]=c*z[k][i]-s*f;
               }
            }
            if (r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]=g;
            e[m]=0.0;
         }
      } while (m != l);
   }
}

void eigen(double **c, long m, long nq)
{
   static double *d=NULL, *e=NULL;

   d=(double *)realloc(d,m*sizeof(double));
   e=(double *)realloc(e,m*sizeof(double));
    
   eig1(c,m,d,e); 
   eig2(d,e,m,c);
   esort(d,c,m,nq);
}

