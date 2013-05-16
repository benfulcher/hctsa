/* certain changes (mainly global names) by Andreas Schmitz (1997) */
/* certain changes (mainly cleanup) by Thomas Schreiber (1998)     */

/**************************************************************************

	Javier Soley, Ph. D,   FJSOLEY @UCRVM2.BITNET
	Escuela de Física y Centro de Investigaciones Geofísicas
	Universidad de Costa Rica

***************************************************************************/

/*	Computes the DISCRETE FOURIER TRANSFORM of very long data series.
 *
 *   This functions are translations from the fortran program in
 *
 *   R. C. Singleton, An algorithm for computing the mixed radix fast
 *   Fourier transform 
 *
 *	IEEE Trans. Audio Electroacoust., vol. AU-17, pp. 93-10, June 1969.
 *   Some features are:
 *
 *		1-) Accepts an order of transform that can be factored not only
 *		    in prime factors such 2 and 4, but also including odd factors
 *		    as 3, 5, 7, 11, etc.
 *		2-) Generates sines and cosines recursively and includes
 *		    corrections for truncation errors.
 *		3-) The original subroutine accepts multivariate data. This 
 *			translation does not implement that option (because I
 *			do not needed right now).
 *
 *	Singleton wrote his subroutine in Fortran and in such a way that it 
 *	could be ported allmost directly to assembly language. I transcribed
 *   it to C with little effort to make it structured. So I apologize to
 *   all those C purists out there!!!!!!!!
 *
 */

				/*  Version 2.0 March/30/92 */
/* Includes */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

/* Defines */
#ifndef M_PI
   #define M_PI   3.14159265358979323846       /* should have been in math.h */
#endif
#ifndef M_PI_2
   #define M_PI_2 1.57079632679489661923       /* should have been in math.h */
#endif
#define	  TWO_PI	((double)2.0 * M_PI)
#define   MAXF  20000
#define   MAXP  20000 

/* Globals */

long nn_gl,m_gl,flag_gl,jf_gl,jc_gl,kspan_gl,ks_gl, kt_gl,nt_gl,kk_gl,i_gl;
double c72_gl, s72_gl, s120_gl, cd_gl, sd_gl, rad_gl, radf_gl;
double at_gl[MAXF], bt_gl[MAXF];
long nfac_gl[MAXF];   
int inc_gl;
long np_gl[MAXP];


/* The functions */

void radix_2(double *a, double *b)
{ 
  long k1, k2;
  double ak, bk, c1, s1;
  kspan_gl >>= 1;
  k1 = kspan_gl +2;
  
  do {
    do  {
      k2 = kk_gl + kspan_gl;
      ak = a[k2-1];
      bk = b[k2-1];
      a[k2-1] = a[kk_gl-1] -ak;
      b[k2-1] = b[kk_gl-1] -bk;
      a[kk_gl-1] += ak;
      b[kk_gl-1] += bk;
      kk_gl = k2 + kspan_gl;
    } while ( kk_gl <= nn_gl);
    kk_gl = kk_gl - nn_gl;
  } while ( kk_gl <= jc_gl);
  
  if ( kk_gl > kspan_gl) 
    flag_gl = 1; 
  else
    {
      do {
	c1 = 1.0 - cd_gl;
	s1 = sd_gl;
	do {
	  do  {
	    do {
	      k2 = kk_gl + kspan_gl;
	      ak = a[kk_gl-1]- a[k2-1];
	      bk = b[kk_gl-1]- b[k2-1];
	      a[kk_gl-1] += a[k2-1];
	      b[kk_gl-1] += b[k2-1];
	      a[k2-1]  = c1*ak - s1*bk;
	      b[k2-1]  = s1*ak + c1*bk;
	      kk_gl = k2 + kspan_gl;
	    } while ( kk_gl < nt_gl );
	    k2 = kk_gl - nt_gl;
	    c1 = -c1;
	    kk_gl = k1 - k2;
	  } while ( kk_gl > k2 );
	  ak = c1- (cd_gl*c1+sd_gl*s1);
	  s1 = (sd_gl*c1-cd_gl*s1) +s1;
	  
	  /***** Compensate for truncation errors   *****/
	  
	  c1 = 0.5/(ak*ak+s1*s1)+0.5;
	  s1 *= c1;
	  c1 *= ak;
	  kk_gl += jc_gl;
	} while ( kk_gl < k2);
	k1 = k1 + inc_gl + inc_gl;
	kk_gl = (k1- kspan_gl) /2 + jc_gl;
      } while ( kk_gl <= jc_gl + jc_gl );
    }
}

void radix_4(int isn, double *a, double *b)
{
  long	k1, k2, k3;
  double  akp, akm, ajm, ajp, bkm, bkp, bjm, bjp;
  double  c1, s1, c2, s2, c3, s3;

  kspan_gl /= 4;
 cuatro_1: 
  c1 = 1.0;
  s1 = 0;
  do {
    do {
      do {
	k1  =  kk_gl + kspan_gl;
	k2  =  k1 + kspan_gl;
	k3  =  k2 + kspan_gl;
	akp =  a[kk_gl-1] + a[k2-1];
	akm =  a[kk_gl-1] - a[k2-1];
	ajp =  a[ k1-1] + a[k3-1];
	ajm =  a[ k1-1] - a[k3-1];
	a[kk_gl-1] = akp + ajp;
	ajp = akp - ajp;
	bkp = b[kk_gl-1] + b[k2-1];
	bkm = b[kk_gl-1] - b[k2-1];
	bjp = b[k1-1] + b[k3-1];
	bjm = b[k1-1] - b[k3-1];
	b[kk_gl-1] = bkp + bjp;
	bjp = bkp - bjp;
	if ( isn < 0) goto cuatro_5;
	akp = akm - bjm;
	akm = akm + bjm;
	bkp = bkm + ajm;
	bkm = bkm - ajm;
	if (s1 == 0.0) goto cuatro_6;
      cuatro_3:		a[ k1-1] = akp*c1 - bkp*s1;
	b[ k1-1] = akp*s1 + bkp*c1;
	a[ k2-1] = ajp*c2 - bjp*s2;
	b[ k2-1] = ajp*s2 + bjp*c2;
	a[ k3-1] = akm*c3 - bkm*s3;
	b[ k3-1] = akm*s3 + bkm*c3;
	kk_gl = k3 + kspan_gl;
      }  while ( kk_gl <= nt_gl);
      
    cuatro_4: 
      c2 = c1 - (cd_gl*c1 + sd_gl*s1);
      s1 = (sd_gl*c1 - cd_gl*s1) + s1;
      
      /***** Compensate for truncation errors *****/
      
      c1 = 0.5 / (c2*c2 + s1*s1) +0.5;
      s1 = c1 * s1;
      c1 = c1 * c2;
      c2 = c1*c1 - s1*s1;
      s2 = 2.0 * c1 *s1;
      c3 = c2*c1 - s2*s1;
      s3 = c2*s1 + s2*c1;
      kk_gl = kk_gl -nt_gl + jc_gl;
    } while ( kk_gl <= kspan_gl);
    
    kk_gl = kk_gl - kspan_gl + inc_gl;
    if ( kk_gl <= jc_gl)  goto cuatro_1;
    if ( kspan_gl == jc_gl)  flag_gl =1;
    s1 = s1 + 0.0;
    return;
  cuatro_5:
    akp = akm + bjm;
    akm = akm - bjm;
    bkp = bkm - ajm;
    bkm = bkm + ajm;
    if (s1 != 0.0) goto cuatro_3;
  cuatro_6:
    a[k1-1] = akp;
    b[k1-1] = bkp;
    b[k2-1] = bjp;
    a[k2-1] = ajp;
    a[k3-1] = akm;
    b[k3-1] = bkm;
    kk_gl = k3 + kspan_gl;
  } while ( kk_gl <= nt_gl);
  goto cuatro_4;
}

/* Find prime factors of n */

void fac_des(long  n)
{
  long k, j, jj;
  k  = n;
  m_gl = 0;
  while ( k-(k / 16)*16 == 0 ) {
    m_gl++;
    nfac_gl[m_gl-1] = 4;
    k /= 16;
  } 
  j  = 3;
  jj = 9;
  do {
    while (k % jj == 0) {
      m_gl++;
      nfac_gl[m_gl-1] = j;
      k /= jj;
    }
    j += 2;
    jj = j * j;
  } while ( jj <= k);
  if (k <= 4) {
    kt_gl = m_gl;
    nfac_gl[m_gl] = k;
    if (k != 1) m_gl++;
  }
  else {
    if (k-(k / 4)*4 == 0) {
      m_gl++;
      nfac_gl[m_gl-1] = 2;
      k /= 4;
    }
    kt_gl = m_gl;
    j = 2;
    do {
      if (k % j == 0 ) {
	m_gl++;
	nfac_gl[m_gl-1] = j;
	k /= j;
      }
      j = ((j+1)/ 2)*2 + 1;
    } while ( j <= k);
  }
  if (kt_gl != 0) {
    j = kt_gl;
    do {
      m_gl++;
      nfac_gl[m_gl-1] = nfac_gl[j-1];
      j--;
    } while ( j != 0);
  }
}    

/* Permute the results to normal order  */

void permute(long n, double *a, double *b)
{
  long  k, j, k1, k2, k3, kspnn, maxf;
  
  double ak, bk;
  long  ii, jj;
  maxf = MAXF;
  np_gl[0] = ks_gl;
  if (kt_gl != 0) {
    k = kt_gl +kt_gl +1;
    if (m_gl < k) k--;
    j = 1;
    np_gl[k] = jc_gl;
    do {
      np_gl[j]   = np_gl[j-1] / nfac_gl[j-1];
      np_gl[k-1] = np_gl[k]   * nfac_gl[j-1];
      j++;
      k--;
    } while (j < k);
    k3 = np_gl[k];   
    kspan_gl = np_gl[1];
    kk_gl = jc_gl+1;
    k2 = kspan_gl + 1;  
    j = 1;
    
    do {
    	do {
    	  ak      = a[kk_gl-1];
    	  a[kk_gl-1] = a[k2-1];
    	  a[k2-1] = ak;
    	  bk      = b[kk_gl-1];
    	  b[kk_gl-1] = b[k2-1];
    	  b[k2-1] = bk;
    	  kk_gl     += inc_gl;
    	  k2     += kspan_gl;
    	} while ( k2 < ks_gl);
    ocho_30:  		do {
    	k2 -= np_gl[j-1];
    	j++;
    	k2 += np_gl[j];
    } while (k2 > np_gl[j-1]);
    	j = 1;
    ocho_40:  		j = j + 0;
    } while (kk_gl < k2);
    kk_gl += inc_gl;
    k2 += kspan_gl;
    if (k2 < ks_gl)  goto ocho_40;
    if (kk_gl < ks_gl)  goto ocho_30;
    jc_gl = k3;
  }
  
  if ( (2*kt_gl +1) < m_gl) {
    kspnn = np_gl[kt_gl];
    /* Permutation of square-free factors of n */
    j = m_gl - kt_gl;
    nfac_gl[j] = 1;
    do {
      nfac_gl[j-1] *= nfac_gl[j];
      j--;
    } while (j != kt_gl);
    kt_gl++;
    nn_gl = nfac_gl[kt_gl-1] -1;
    if (nn_gl > MAXP) {
      printf("product of square free factors exceeds allowed limit\n");
      exit(2);
    }
    jj =0;
    j=0;
    goto nueve_06;
  nueve_02:
    jj -= k2;
    k2 = kk_gl;
    k++;
    kk_gl = nfac_gl[k-1];
    do {
      jj += kk_gl;
      if ( jj >= k2 ) goto nueve_02;
      np_gl[j-1] = jj;
    nueve_06:
      k2 = nfac_gl[kt_gl-1];
      k = kt_gl+1;
      kk_gl = nfac_gl[k-1];
      j++;
    } while (j <= nn_gl);
    /* determine the permutation cycles  of length greater then 1 */
    j =0;
    goto nueve_14;
    do {
      do {		
	k = kk_gl;
	kk_gl = np_gl[k-1];
	np_gl[k-1] = -kk_gl;
      } while ( kk_gl != j);
      k3 = kk_gl;
    nueve_14:
      do {
	j++;
	kk_gl = np_gl[j-1];
      } while (kk_gl <0);
    } while ( kk_gl != j);
    np_gl[j-1] = -j;
    if (j != nn_gl) goto nueve_14;
    maxf *= inc_gl;				
    /* Reorder a and b following the permutation cycles */
    goto nueve_50;
    do {
      do {
	do { j--;} while (np_gl[j-1] <0);
	jj = jc_gl;
	do {
	  kspan_gl = jj;
	  if ( jj > maxf) kspan_gl = maxf;
	  jj -= kspan_gl;
	  k = np_gl[j-1];	
	  kk_gl = jc_gl*k + ii + jj;
	  k1 = kk_gl + kspan_gl;
	  k2 =0;
	  do {
	    k2++;
	    at_gl[k2-1] = a[k1-1];
	    bt_gl[k2-1] = b[k1-1];
	    k1 -= inc_gl;
	  } while (k1 != kk_gl);
	  do {
	    k1 = kk_gl + kspan_gl;
	    k2 = k1 - jc_gl*(k + np_gl[k-1]);
	    k = -np_gl[k-1];
	    do {
	      a[k1-1] = a[k2-1];
	      b[k1-1] = b[k2-1];
	      k1 -= inc_gl;
	      k2 -= inc_gl;
	    } while (k1 != kk_gl);
	    kk_gl = k2;
	  } while ( k != j);
	  k1 = kk_gl + kspan_gl;
	  k2 = 0;
	  do {
	    k2++;
	    a[k1-1] = at_gl[k2-1];
	    b[k1-1] = bt_gl[k2-1];
	    k1 -= inc_gl;
	  } while ( k1 != kk_gl);							
	} while ( jj != 0 );
      } while (j !=1);
    nueve_50:
      j = k3+1;
      nt_gl -= kspnn;
      ii = nt_gl - inc_gl +1;
    } while ( nt_gl >= 0);
    k = k + 0;
  }
}


/**************************************************************************
  Functions for prime factor radix
 ***************************************************************************/

void radix_3(double *a, double *b)
{
  long	k1, k2;
  double  ak, bk, aj, bj;
  
  do {
    do {
      
      k1 = kk_gl + kspan_gl;
      k2 = k1+ kspan_gl;
      ak = a[kk_gl-1];
      bk = b[kk_gl-1];
      aj = a[k1-1] + a[k2-1];
      bj = b[k1-1] + b[k2-1];
      a[kk_gl-1] = ak + aj;
      b[kk_gl-1] = bk + bj;
      ak = -0.5*aj + ak;
      bk = -0.5*bj + bk;
      aj = (a[k1-1]-a[k2-1])*s120_gl;
      bj = (b[k1-1]-b[k2-1])*s120_gl;
      a[k1-1] = ak - bj;
      b[k1-1] = bk + aj;
      a[k2-1] = ak + bj;
      b[k2-1] = bk - aj;
      kk_gl = k2 + kspan_gl;
    } while ( kk_gl < nn_gl);
    kk_gl = kk_gl - nn_gl;
  } while ( kk_gl <= kspan_gl);
}

void radix_5(double *a, double *b)
{
  long k1, k2, k3, k4;
  double ak, aj, bk, bj, akp, akm, ajm, ajp, aa, bkp, bkm, bjm, bjp, bb;
  double c2, s2;
  c2 = c72_gl*c72_gl - s72_gl*s72_gl;
  s2 = 2 * c72_gl * s72_gl;
  do {
    do {
      k1 = kk_gl + kspan_gl;
      k2 = k1 + kspan_gl;
      k3 = k2 + kspan_gl;
      k4 = k3 + kspan_gl;
      akp = a[k1-1] + a[k4-1];
      akm = a[k1-1] - a[k4-1];
      bkp = b[k1-1] + b[k4-1];
      bkm = b[k1-1] - b[k4-1];
      ajp = a[k2-1] + a[k3-1];
      ajm = a[k2-1] - a[k3-1];
      bjp = b[k2-1] + b[k3-1];
      bjm = b[k2-1] - b[k3-1];
      aa  = a[kk_gl-1];
      bb  = b[kk_gl-1];
      a[kk_gl-1] = aa + akp + ajp;
      b[kk_gl-1] = bb + bkp + bjp;
      ak = akp*c72_gl + ajp*c2 + aa;
      bk = bkp*c72_gl + bjp*c2 + bb;
      aj = akm*s72_gl + ajm*s2;
      bj = bkm*s72_gl + bjm*s2;
      a[k1-1] = ak - bj;
      a[k4-1] = ak + bj;
      b[k1-1] = bk + aj;
      b[k4-1] = bk - aj;
      ak = akp*c2 + ajp*c72_gl + aa;
      bk = bkp*c2 + bjp*c72_gl + bb;
      aj = akm*s2 - ajm*s72_gl;
      bj = bkm*s2 - bjm*s72_gl;
      a[k2-1] = ak - bj;
      a[k3-1] = ak + bj;
      b[k2-1] = bk + aj;
      b[k3-1] = bk - aj;
      kk_gl = k4 + kspan_gl;
    } while ( kk_gl < nn_gl);
    kk_gl -= nn_gl;
  } while ( kk_gl <= kspan_gl);
}

void fac_imp(double *a, double *b)
{    
  long k, kspnn, j, k1, k2, jj;
  double ak, bk, aa, bb, aj, bj;
  double c1, s1, c2, s2;
  double ck[MAXF], sk[MAXF];

  k = nfac_gl[i_gl-1];
  kspnn = kspan_gl;
  kspan_gl /= k;
  if (k==3) radix_3(a, b);
  if (k==5) radix_5(a, b);
  if ((k==3) || (k==5)) goto twi;
  if (k!=jf_gl) {
    jf_gl = k;
    s1 = rad_gl/k;
    c1 = cos(s1);
    s1 = sin(s1);
    ck[jf_gl-1] = 1.0;
    sk[jf_gl-1] = 0.0;
    j = 1;
    do {
      ck[j-1] = ck[k-1]*c1 + sk[k-1]*s1;
      sk[j-1] = ck[k-1]*s1 - sk[k-1]*c1;
      k--;
      ck[k-1] =  ck[j-1];
      sk[k-1] = -sk[j-1];
      j++;
    } while ( j<k);
  }
  do {
    do {
      k1 = kk_gl;
      k2 = kk_gl + kspnn;
      aa = a[kk_gl-1];
      bb = b[kk_gl-1];
      ak = aa;
      bk = bb;
      j  = 1;
      k1 = k1 + kspan_gl;
      do {
	k2 -= kspan_gl;
	j++;
	at_gl[j-1] = a[k1-1] + a[k2-1];
	ak += at_gl[j-1];
	bt_gl[j-1] = b[k1-1] + b[k2-1];
	bk += bt_gl[j-1];
	j++;
	at_gl[j-1]  = a[k1-1] - a[k2-1];
	bt_gl[j-1] = b[k1-1] - b[k2-1];
	k1 += kspan_gl;
      } while ( k1 < k2);
      a[kk_gl-1] = ak;
      b[kk_gl-1] = bk;
      k1 = kk_gl;
      k2 = kk_gl + kspnn;
      j = 1;
      do {
	k1 += kspan_gl;
	k2 -= kspan_gl;
	jj = j; 
	ak = aa;
	bk = bb;
	aj = 0.0;
	bj = 0.0;
	k = 1;
	do {
	  k++;
	  ak = at_gl[k-1]*ck[jj-1] + ak;
	  bk = bt_gl[k-1]*ck[jj-1] + bk;
	  k++;
	  aj = at_gl[k-1]*sk[jj-1] + aj;
	  bj = bt_gl[k-1]*sk[jj-1] + bj;
	  jj += j;
	  if (jj>jf_gl)  jj-=jf_gl;
	} while ( k<jf_gl);
	k = jf_gl - j;
	a[k1-1] = ak - bj;
	b[k1-1] = bk + aj;
	a[k2-1] = ak + bj;
	b[k2-1] = bk - aj;
	j++;
      } while (j<k);
      kk_gl += kspnn;
    } while (kk_gl <= nn_gl);
    kk_gl -= nn_gl;
  } while ( kk_gl <= kspan_gl);
  
  /***** Multiply by twiddle factors  *****/
  
 twi:
  if (i_gl==m_gl) flag_gl = 1;
  else {
    kk_gl = jc_gl + 1;
    do {
      c2 = 1.0 - cd_gl;
      s1 = sd_gl;
      do {
	c1  = c2;
	s2  = s1;
	kk_gl += kspan_gl;
	do {
	  do {
	    ak = a[kk_gl-1];
	    a[kk_gl-1] = c2*ak - s2*b[kk_gl-1];
	    b[kk_gl-1] = s2*ak + c2*b[kk_gl-1];
	    kk_gl += kspnn;
	  } while ( kk_gl <= nt_gl);
	  ak = s1 * s2;
	  s2 = s1*c2 + c1*s2;
	  c2 = c1*c2 - ak;
	  kk_gl = kk_gl - nt_gl + kspan_gl;
	} while ( kk_gl <= kspnn);
	c2 = c1 - (cd_gl*c1 + sd_gl*s1);
	s1 = s1 + (sd_gl*c1 - cd_gl*s1);
	
	/***** Compensate for truncation errors *****/
	
	c1 = 0.5/(c2*c2 + s1*s1) + 0.5;
	s1 *= c1;
	c2 *= c1;
	kk_gl  = kk_gl - kspnn + jc_gl;
      } while ( kk_gl <= kspan_gl);		  
      kk_gl = kk_gl - kspan_gl + jc_gl + inc_gl;
    } while ( kk_gl <= (jc_gl+jc_gl));	  
    
  } 
}

void  sing(double *a, double *b, long n, long isn)
{
  if ( n < 2 ) return;
  inc_gl  = isn;
  rad_gl  = TWO_PI;
  s72_gl  = rad_gl / 5.0;
  c72_gl  = cos(s72_gl);
  s72_gl  = sin(s72_gl);
  s120_gl = sqrt(0.75);
  if (isn < 0) { 
    s72_gl  = -s72_gl;
    s120_gl = -s120_gl;
    rad_gl  = -rad_gl;
    inc_gl  = -inc_gl;
  } 
  nt_gl    = inc_gl * n;
  ks_gl    = inc_gl * n;
  kspan_gl = ks_gl;
  nn_gl    = nt_gl - inc_gl;
  jc_gl    = ks_gl / n;
  radf_gl  = rad_gl * jc_gl * 0.5;
  i_gl     = 0;
  jf_gl    = 0;
  flag_gl  = 0;
  fac_des (n );  
  do   {
    sd_gl = radf_gl / kspan_gl;
    cd_gl = 2.0 * sin(sd_gl) * sin(sd_gl);
    sd_gl = sin(sd_gl+sd_gl);
    kk_gl = 1;
    i_gl  = i_gl + 1;
    if (nfac_gl[i_gl-1]==2)
      radix_2( a, b);
    if (nfac_gl[i_gl-1]==4) 
      radix_4(isn, a, b); 
    if ( (nfac_gl[i_gl-1]!=2) && (nfac_gl[i_gl-1]!=4)) 
      fac_imp(a, b);	
  } while ( flag_gl != 1);
  permute(n, a, b);
} 

/* Calculates the the Fourier transform of 2*half_length real values.
   Original data values are stored alternately in arrays a and b.
   The cosine coefficients are in a[0], a[1] .........a[half_length} and
   the sine coefficients are in b[0], b[1] .........b[half_length}.
   The coeffcients must be scaled by 1/(4*half_length) in the calling
   procedure.  */

/* April/1/92 Tried Singleton's subroutine and it does not seem to work.
   I am modifying it and folowing the procedure of Cooley, Lewis and Welch
   J. Sound Vib., vol. 12, pp 315-337, July  1970. Will extend the 
   procedure so half_length can be odd also. */	

 /*	Assume we have the transform A(n) of x(even) + i x(odd)
	1-)  A1(n) = (1/2) [ Ac(-n) + A(n)] =(1/2) [Ac(N-n) + A(n)]  
		(Ac = complex of A)    (for N even or odd)
	2-)  A2(n) = (i/2)[Ac(-n)-A(n)] =(i/2) [Ac(N-n)-A(n)]  
		(for N even or odd)
	3-)  C(n) = (1/2)[A1(n) + A2(n)*W2N**(-n)]   (1)
		0,1,2,3..... N   N even
		0,1,2,3..... N-1 N odd						    
		Use the simmetry of the A1 and A2 sequences
		C(N-n) = (1/2) [A1(N-n) + A2(N-n)*W2N**(-N+n)] 
			  = (1/2) [A1c(n) - A2c(n)*W2N**(n)]   (2)
		Evaluate (1) and (2) for n=0,1,2,...N/2 -1 and (1) for N/2
		if N is even  ( (2) is also good
		Evaluate (1) for n =0 and
			    (1) and  (2) for n=1,2,...(N-1)/2  
		if N is odd
				
		Let the factors of two be taken care in the normalization
		outside this procedure, ie, the coefficients will be
		four times larger. */							    

/* 4/april/1992 Everything is working fine. Singleton's Realtr
   works ok. */

void  realtr(double *a, double *b, long half_length, int isn)
{
  long  nh, j, k;
  double sd, cd, sn, cn, a1r, a1i, a2r, a2i, re, im;
  nh = half_length >> 1;    /* Should work for even and odd */
  sd = M_PI_2 /half_length;
  cd = 2.0 * sin(sd) * sin(sd);
  sd = sin(sd+sd);	
  sn = 0;
  if ( isn <0) {
    cn = 1 ;
    a[half_length] = a[0];
    b[half_length] = b[0];  
  }
  else {
    cn = -1;
    sd = -sd;
  }
  
  for (j=0; j <= nh; j++) {
    k = half_length-j;
    a1r = a[j] + a[k];
    a2r = b[j] + b[k];
    a1i = b[j] - b[k];
    a2i = -a[j] + a[k];
    re = cn*a2r + sn*a2i;
    im =-sn*a2r + cn*a2i; 
    b[k] = +im - a1i;  
    b[j] = im + a1i;
    a[k] = a1r - re;
    a[j] = a1r + re;
    a1r = cn - (cd*cn + sd*sn);
    sn = (sd*cn - cd*sn) + sn;
    /* compensate for truncation error */
    cn = 0.5/(a1r*a1r + sn*sn) + 0.5;
    sn *= cn;
    cn *= a1r;
  }
}

#undef   TWO_PI
#undef   MAXF
#undef   MAXP
