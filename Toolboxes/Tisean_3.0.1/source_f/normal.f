c===========================================================================
c
c   This file is part of TISEAN
c 
c   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
c 
c   TISEAN is free software; you can redistribute it and/or modify
c   it under the terms of the GNU General Public License as published by
c   the Free Software Foundation; either version 2 of the License, or
c   (at your option) any later version.
c
c   TISEAN is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with TISEAN; if not, write to the Free Software
c   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
c
c===========================================================================
c   utilities for normalisation of time series
c   author T. Schreiber (1998)
c===========================================================================
      subroutine rms(nmax,x,sc,sd)
c  return mean sc and rms amplitude sd
      dimension x(nmax)

      sc=0.
      do 10 n=1,nmax
 10      sc=sc+x(n)
      sc=sc/nmax
      sd=0.
      do 20 n=1,nmax
 20      sd=sd+(x(n)-sc)**2
      sd=sqrt(sd/nmax)
      end

      subroutine normal(nmax,x,sc,sd)
c  subtract mean, return mean sc and rms amplitude sd
      dimension x(nmax)

      call rms(nmax,x,sc,sd)
      do 10 n=1,nmax
 10      x(n)=x(n)-sc
      end

      subroutine normal1(nmax,x,sc,sd)
c  subtract mean, rescale to unit variance, 
c  return mean sc and rms amplitude sd
      dimension x(nmax)

      call rms(nmax,x,sc,sd)
      if(sd.eq.0.) stop 
     .   "normal1: zero variance, cannot normalise"
      do 10 n=1,nmax
 10      x(n)=(x(n)-sc)/sd
      end

      subroutine minmax(nmax,x,xmin,xmax)
c  obtain smallest and  largest value in x
      dimension x(nmax)

      xmin=x(1)
      xmax=x(1)
      do 10 n=2,nmax
         xmin=min(x(n),xmin)
 10      xmax=max(x(n),xmax)
      end

