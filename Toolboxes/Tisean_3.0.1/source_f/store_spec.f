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
c   store data periodogram of x 
c   if iback.ne.0 transform back to get autocorrelation instead
C   author Thomas Schreiber (1998)
c===========================================================================
      subroutine store_spec(nmax,x,iback)
      parameter(nx=1000000)
      dimension x(nmax), w1(nx), w2(nx), iw(15)
      save w2, iw

      if(nmax.gt.nx) stop "store_spec: make nx larger."
      call rffti1(nmax,w2,iw)  
      call rfftf1(nmax,x,w1,w2,iw)
      do 10 n=1,nmax
 10      x(n)=x(n)/real(nmax)
      x(1)=x(1)**2
      do 20 n=2,(nmax+1)/2
         amp=x(2*n-2)**2+x(2*n-1)**2
         pha=atan2(x(2*n-1),x(2*n-2))
         x(2*n-2)=amp
 20      x(2*n-1)=pha
      if(mod(nmax,2).eq.0) x(nmax)=x(nmax)**2
      if(iback.eq.0) return
      do 30 n=1,nmax
 30      x(n)=x(n)*nmax
      do 40 n=2,(nmax+1)/2
 40      x(2*n-1)=0
      call rfftb1(nmax,x,w1,w2,iw)
      end
