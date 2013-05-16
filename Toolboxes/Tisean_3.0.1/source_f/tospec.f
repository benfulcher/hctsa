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
c   Wiener filter x according to spectrum a
c   (routine needed by surrogates.f)
C   author Thomas Schreiber (1998)
c===========================================================================
      function tospec(nmax,a,x,ibin)
      parameter(nx=1000000)
      dimension x(nmax), a(nmax), w(nx), w1(nx), w2(nx), iw(15)
      save w2, iw

      if(nmax.gt.nx) stop "tospec: make nx larger."
      do 10 n=1,nmax
 10      w(n)=x(n)
      call rffti1(nmax,w2,iw)  
      call rfftf1(nmax,x,w1,w2,iw)
      do 20 n=1,nmax
 20      x(n)=x(n)/real(nmax)
      x(1)=x(1)*(a(1)/x(1)**2)
      do 30 i=2+ibin,(nmax+1)/2-ibin,2*ibin+1
         p=0
         do 40 ib=i-ibin,i+ibin
 40         p=p+x(2*ib-2)**2+x(2*ib-1)**2
         ab=a(2*i-2)/p
         do 30 ib=i-ibin,i+ibin
            x(2*ib-2)=x(2*ib-2)*ab
 30         x(2*ib-1)=x(2*ib-1)*ab
      if(mod(nmax,2).eq.0) x(nmax)=x(nmax)*(a(nmax)/x(nmax)**2)
      call rfftb1(nmax,x,w1,w2,iw)
      tospec=0
      do 50 n=1,nmax
 50      tospec=tospec+(x(n)-w(n))**2
      tospec=sqrt(tospec/nmax)
      end
