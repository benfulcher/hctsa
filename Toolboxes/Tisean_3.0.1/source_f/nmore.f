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
c   utilities for TISEAN f-sources
c
      function nmore(n)
c find smallest factorisable number .ge.n

      nmore=n
 1    if(isfact(nmore).eq.1) return
      nmore=nmore+1
      goto 1
      end

      function nless(n)
c find largest factorisable number .le.n

      nless=n
 1    if(isfact(nless).eq.1) return
      nless=nless-1
      goto 1
      end

      function isfact(n)
c determine if n is factorisable using the first nprimes primes
      parameter(nprimes=3)
      dimension iprime(nprimes)
      data iprime/2,3,5/

      isfact=1
      ncur=n
 1    if(ncur.eq.1) return
      do 10 i=1,nprimes
         if(mod(ncur,iprime(i)).eq.0) then
            ncur=ncur/iprime(i)
            goto 1
         endif
 10      continue
      isfact=0
      end
