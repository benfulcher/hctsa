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
c   box assisted sorting/ranking utilities 
c   author T. Schreiber (1998) based on earlier versions
c===========================================================================
      subroutine rank(nmax,x,list)
c  rank points in x
      parameter(nptr=100000)
      dimension x(nmax), list(nmax), jptr(0:nptr)

      call minmax(nmax,x,xmin,xmax)
      if(xmin.eq.xmax) then
         do 10 n=1,nmax
 10         list(n)=n
         return
      endif
      nl=min(nptr,nmax/2)
      sc=(nl-1)/(xmax-xmin)
      do 20 i=0,nl
 20      jptr(i)=0
      do 30 n=1,nmax
         xn=x(n)
         i=int((xn-xmin)*sc)
         ip=jptr(i)
         if ((ip.eq.0).or.(xn.le.x(ip))) then
            jptr(i)=n
         else
 1          ipp=ip
            ip=list(ip)
            if ((ip.gt.0).and.(xn.gt.x(ip))) goto 1
            list(ipp)=n
         endif
 30      list(n)=ip
      n=0
      do 40 i=0,nl
         ip=jptr(i)
 2       if (ip.eq.0) goto 40
         n=n+1
         ipp=ip
         ip=list(ip)
         list(ipp)=n
         goto 2
40       continue
      end

      subroutine indexx(nmax,x,list)
c make index table using rank
      dimension x(nmax), list(nmax)
      
      call rank(nmax,x,list)
      call rank2index(nmax,list)
      end

      subroutine rank2index(nmax,list)
c converts a list of ranks into an index table (or vice versa) in place
      integer list(nmax)

      do 10 n=1,nmax
 10      list(n)=-list(n)
      do 20 n=1,nmax
         if(list(n).gt.0) goto 20               ! has been put in place already
         ib=n
         im=-list(n)
 1       it=-list(im)
         list(im)=ib
         if(it.ne.n) then
            ib=im
            im=it
            goto 1
         else
            list(n)=im
         endif
 20      continue
      end

      subroutine sort(nmax,x,list)
c sort using rank and rank2sort
      dimension x(nmax), list(nmax)

      call rank(nmax,x,list)
      call rank2sort(nmax,x,list)
      end

      subroutine rank2sort(nmax,x,list)
c sort x using list of ranks
      dimension x(nmax), list(nmax)

      do 10 n=1,nmax
 10      list(n)=-list(n)
      do 20 n=1,nmax
         if(list(n).gt.0) goto 20               ! has been put in place already
         ib=n
         hb=x(n)
 1       it=-list(ib)
         list(ib)=it
         ht=x(it)
         x(it)=hb
         if(it.ne.n) then
            ib=it
            hb=ht
            goto 1
         endif
 20      continue
      end

      subroutine index2sort(nmax,x,list)
c sort x using list of indices
      dimension x(nmax), list(nmax)

      do 10 n=1,nmax
 10      list(n)=-list(n)
      do 20 n=1,nmax
         if(list(n).gt.0) goto 20               ! has been put in place already
         ib=n
         h=x(n)
 1       it=-list(ib)
         list(ib)=it
         if(it.ne.n) then
            x(ib)=x(it)
            ib=it
            goto 1
         else
            x(ib)=h
         endif
 20      continue
      end

      function which(nmax,x,k,list)
      dimension x(nmax), list(nmax)

      call indexx(nmax,x,list)
      which=x(list(k))
      end

