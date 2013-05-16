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
c   part of the TISEAN randomize package for constraint surrogates
c   permutation scheme that swaps to randomly chosen data points
c   this may also be used as a template for your own attempts
c   author T. Schreiber (1999)
c
c-------------------------------------------------------------------
c get permutation specific options
c
      subroutine opts_permute()
      parameter(nx=100000)
      dimension nxclu(nx)
      character*80 filex
      common /permutecom/ mxclu, nxclu

      call stcan('X',filex,' ')
      mxclu=0
      if(filex.eq." ") return
      open(10,file=filex,status="old",err=999)
 1    read(10,*,err=999,end=998) nn
      mxclu=mxclu+1
      nxclu(mxclu)=nn
      goto 1
 998  return
 999  write(istderr(),'(a)') "permute: cannot open "//filex
      stop
      end

c-------------------------------------------------------------------
c print version information on permutation scheme
c
      subroutine what_permute()
      call ptext("Permutation scheme: random pairs")
      end

c-------------------------------------------------------------------
c print permutation specific usage message
c
      subroutine usage_permute()
      call ptext("Permutation options: [-X xfile]")
      call popt("X", "list of indices excluded from permutation")
      end

c-------------------------------------------------------------------
c initialise all that is needed for permutation scheme 
c
      subroutine permute_init()
      parameter(nx=100000)
      dimension x(nx)
      common nmax,cost,temp,cmin,rate,x
      
      if(nmax.gt.nx) stop "permute: make nx larger."
      do 10 i=1,nmax
         call permute(n1,n2)
 10      call exch(n1,n2)
      end

c-------------------------------------------------------------------
c find two indices n1, n2 to be exchanged, maybe using a parameter 
c par provided by the cooling schedule
c
      subroutine permute(n1,n2)
      parameter(nx=100000)
      dimension nxclu(nx)
      common /permutecom/ mxclu, nxclu
      common nmax
      external rand

 1    n1=min(int(rand(0.0)*nmax)+1,nmax)
      do 10 n=1,mxclu
 10      if(n1.eq.nxclu(n)) goto 1
 2    n2=min(int(rand(0.0)*nmax)+1,nmax)
      if(n2.eq.n1) goto 2
      do 20 n=1,mxclu
 20      if(n2.eq.nxclu(n)) goto 2
      end

c-------------------------------------------------------------------
c given two indices n1, n2, actually perform the exchange
c
      subroutine exch(n1,n2)
      parameter(nx=100000)
      dimension x(nx)
      common nmax,cost,temp,cmin,rate,x

      h=x(n1)
      x(n1)=x(n2)
      x(n2)=h
      end
