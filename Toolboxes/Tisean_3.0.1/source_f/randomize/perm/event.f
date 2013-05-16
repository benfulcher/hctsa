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
C   permutation scheme for event times
c   one event time is changed such that the two adjacent intervals are swapped
c   author T. Schreiber (1999)
c
c-------------------------------------------------------------------
c get permutation specific options
c
      subroutine opts_permute()
      end

c-------------------------------------------------------------------
c print version information on permutation scheme
c
      subroutine what_permute()
      call ptext("Permutation scheme: event time preserving intervals")
      end

c-------------------------------------------------------------------
c print permutation specific usage message
c
      subroutine usage_permute()
      end

c-------------------------------------------------------------------
c initialise all that is needed for permutation scheme 
c
      subroutine permute_init()
      parameter(nx=100000) 
      dimension x(nx)
      common nmax,cost,temp,cmin,rate,x

      do 10 n=1,nmax*log(nmax*1.)
         call permute(n1,n2)
 10      call exch(n1,n2)
      end

c-------------------------------------------------------------------
c find two indices n1, n2 to be exchanged, maybe using a parameter 
c par provided by the cooling schedule
c
c here, n2 is not used at all; event 1 and nmax are never changed
c
      subroutine permute(n1,n2)
      common nmax
      external rand

      n1=min(int(rand(0.0)*nmax)+2,nmax-1)
      end

c-------------------------------------------------------------------
c given two indices n1, n2, actually perform the exchange
c
      subroutine exch(n1,n2)
      parameter(nx=100000)
      dimension x(nx)
      common nmax,cost,temp,cmin,rate,x

      x(n1)=x(n1-1)+x(n1+1)-x(n1)
      end
