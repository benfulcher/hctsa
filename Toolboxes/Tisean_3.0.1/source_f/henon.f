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
c   henon.f
c   iterate Henon map
c   author Thomas Schreiber (1998)
c===========================================================================
      double precision xo, yo, xn, yn, a, b
      character*72 fout
      data a/1.4/,b/0.3/,ntrans/10000/,xo/.68587/,yo/.65876/
      data iverb/1/

      call whatido("Henon map",iverb)
      nmax=imust('l')
      ntrans=ican('x',ntrans)
      a=fcan('A',real(a))
      b=fcan('B',real(b))
      xo=fcan('X',real(xo))
      yo=fcan('Y',real(yo))
      isout=igetout(fout,iverb)

      if(isout.eq.1) fout="henon.dat"
      call outfile(fout,iunit,iverb)
      n=-ntrans
 1    n=n+1
      xn=1.-a*xo**2+b*yo
      yn=xo
      xo=xn
      yo=yn
      if(n.lt.1) goto 1
      write(iunit,*) real(xn), real(yn)
      if(nmax.eq.0.or.n.lt.nmax) goto 1
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-l# [-A# -B# -X# -Y# -o outfile -x# -V# -h]")
      call popt("l","number of points x,y (l=0: infinite)")
      call popt("A","parameter a (1.4)")
      call popt("B","parameter b (0.3)")
      call popt("X","initial x")
      call popt("Y","initial y")
      call popt("x","number of transients discarded (10000)")
      call pout("henon.dat")
      call pall()
      stop
      end


