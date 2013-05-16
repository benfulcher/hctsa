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
c   compare.f
c   compare two data sets
c   author T. Schreiber
c===========================================================================

      parameter(nx=1000000,mx=2)
      character*72 file
      dimension x(nx,mx), icol(mx)
      data iverb/1/

      call whatido("compare time series in RMS sense",iverb)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      call columns(mc,mx,icol)
      mcmax=mx
      if(nstrings().ne.1) call usage()
      call nthstring(1,file)

      nmax=nmaxx
      call xreadfile(nmax,mcmax,nx,x,nexcl,icol,file,iverb)
      if(file.eq."-") file="stdin"

      call rms(nmax,x(1,1),sc1,sd1)
      call rms(nmax,x(1,2),sc2,sd2)
      do 10 n=1,nmax
 10      x(n,1)=x(n,2)-x(n,1)
      call rms(nmax,x(1,1),scd,sdd)

      write(istderr(),*)
      write(istderr(),*) "col ", icol(1), ": Mean ", sc1, 
     .   ", standard deviation ", sd1
      write(istderr(),*) "col ", icol(2), ": Mean ", sc2, 
     .   ", standard deviation ", sd2
      write(istderr(),*)
      write(istderr(),*) "mean difference              ", scd 
      write(istderr(),*)  
     .   "root mean squared difference ", sqrt(sdd**2+scd**2) 
      write(istderr(),*) "standard deviation           ", sdd
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-l# -x# -c#[,#] -V# -h] file")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","columns to be read (1,2)")
      call pall()
      stop
      end
