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
c   notch filter in the time domain
c   author T. Schreiber
c===========================================================================

      parameter(nx=1000000)
      dimension x(nx), y(nx)
      character*72 file, fout
      data h/1./, w/0.01/, pi/3.1415926/
      data iverb/1/

      call whatido("notch filter",iverb)
      f=fmust("X")
      h=fcan("f",h)
      w=fcan("w",w)
      fw=tan(pi*f/h)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         if(file.eq."-") file="stdin"
         d=fnotch(nmax,x,y,fw,w)
         if(isout.eq.1) call addsuff(fout,file,"_notch")
 10      call writefile(nmax,y,fout,iverb)
      end

      function fnotch(nmax,x,y,fw,w)
      dimension x(nmax), y(nmax)

      a=(1+w*fw)**2+fw**2
      c0=   (1+fw**2)/a
      c1=-2*(1-fw**2)/a
      c2=c0
      d1= 2*(1-w**2*fw**2-fw**2)/a
      d2=  -((1-w*fw)**2+fw**2)/a

      y(1)=c0*x(1)
      y(2)=c0*x(2)+c1*x(1)+d1*y(1)
      do 10 n=3,nmax
 10      y(n)=c0*x(n)+c1*x(n-1)+c2*x(n-2)+d1*y(n-1)+d2*y(n-2)
      fnotch=0
      do 20 n=1,nmax
 20      fnotch=fnotch+(x(n)-y(n))**2
      fnotch=sqrt(fnotch/nmax)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-X# [-f# -w# -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("X","frequency to be cancelled")
      call popt("f","sampling rate of data (1)")
      call popt("w","width of filter (f/100)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_notch")
      call pall()
      stop
      end
