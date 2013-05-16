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
c   events.f
c   convert inter-event intervals to event times
c   author T. Schreiber (1999)
c===========================================================================
      parameter(nx=1000000)
      dimension x(nx)
      character*72 file, fout
      data iverb/1/

      call whatido("interval to event time conversion",iverb)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         nmax=nmax+1
         do 20 n=nmax,2,-1
 20         x(n)=x(n-1)
         x(1)=0
         do 30 n=2,nmax
 30         x(n)=x(n)+x(n-1)
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_st")
 10      call writefile(nmax,x,fout,iverb)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_st")
      call pall()
      stop
      end


