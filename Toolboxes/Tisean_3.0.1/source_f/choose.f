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
c   choose.f
c   Choose columns and sub-sequences from a file
c   author T. Schreiber (1999)
c===========================================================================

      parameter(nx=1000000,mx=5)
      dimension x(nx,mx), icol(mx)
      character*72 file, fout
      data iverb/15/

      call whatido("Choose columns and sub-sequences from a file",iverb)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      mcmax=ican("m",0)
      call columns(mc,mx,icol)
      if(mcmax.eq.0) mcmax=max(1,mc)
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call xreadfile(nmax,mcmax,nx,x,nexcl,icol,file,iverb)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_select")
      call outfile(fout,iunit,iverb)
      call xwritefile(nmax,mcmax,nx,x,fout,iverb)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-o outfile -l# -x# -m# -c#[,#] -V# -h] file")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("m","number of columns to be read (1)")
      call popt("c","columns to be read (1)")
      call pout("file_select")
      call pall()
      stop
      end
