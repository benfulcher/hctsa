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
c   Wiener filter (2): filter using periodogram in file
c   author T. Schreiber (1998)
c===========================================================================
      parameter(nx=1000000)
      dimension x(nx), a(nx)
      character*72 file, fout, ffin
      data h/1./, dh/0./
      data iverb/1/

      call whatido("Wiener filter (second part)",iverb)
      h=fcan("f",h)
      dh=fcan("w",dh)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isfin=igetfin(ffin,iverb)
      isout=igetout(fout,iverb)
      call nthstring(1,file)
      if(file.eq."-") stop "wiener2: cannot read stdin"
      call readfile(nmax,x,nexcl,jcol,file,iverb)
      call normal(nmax,x,sc,sd)
      nmaxp=nmore(nmax)

      if(dh.eq.0.) dh=h/nmaxp
      ibin=nmaxp*dh/(2*h)
      if(ibin.gt.0) write(istderr(),*) 
     .   "wiener1: binning", 2*ibin+1," frequencies"
      if(isout.eq.1) call addsuff(fout,file,"_amp")
      if(fout.eq." ") fout="-"
      call infile(fout,iunit,iverb)
      read(iunit,*) dum, a(1)
      do 10 i=2+ibin,(nmaxp+1)/2-ibin,2*ibin+1
 10      read(iunit,*) dum, a(2*i-2)
      if(mod(nmaxp,2).eq.0)  read(iunit,*) 
     .   dum, a(nmaxp)
      d=tospec(nmaxp,a,x,ibin)
      if(iv_io(iverb).eq.1) write(istderr(),*) "rms correction: ", d
      if(isfin.eq.1) call addsuff(ffin,file,"_wc")
      do 20 n=1,nmax
 20      x(n)=x(n)+sc
      call writefile(nmax,x,ffin,iverb)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-f# -w# -o outfile1 -O outfile -l# -x# -c# -V# -h] file")
      call ptext("to provide periodogram, first run:"//
     .   " wiener1 [-f# -w# -o outfile -l# -x# -c# -V# -h] file")
      call ptext("make sure -f# -w# are the same in both wiener calls")
      call popt("f","sampling rate (e.g. in Hz, default 1.)")
      call popt("w","frequency resolution (e.g. in Hz, default 1/N)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call popt("o","output file of wiener1, just -o means file_amp")
      call popt("O","final output file name, just -O means file_wc")
      call pall()
      call ptext("Note: ""-"" not accepted as file")
      write(istderr(),'()') 
      stop
      end

      function igetfin(fout,iverb)
c gets alternate output file name, default " "
c return 1 if fout must be determined from input file name
      character*(*) fout

      igetfin=0
      call stcan("O",fout," ")
      if(fout.ne." ") return
      igetfin=lopt("O",1)
      end
