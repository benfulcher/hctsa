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
c   Fourier power spectrum
c   author T. Schreiber (1998) and earlier
c   modified by H. Kantz 2007
c=========================================================================== 
      parameter(nx=1000000)
      dimension x(nx)
      character*72 file, fout
      data h/1./, dh/0./
      data iverb/1/

      call whatido("Power spectrum by FFT",iverb)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      h=fcan("f",h)
      dh=fcan("w",dh)
      isout=igetout(fout,iverb)

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_sp")
         nmaxp=nless(nmax)
         if(nmaxp.ne.nmax) 
     .      write(istderr(),*) "spectrum: using first ", nmaxp
         if(dh.eq.0.) dh=h/nmaxp
         ibin=nmaxp*dh/(2*h)
         if(ibin.gt.0) write(istderr(),*) 
     .      "spectrum: binning", 2*ibin+1," frequencies"
         call store_spec(nmaxp,x,0)
         call outfile(fout,iunit,iverb)
         write(iunit,*) 0., x(1)
         do 20 i=2+ibin,nmaxp/2+1-ibin,2*ibin+1
            p=0
            do 30 ib=i-ibin,i+ibin
 30            p=p+x(2*ib-2)
 20         write(iunit,*) h*(i-1)/real(nmaxp), p
            if(iunit.eq.istdout()) write(iunit,*)
            if(iunit.eq.istdout()) write(iunit,*)
 10         if(iunit.ne.istdout()) close(iunit)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-f# -w# -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("f","sampling rate, e.g. in Hz [default 1.]")
      call popt("w","frequency resolution, e.g. in Hz, [default 1/N]")
      call popt("l","number of values to be read [all]")
      call popt("x","number of values to be skipped [0]")
      call popt("c","column to be read, [1] or file,#")
      call pout("file_sp")
      call pall()
      stop
      end
