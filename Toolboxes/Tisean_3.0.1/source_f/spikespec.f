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
c   power spectrum of spike trains
c   author T. Schreiber (1999)
c===========================================================================
      parameter(nx=1000000)
      dimension x(nx), lx(nx), sp(nx)
      character*72 file, fout
      data iverb/1/

      call whatido("power spectrum of spike trains",iverb)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      fmax=fcan("F",0)
      nfreq=min(ican("#",0),nx)
      fres=fcan("w",0)
      inter=lopt("i",1)
      isout=igetout(fout,iverb)

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         if(inter.eq.0) goto 1
         do 20 n=2,nmax
 20         x(n)=x(n)+x(n-1)
 1       call sort(nmax,x,lx)
         if(fmax.le.0.) fmax=2*nmax/(x(nmax)-x(1))
         if(nfreq.le.0) nfreq=fmax*(x(nmax)-x(1))/2
         write(istderr(),*) "spikespec: total time covered: ", 
     .      x(nmax)-x(1)
         write(istderr(),*) "spikespec: computing ", nfreq, 
     .      " frequencies up to ", fmax
         do 30 n=1,nfreq
            f=(n*fmax)/nfreq
 30         sp(n)=sspect(nmax,x,f)
         ibin=nfreq*fres/2
         if(ibin.gt.0) write(istderr(),*) 
     .      "spikespec: binning", 2*ibin+1," frequencies"
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_ss")
         call outfile(fout,iunit,iverb)
         do 40 n=1+ibin,nfreq-ibin,2*ibin+1
            f=(n*fmax)/nfreq
            p=0
            do 50 ib=n-ibin,n+ibin
 50            p=p+sp(ib)
 40         write(iunit,*) f, p
            if(iunit.eq.istdout()) write(iunit,*)
            if(iunit.eq.istdout()) write(iunit,*)
 10      if(iunit.ne.istdout()) close(iunit)
      end

      function sspect(nmax,x,f)
      dimension x(nmax)
      data pi/3.1415926/

      omega=2*pi*f
      sr=0
      si=0
      do 10 n=1,nmax
         sr=sr+cos(omega*x(n))
 10      si=si+sin(omega*x(n))
      sspect=sr**2+si**2
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-F# -## -w# -i -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("F","maximal frequency [2*l / total time]")
      call popt("#","number of frequencies [F* total time /2]")
      call popt("w","frequency resolution [0]")
      call popt("i","input data: intervals rather than times")
      call popt("l","number of values to be read [all]")
      call popt("x","number of values to be skipped [0]")
      call popt("c","column to be read [1 or file,#]")
      call pout("file_ss")
      call pall()
      stop
      end



