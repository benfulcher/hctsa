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
c   autocor.f
c   autocorrelation function through FFT
c   author T. Schreiber (1998), H. Kantz (2007)
c===========================================================================

      parameter(nx=1000000)
      dimension x(2*nx)
      character*72 file, fout
      data iverb/1/

      call whatido("autocorrelation function estimated by FFT",iverb)
      ivar=lopt('v',1)
      iper=lopt('p',1)
      iexact=lopt('P',1)
      if(iexact.ne.0) iper=1
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         if(file.eq."-") file="stdin"
         if(ivar.eq.0) call normal(nmax,x,sc,sd)
         if(iper.eq.0) then
            nmaxp=nmore(2*nmax)
            do 20 n=nmax+1,nmaxp
 20            x(n)=0.
            call store_spec(nmaxp,x,1)
            do 30 n=1,nmax
 30            x(n)=x(n)/real(nmax-n+1)
         else
            nmaxp=nmax
            if(iexact.eq.0) then
               nmaxp=nless(nmax)
               if(nmaxp.ne.nmax.and.iv_io(iverb).eq.1)  
     .            write(istderr(),*) "autocor: using", nmaxp
            endif
            call store_spec(nmaxp,x,1)
            do 50 n=1,nmaxp
 50            x(n)=x(n)/real(nmaxp)
         endif
         if(isout.eq.1) call addsuff(fout,file,"_co")
         call outfile(fout,iunit,iverb)
         if(ivar.eq.0) then
            if(sd.eq.0) stop "autocor: cannot normalise - zero variance"
            fsc=1./x(1)
         else
            fsc=1.
         endif
         do 60 n=1,min(nmax,nmaxp)
 60         write(iunit,*) n-1, fsc*x(n)
 10      if(iunit.ne.istdout()) close(iunit)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-v -p -P -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("v","give unnormalised autocovariance")
      call popt("p","assume periodic continuation")
      call popt("P","assume periodic continuation exactly")
      call popt("l","number of values to be read [all]")
      call popt("x","number of values to be skipped [0]")
      call popt("c","column to be read [1 or file,#]")
      call pout("file_co")
      call pall()
      stop
      end



