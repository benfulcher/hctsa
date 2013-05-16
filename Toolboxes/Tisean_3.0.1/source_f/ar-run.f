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
c   ar-run.f
c   iterate AR model, e.g. as fitted by ar-model (Dresden)
c   author T. Schreiber (1999)
c===========================================================================

      parameter(npmax=100)
      character*72 file, fout, fline
      dimension x(-npmax:npmax), a(npmax)
      external rand
      data np/npmax/, ntrans/10000/, iuni/0/
      data iverb/1/

      call whatido("iterate AR model, e.g. as fitted by ar-model",iverb)
      np=ican("p",np)
      if(np.gt.npmax) stop "ar-run: make npmax larger."
      nmax=imust('l')
      ntrans=ican("x",ntrans)
      if(lopt("u",1).eq.1) iuni=1
      r=rand(sqrt(abs(fcan("I",0.0))))
      isout=igetout(fout,iverb)

      do 10 n=1,npmax
         x(-n)=0.
 10      x(n)=0.
      call nthstring(1,file)
      call infile(file,iunit,iverb)
      read(iunit,'(a)') fline
      if(fline(1:1).eq."#") then
         read(fline(18:72),'(f20.0)',err=999) var
         do 20 j=1,np
            read(iunit,'(a1,f20.0)',err=999) fline(1:1), a(j)
 20         if(fline(1:1).ne."#") goto 1
      else
         read(fline(1:72),'(f20.0)',err=999) var
         do 30 j=1,np
 30         read(iunit,'(f20.0)',err=999,end=1) a(j)
      endif
 1    np=j-1
      if(iv_echo(iverb).eq.1) then
         write(istderr(),*) 'coefficients:      ', (a(i),i=1,np)
         write(istderr(),*) 'driving amplitude: ', var
      endif
      if(isout.eq.1) fout="ar.dat"
      call outfile(fout,iunit,iverb)
      n=-ntrans
 2    n=n+1
      nn=mod(n+ntrans,np)+1
      xx=rgauss(0.0,var)
      do 40 j=1,np
 40      xx=xx+a(j)*x(nn-j)
      x(nn)=xx
      x(nn-np)=xx
      if(n.lt.1) goto 2
      write(iunit,*) xx
      if(nmax.eq.0.or.n.lt.nmax) goto 2
      stop

 999  write(istderr(),'(a)') "wrong input format! try:"
      write(istderr(),'(a)') "(rms of increments)"
      write(istderr(),'(a)') "a(1)"
      write(istderr(),'(a)') "a(2)"
      write(istderr(),'(a)') "..."
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-l# [-p# -I# -o outfile -x# -V# -h] file")
      call popt("l","number of iterations (l=0: infinite)")
      call popt("p","order of AR-model (default determined from input)")
      call popt("I","seed for random numbers")
      call popt("x","number of transients discarded (10000)")
      call pout("ar.dat")
      call pall()
      stop
      end

