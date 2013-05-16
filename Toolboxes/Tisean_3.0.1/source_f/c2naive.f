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
c   c2naive.f
c   correlation integral c2, no fast neighbour search
c   complete direct neighbour search
c   author T. Schreiber (1998)
c   modified H. Kantz, Feb. 2007
c===========================================================================

      parameter(nx=1000000,me=30,meps=800)
      dimension x(nx), c(0:meps,me)
c j=a*(log(d)-log(xmax-xmin))    d=xmax-xmin -> j=0
c a=-rs/log(2.)                  d=s*2**(-j/res)
      character*72 file, fout
      data res/2./
      data iverb/1/

      call whatido("correlation sum, complete naive neighbour search,
     .univariate data only"
     .             ,iverb)
      call whatido("univariate data only",iverb)
      id=imust("d")
      mmax=imust("M")
      mmin=ican("m",1)
      ntmin=imust("t")
      ntmax=ican("T",nx)
      res=fcan("#",res)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)
      if(fout.eq." ") isout=1

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_c2")
         call minmax(nmax,x,xmin,xmax)
         sc=xmax-xmin
         a=-res/log(2.)
         do 20 m=mmin,mmax
            do 20 j=0,meps
 20            c(j,m)=0.
         call d2naive(nmax,x,id,mmin,mmax,c,meps,log(sc),a,ntmin,ntmax)
         call outfile(fout,iunit,iverb)
         do 30 m=mmin,mmax
            write(iunit,'(4h#m= ,i5)') m
            do 40 j=meps-1,0,-1
 40            c(j,m)=c(j,m)+c(j+1,m)
            do 50 j=0,meps
               if(c(j,m).eq.0.) goto 1
 50            write(iunit,*) sc*2**(-j/res), c(j,m)/c(0,m)
 1          write(iunit,'()')
 30         write(iunit,'()')
         close(iunit)
 10      continue
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# -M# -t# [-m# -##"//
     .   " -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("d","delay")
      call popt("M","maximal embedding dimension")
      call popt("t","minimal time separation")
      call popt("m","minimal embedding dimension [1]")
      call popt("#","resolution, values per octave [2]")
c      call popt("T","for Guido")
      call popt("l","number of values to be read [all]")
      call popt("x","number of values to be skipped [0]")
      call popt("c","column to be read [1 or file,#]")
      call pout("file_c2")
      call pall()
      stop
      end

      subroutine d2naive(nmax,x,id,mmin,mmax,c,meps,scl,a,ntmin,ntmax)
      parameter(nx=1000000,tiny=1e-30)
      dimension x(nmax),c(0:meps,mmax),d(nx)

      if(nmax.gt.nx) stop "d2naive: make nx larger."
      nlast=min(nmax-(mmax-1)*id-1,ntmax)
      do 10 ndt=ntmin,nlast
         do 20 n=ndt+1,nmax
 20         d(n)=max(abs(x(n)-x(n-ndt)),tiny)
         do 10 n=ndt+1+(mmax-1)*id,nmax
            dmax=d(n)
            do 30 m=2,mmin-1
 30            dmax=max(dmax,d(n-(m-1)*id))
            j=int(a*(log(dmax)-scl))
            do 10 m=mmin,mmax
               if(d(n-(m-1)*id).gt.dmax) then
                  dmax=d(n-(m-1)*id)
                  j=int(a*(log(dmax)-scl))
               endif
 10            c(j,m)=c(j,m)+1
      end
