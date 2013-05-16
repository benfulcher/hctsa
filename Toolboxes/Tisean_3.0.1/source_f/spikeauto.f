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
c   spike train autocorrelation function
c   author T. Schreiber (1998) based on earlier versions
c===========================================================================
      parameter(nx=1000000, nhist=100000)
      dimension x(nx), lx(nx), ihist(nhist)
      character*72 file, fout
      data iverb/1/

      call whatido("spike train autocorrelation function",iverb)
      bin=fmust("d")
      totbin=fmust("D")
      nbin=int(totbin/bin)+1
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
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
         do 30 i=1,nbin
 30         ihist(i)=0
         do 40 n1=1,nmax
            do 50 n2=n1+1,nmax
               il=int((x(n2)-x(n1))/bin)+1
               if(il.gt.nbin) goto 40
 50            ihist(il)=ihist(il)+1
 40         continue
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_sco")
         call outfile(fout,iunit,iverb)
         do 60 i=1,nbin
 60         write(iunit,*) (i-0.5)*bin, ihist(i)
 10      if(iunit.ne.istdout()) close(iunit)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# -D# [-i -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("d","time span of one bin")
      call popt("D","total time spanned")
      call popt("i","expect intervals rather than times")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_sco")
      call pall()
      stop
      end



