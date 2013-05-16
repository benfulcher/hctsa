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
c   space time separation plot
c   see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997,2004)
c   author T. Schreiber (1998) based on earlier version
c==========================================================================
      parameter(nx=1000000,mdt=500,mfrac=100)
      dimension x(nx), stp(mfrac,mdt)
      character*72 file, fout
      data idt/1/, perc/0.05/, ndt/100/
      data iverb/1/

      call whatido("space-time separation plot",iverb)
      id=imust("d")
      m=imust("m")
      idt=ican("#",idt)
      ndt=min(ican("t",ndt),mdt)
      perc=fcan("%",perc)
      nfrac=min(mfrac,int(1/perc))
      perc=1./real(nfrac)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)
      if(iv_io(iverb).eq.1) write(istderr(),*) "computing ", nfrac, 
     .   " levels at fractions ", perc, 2*perc, "..."

      call nthstring(1,file)
      call readfile(nmax,x,nexcl,jcol,file,iverb)
      call minmax(nmax,x,xmin,xmax)
      call stplot(nmax,x,id,m,xmax-xmin,stp,nfrac,ndt,idt)
      if(isout.eq.1) call addsuff(fout,file,"_stp")
      call outfile(fout,iunit,iverb)
      do 10 iper=1,mfrac
         do 20 it=1,ndt
 20         write(iunit,*) it*idt, stp(iper,it)
 10      write(iunit,'()')
      end

      subroutine usage()
c usage message
      
      call whatineed(
     .   " -d# -m# [-## -t# -%# -o outfile -l# -x# -c# -V# -h] file")
      call popt("d","delay")
      call popt("m","embedding dimension")
      call popt("#","time resolution (1)")
      call popt("t","time steps (100, <500)")
      call popt("%","fraction at wich to create levels (0.05, >0.01)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_stp")
      call pall()
      stop
      end

      subroutine stplot(nmax,y,id,m,epsmax,stp,nfrac,mdt,idt)
      parameter(meps=1000,mfrac=100)
      dimension y(nmax),stp(mfrac,mdt),ihist(meps)

      do 10 it=1,mdt
         do 20 ieps=1,meps
 20         ihist(ieps)=0
         do 30 n=it*idt+(m-1)*id+1,nmax
            dis=0                            ! compute distance in m dimensions
            do 40 me=0,m-1
 40            dis=max(dis,abs(y(n-me*id)-y(n-me*id-it*idt)))
            ih=min(int(meps*dis/epsmax)+1,meps)
 30         ihist(ih)=ihist(ih)+1
         do 10 ifrac=1,nfrac
            need=(nmax-it*idt-(m-1)*id)*ifrac/real(nfrac)
            is=0
            do 50 ieps=1,meps
               is=is+ihist(ieps)
 50            if(is.ge.need) goto 1
 1          stp(ifrac,it)=ieps*epsmax/meps
 10      continue
      end





