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
c   part of the TISEAN randomize package for constraint surrogates
c   cost function
c   binned spike train autocorrelation function
c   author T. Schreiber (1999)
c
c-------------------------------------------------------------------
c get cost function specific options
c
      subroutine opts_cost(ncol)
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight

      iweight=ican('W',0)
      bininv=1./fmust("d")
      totbin=fmust("D")
      nbin=min(int(totbin*bininv)+1,nhist)
      inter=lopt("i",1)
      ncol=1
      end

c-------------------------------------------------------------------
c print version information on cost function
c
      subroutine what_cost()
      call ptext("Cost function: spike train autocorrelation function")
      end

c-------------------------------------------------------------------
c print cost function specific usage message
c
      subroutine usage_cost()
      call ptext("Cost function options: -d# -D# [-i -W#]")
      call popt("d","time span of one bin")
      call popt("D","total time spanned")
      call popt("i","expect intervals rather than times")
      call popt("W",
     .   "average: 0=max(c) 1=|c|/lag 2=(c/lag)**2 (0)")
      end

c-------------------------------------------------------------------
c initialise all that is needed for cost function
c
      subroutine cost_init()
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight

      call sauto(nbin,bininv,ihist0)
      end

c-------------------------------------------------------------------
c initial transformation on time series and its inverse
c here: series internally stored as intervals 
c
      subroutine cost_transform(nmax,mcmax,nxdum,x)
      parameter(nx=100000)
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight
      dimension nxclu(nx)
      common /permutecom/ mxclu, nxclu
      dimension x(*), lx(nx)

      if(inter.eq.1) return
      call sort(nmax,x,lx)
      do 10 n=nmax,2,-1
 10      x(n)=x(n)-x(n-1)
      mxclu=mxclu+1
      nxclu(mxclu)=1
      end

      subroutine cost_inverse(nmax,mcmax,nxdum,x,y)
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight
      dimension x(*), y(*)
      
      do 10 n=1,nmax
 10      y(n)=x(n)
      if(inter.eq.1) return
      do 20 n=2,nmax
 20      y(n)=y(n)+y(n-1)
      end

c-------------------------------------------------------------------
c compute full cost function from scratch
c
      subroutine cost_full(iv)
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight
      common nmax,cost

      call sauto(nbin,bininv,ihist)
      cost=aver(ihist0,ihist)
      if(iv.ne.0) call dump()
      end

c-------------------------------------------------------------------
c compute changed cost function on exchange of n1 and n2 
c
      subroutine cost_update(nn1,nn2,cmax,iaccept,iv)
      parameter(nx=100000)
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight
      dimension ihcop(nhist), x(nx)
      common nmax,cost,temp,cmin,rate,x

      n1=min(nn1,nn2)
      n2=max(nn1,nn2)
      comp=0
      iaccept=0
      do 10 i=1,nbin
 10      ihcop(i)=ihist(i)
      dx=0
      do 20 nn=n1,1,-1
         if(nn.lt.n1) dx=dx+x(nn)
         if(int(dx*bininv)+1.gt.nbin) goto 1
         dxx=dx
         do 30 nnn=n1,n2-1
            dxx=dxx+x(nnn)
            il=int(dxx*bininv)+1
            if(il.gt.nbin) goto 20
 30         ihcop(il)=ihcop(il)-1
 20      continue
 1    dx=0
      do 40 nn=n2,1,-1
         if(nn.lt.n2) dx=dx+x(nn)
         if(int(dx*bininv)+1.gt.nbin) goto 2
         dxx=dx
         do 50 nnn=n2,nmax
            dxx=dxx+x(nnn)
            il=int(dxx*bininv)+1
            if(il.gt.nbin) goto 40
 50         ihcop(il)=ihcop(il)-1
 40      continue
 2    call exch(n1,n2)
      dx=0
      do 60 nn=n1,1,-1
         if(nn.lt.n1) dx=dx+x(nn)
         if(int(dx*bininv)+1.gt.nbin) goto 3
         dxx=dx
         do 70 nnn=n1,n2-1
            dxx=dxx+x(nnn)
            il=int(dxx*bininv)+1
            if(il.gt.nbin) goto 60
 70         ihcop(il)=ihcop(il)+1
 60      continue
 3    dx=0
      do 80 nn=n2,1,-1
         if(nn.lt.n2) dx=dx+x(nn)
         if(int(dx*bininv)+1.gt.nbin) goto 4
         dxx=dx
         do 90 nnn=n2,nmax
            dxx=dxx+x(nnn)
            il=int(dxx*bininv)+1
            if(il.gt.nbin) goto 80
 90         ihcop(il)=ihcop(il)+1
 80      continue
 4    comp=aver(ihist0,ihcop)
      if(comp.ge.cmax) then
         call exch(n1,n2)
         return
      endif
      cost=comp  ! if got here: accept
      iaccept=1
      if(iv.ne.0) call panic(ihcop)
      do 100 i=1,nbin
 100     ihist(i)=ihcop(i)
      end

c-------------------------------------------------------------------
c compute autocorrealtion from scratch
c
      subroutine sauto(nbin,bininv,ihist)
      parameter(nx=100000)
      dimension ihist(*)
      common nmax,cost,temp,cmin,rate,x
      dimension x(nx)

      do 10 i=1,nbin
 10      ihist(i)=0
      do 20 n1=1,nmax
         dx=0
         do 30 n2=n1,nmax
            dx=dx+x(n2)
            il=int(dx*bininv)+1
            if(il.gt.nbin) goto 20
 30         ihist(il)=ihist(il)+1
 20      continue
      end

c-------------------------------------------------------------------
c weighted average of autocorrelation 
c
      function aver(ih1,ih2)
      parameter(nhist=100000)
      dimension ih1(nhist), ih2(nhist)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight

      aver=0
      if(iweight.eq.0) then
         do 10 i=1,nbin
 10         aver=max(aver,real(abs(ih1(i)-ih2(i))))
      else if(iweight.eq.1) then
         do 20 i=1,nbin
 20         aver=aver+real(abs(ih1(i)-ih2(i)))/real(i)
      else if(iweight.eq.2) then
         do 30 i=1,nbin
 30         aver=aver+(ih1(i)-ih2(i))**2/real(i)
      endif
      end

c-------------------------------------------------------------------
c diagnostic output
c
      subroutine dump()
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight

      write(istderr(),'(5hgoal ,4i12)') (ihist0(n),n=1,min(4,nbin))
      write(istderr(),'(5his   ,4i12)') (ihist(n),n=1,min(4,nbin))
      write(istderr(),'(5hmiss ,4i12)') 
     .   (abs(ihist0(n)-ihist(n)),n=1,min(4,nbin))
      write(istderr(),'()')
      end

      subroutine panic(ihcop)
      parameter(nhist=100000)
      dimension ihist0(nhist), ihist(nhist)
      common /costcom/ inter, bininv, nbin, ihist0, ihist, iweight
      dimension ihcop(*)

      call cost_full(0)
      write(istderr(),'(7hupdate ,4i12)') (ihcop(n),n=1,min(4,nbin))
      write(istderr(),'(7hfresh  ,4i12)') (ihist(n),n=1,min(4,nbin))
      write(istderr(),'(7hdiscr  ,4i12)') 
     .   (abs(ihcop(n)-ihist(n)),n=1,min(4,nbin))
      write(istderr(),'()')
      end
