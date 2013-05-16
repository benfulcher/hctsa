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
c   spike train power spectrum
c   author T. Schreiber (1999)
c
c-------------------------------------------------------------------
c get cost function specific options
c
      subroutine opts_cost(ncol)
      parameter(mfreq=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight

      iweight=ican('W',0)
      fmax=fcan("F",0)
      nfreq=ican("#",0)
      inter=lopt("i",1)
      ncol=1
      end

c-------------------------------------------------------------------
c print version information on cost function
c
      subroutine what_cost()
      call ptext("Cost function: spike train power spectrum")
      end

c-------------------------------------------------------------------
c print cost function specific usage message
c
      subroutine usage_cost()
      call ptext("Cost function options: [-F# -## -w# -i]")
      call popt("W",
     .   "average: 0=max(s) 1=|s|/f 2=(s/f)**2 3=|s| (0)")
      call popt("F","maximal frequency (2*l / total time)")
      call popt("#","number of frequencies (F* total time /2)")
      call popt("w","frequency resolution (0)")
      call popt("i","expect intervals rather than times")
      end

c-------------------------------------------------------------------
c initialise all that is needed for cost function
c
      subroutine cost_init()
      parameter(nx=100000,mfreq=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      dimension x(nx)
      common nmax,cost,temp,cmin,rate,x

      if(fmax.le.0.) fmax=2*nmax/(x(nmax)-x(1))
      if(nfreq.le.0) nfreq=fmax*(x(nmax)-x(1))/2
      if(nfreq.gt.mfreq) write(istderr(),'(a)') 
     .   "truncated to ", mfreq," frequencies"
      nfreq=min(mfreq,nfreq)
      write(istderr(),*) "randomize_spikespec: total time covered: ", 
     .   x(nmax)-x(1)
      write(istderr(),*) "randomize_spikespec: computing ", nfreq, 
     .   " frequencies up to ", fmax
      call sspect(nfreq,fmax/nfreq,sp0r,sp0i,sp0)
      end

c-------------------------------------------------------------------
c initial transformation on time series and its inverse
c
      subroutine cost_transform(nmax,mcmax,nxdum,x)
      parameter(mfreq=100000,nx=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      dimension x(nx), lx(nx)

      if(inter.eq.0) goto 1
      do 10 n=2,nmax
 10      x(n)=x(n)+x(n-1)
 1    call sort(nmax,x,lx)
      end

      subroutine cost_inverse(nmax,mcmax,nxdum,x,y)
      parameter(mfreq=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      dimension x(nmax), y(nmax)
      
      do 10 n=1,nmax
 10      y(n)=x(n)
      if(inter.ne.1) return
      do 20 n=nmax,2,-1
 20      y(n)=y(n)-y(n-1)
      end

c-------------------------------------------------------------------
c compute full cost function from scratch
c
      subroutine cost_full(iv)
      parameter(mfreq=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      common nmax,cost

      call sspect(nfreq,fmax/nfreq,spr,spi,sp)

      cc=0
      do 10 i=1,nfreq
 10      call aver(cc,sp(i)-sp0(i),i)
      cost=cc
      if(iv.ne.0) call dump()
      end

c-------------------------------------------------------------------
c compute changed cost function on exchange of n1 and n2 
c
      subroutine cost_update(n1,n2,cmax,iaccept,iv)
      parameter(mfreq=100000,nx=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      dimension sprcop(mfreq), spicop(mfreq), spcop(mfreq), x(nx)
      common nmax,cost,temp,cmin,rate,x
      data pi/3.1415926/

      comp=0
      iaccept=0
      do 10 i=1,nfreq
         f=i*(fmax/nfreq)
         omega=2*pi*f
         xx=x(n1-1)+x(n1+1)-x(n1)
         sprcop(i)=spr(i)-cos(omega*x(n1))+cos(omega*xx)
         spicop(i)=spi(i)-sin(omega*x(n1))+sin(omega*xx)
         spcop(i)=sprcop(i)**2+spicop(i)**2
         call aver(comp,sp0(i)-spcop(i),i)
 10      if(comp.ge.cmax) return
      cost=comp  ! if got here: accept
      iaccept=1
      call exch(n1,n2)
      if(iv.ne.0) call panic(spcop)
      do 20 i=1,nfreq
         spr(i)=sprcop(i)
         spi(i)=spicop(i)
 20      sp(i)=spcop(i)
      end

c-------------------------------------------------------------------
c compute spectrum from scratch
c
      subroutine sspect(nfreq,fres,spr,spi,sp)
      parameter(nx=100000)
      dimension spr(*), spi(*), sp(*), x(nx)
      common nmax,cost,temp,cmin,rate,x
      data pi/3.1415926/

      do 10 i=1,nfreq
         f=i*fres
         omega=2*pi*f
         sr=0
         si=0
         do 20 n=1,nmax
            sr=sr+cos(omega*x(n))
 20         si=si+sin(omega*x(n))
         spr(i)=sr
         spi(i)=si
 10      sp(i)=sr**2+si**2
      end

c-------------------------------------------------------------------
c weighted average of autocorrelation 
c
      subroutine aver(cav,dc,n)
      parameter(mfreq=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight

      if(iweight.eq.0) then
         cav=max(cav,abs(dc))    ! max (L-infinity) norm
      else if(iweight.eq.1) then
         cav=cav+abs(dc)/n       ! L-1 norm (sum of moduli), weighted by freq.
      else if(iweight.eq.2) then
         cav=cav+(dc/n)**2       ! L-2 norm (sum of squares), weighted by freq.
      else if(iweight.eq.2) then
         cav=cav+abs(dc)         ! L-1 norm (sum of moduli)
      endif
      end

c-------------------------------------------------------------------
c diagnostic output
c
      subroutine dump()
      parameter(mfreq=100000)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      common nmax

      write(istderr(),'(5hgoal ,4f9.3)') (sp0(n),n=1,min(4,nfreq))
      write(istderr(),'(5his   ,4f9.3)') (sp(n),n=1,min(4,nfreq))
      write(istderr(),'(5hmiss ,4f9.3)') 
     .   ((sp0(n)-sp(n)),n=1,min(4,nfreq))
      write(istderr(),'()')
      end

      subroutine panic(spcop)
      parameter(mfreq=100000)
      dimension spcop(*)
      dimension sp0r(mfreq), sp0i(mfreq), spr(mfreq), spi(mfreq),
     .   sp0(mfreq), sp(mfreq)
      common /costcom/ nfreq, fmax, inter,
     .    sp0r, sp0i, sp0, spr, spi, sp, sd, sc, iweight
      common nmax

      call cost_full(0)
      write(istderr(),'(7hupdate ,4f9.3)') 
     .   (spcop(n),n=1,min(4,nfreq))
      write(istderr(),'(7hfresh  ,4f9.3)') (sp(n),n=1,min(4,nfreq))
      write(istderr(),'(7hdiscr  ,4f9.3)') 
     .   ((spcop(n)-sp(n)),n=1,min(4,nfreq))
      write(istderr(),'()')
      end
