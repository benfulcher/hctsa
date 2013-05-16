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
c   autocorrelation function with periodic continuation
c   author T. Schreiber (1999)
c
c-------------------------------------------------------------------
c get cost function specific options
c
      subroutine opts_cost(ncol)
      parameter(mlag=100000)
      dimension c0(mlag), c(mlag)
      common /costcom/ nlag, c0, c, sd, sc, iweight
      
      nlag=imust('D')
      iweight=ican('W',0)
      ncol=1
      end

c-------------------------------------------------------------------
c print version information on cost function
c
      subroutine what_cost()
      call ptext("Cost function: periodic autocorrelation")
      end

c-------------------------------------------------------------------
c print cost function specific usage message
c
      subroutine usage_cost()
      call ptext("Cost function options: -D# [-W#]")
      call popt("D","number of lags")
      call popt("W",
     .   "average: 0=max(c) 1=|c|/lag 2=(c/lag)**2 3=max(c)/lag (0)")
      end

c-------------------------------------------------------------------
c initialise all that is needed for cost function
c
      subroutine cost_init()
      parameter(mlag=100000)
      dimension c0(mlag), c(mlag)
      common /costcom/ nlag, c0, c, sd, sc, iweight

      if(nlag.gt.mlag) write(istderr(),'(a)') 
     .   "truncated to ", mlag," lags"
      nlag=min(mlag,nlag)
      call auto(nlag,c0)
      end

c-------------------------------------------------------------------
c initial transformation on time series and its inverse
c
      subroutine cost_transform(nmax,mcmax,nxdum,x)
      dimension x(nmax)
      parameter(mlag=100000)
      dimension c0(mlag), c(mlag)
      common /costcom/ nlag, c0, c, sd, sc, iweight

      call normal1(nmax,x,sc,sd)
      end

      subroutine cost_inverse(nmax,mcmax,nxdum,x,y)
      dimension x(nmax), y(nmax)
      parameter(mlag=100000)
      dimension c0(mlag), c(mlag)
      common /costcom/ nlag, c0, c, sd, sc, iweight
      
      do 10 n=1,nmax
 10      y(n)=x(n)*sd+sc
      end

c-------------------------------------------------------------------
c compute full cost function from scratch
c
      subroutine cost_full(iv)
      parameter(mlag=100000)
      dimension c0(mlag), c(mlag)
      common /costcom/ nlag, c0, c, sd, sc, iweight
      common nmax,cost

      call auto(nlag,c)
      cc=0
      do 10 n=1,nlag
 10      call aver(cc,c0(n)-c(n),n)
      cost=cc
      end

c-------------------------------------------------------------------
c compute changed cost function on exchange of n1 and n2 
c
      subroutine cost_update(n1,n2,cmax,iaccept,iv)
      parameter(mlag=100000,nx=100000)
      dimension c0(mlag), c(mlag), ccop(mlag), x(nx)
      common /costcom/ nlag, c0, c, sd, sc, iweight
      common nmax,cost,temp,cmin,rate,x

      comp=0
      iaccept=0
      do 10 n=1,nlag
         cc=c(n)
         dx=x(n2)-x(n1)
         nd1=n1-n
         if(nd1.lt.1) nd1=nd1+nmax
         if(nd1.ne.n2) cc=cc+dx*x(nd1)
         nu1=n1+n
         if(nu1.gt.nmax) nu1=nu1-nmax
         if(nu1.ne.n2) cc=cc+dx*x(nu1)
         nd2=n2-n
         if(nd2.lt.1) nd2=nd2+nmax
         if(nd2.ne.n1) cc=cc-dx*x(nd2)
         nu2=n2+n
         if(nu2.gt.nmax) nu2=nu2-nmax
         if(nu2.ne.n1) cc=cc-dx*x(nu2)
         call aver(comp,c0(n)-cc,n)
         if(comp.ge.cmax) return
 10      ccop(n)=cc
      cost=comp  ! if got here: accept
      iaccept=1
      call exch(n1,n2)
      do 20 n=1,nlag
 20      c(n)=ccop(n)
      end

c-------------------------------------------------------------------
c compute autocorrelation from scratch
c
      subroutine auto(nlag,c)
      parameter(nx=100000)
      dimension c(*), x(nx)
      common nmax,cost,temp,cmin,rate,x
      
      do 10 n=1,nlag
         cc=0
         do 20 i=1,nmax
            ii=i-n
            if(ii.lt.1) ii=ii+nmax
 20         cc=cc+x(ii)*x(i)
 10      c(n)=cc
      end

c-------------------------------------------------------------------
c weighted average of autocorrelation 
c
      subroutine aver(cav,dc,n)
      parameter(mlag=100000)
      dimension c0(mlag), c(mlag)
      common /costcom/ nlag, c0, c, sd, sc, iweight
      common nmax

      if(iweight.eq.0) then
         cav=max(cav,abs(dc)/real(nmax))
      else if(iweight.eq.1) then
         cav=cav+abs(dc)/real(nmax*n)
      else if(iweight.eq.2) then
         cav=cav+(dc/real(nmax*n))**2
      else
         cav=max(cav,abs(dc)/real(nmax*n))
      endif
      end

