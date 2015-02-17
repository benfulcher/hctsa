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
c cross-recurrence plot
c authors H. Kantz & T. Schreiber (2004)
c modified by H. Kantz Feb 2007 (multivariate version)
      program xrecur
      parameter(nx=1000000,me=30,meps=1000,mx=10,nx1=10000000)
      dimension x(nx,mx), y(nx,mx)
      dimension x1(nx1),y1(nx1)
      integer mlist(2), inot(nx1), icol(mx)
      character*72 file1, file2, fout
      data eps/1e-3/, id0/1/
      data iverb/1/, xperc/100./

      call whatido(
     ."cross recurrence plot of two scalar or vector valued data sets"
     .             ,iverb)

c ================================================================
c assume two input files (file1, file2) 
c - with identical structure in terms of colums
c - with identical embedding spaces 
c - with individual length and exclusions: l1,l2, x1,x2
c - univariate data: maximum time series length nx1
c - multivariate data: maximum length nx, maximum dimension mx
c
c either: fixed epsilon, plot certain percentage of all 
c         neighbours found
c or:     fix the numbers of points of the y time series to be 
c         found as neighbours of the x time series: non-symmetric! 
c norm: max-norm after rescaling of the individual components
c       unless rescaling is switched off
c
c=================================================================

      id=ican("d",id0)

      mmax=2
      mdim=1

      call imcan("m",2,mc,mlist)
      if (mc.ge.2) then
       mmax=mlist(2)
       mdim=mlist(1)
       if (mc.gt.2) print*,'extra arguments of -m ignored'
      endif
      ntmin=0
      kmin=ican("k",0)
      idown1=ican("s",1)
      idown2=ican("S",1)
      eps=fcan("r",eps)
      xperc=fcan("%",xperc)
      nmaxx=ican("l",nx)
      nmaxy=ican("L",nx)
      nexcl1=ican("x",0)
      nexcl2=ican("X",0)
      iscal=0
      iscal=lopt("n",1)

      call columns(mc,mx,icol)
      if (mc.gt.0.and.mc.ne.mdim) stop 'improper number of columns'
      isout=igetout(fout,0)
      if(fout.eq." ") isout=1
c     
      call nthstring(1,file1)
      if(file1.eq."-") stop 'missing input filename'
      if (mdim.gt.1) 
     .  call xreadfile(nmaxx,mdim,nx,x,nexcl1,icol,file1,iverb)
      if (mdim.eq.1) 
     .  call xreadfile(nmaxx,1,nx1,x1,nexcl1,icol,file1,iverb)
c     
      call nthstring(2,file2)
      if(file2.eq."-") stop 'missing second input filename'
      if (mdim.gt.1) 
     . call xreadfile(nmaxy,mdim,nx,y,nexcl2,icol,file2,iverb)
      if (mdim.eq.1) 
     . call xreadfile(nmaxy,1,nx1,y1,nexcl2,icol,file2,iverb)

      if(isout.eq.1) then 
        call addsuff(fout,file1,"_")
        call addsuff(fout,fout,file2)
        call addsuff(fout,fout,"_xrec")
      endif

c     rescale data if flag -n is not set (iscal=1)
      if (iscal.ne.1) then 
       print*,'normalizing data to unit interval'

       if (mdim.gt.1) then

c      rescaling each component of file 1
        do imx=1,mdim
         xmin=x(1,imx)
         xmax=x(1,imx)
         do i1=2,nmaxx
          xmax=max(xmax,x(i1,imx))
          xmin=min(xmin,x(i1,imx))
         enddo
         scal=.9999d0/(xmax-xmin)
         do i1=1,nmaxx
          x(i1,imx)=(x(i1,imx)-xmin)*scal
         enddo
        enddo

c     rescaling each component of file 2

        do imx=1,mdim
         xmin=y(1,imx)
         xmax=y(1,imx)
         do i1=2,nmaxy
          xmax=max(xmax,y(i1,imx))
          xmin=min(xmin,y(i1,imx))
         enddo
         scal=.9999d0/(xmax-xmin)
         do i1=1,nmaxy
          y(i1,imx)=(y(i1,imx)-xmin)*scal
         enddo
        enddo

       else

c     rescaling the single component of file 1
         xmin=x1(1)
         xmax=x1(1)
         do i1=2,nmaxx
          xmax=max(xmax,x1(i1))
          xmin=min(xmin,x1(i1))
         enddo
         scal=.9999d0/(xmax-xmin)
         do i1=1,nmaxx
          x1(i1)=(x1(i1)-xmin)*scal
         enddo

c     rescaling the single component of file 2

         xmin=y1(1)
         xmax=y1(1)
         do i1=2,nmaxy
          xmax=max(xmax,y1(i1))
          xmin=min(xmin,y1(i1))
         enddo
         scal=.9999d0/(xmax-xmin)
         do i1=1,nmaxy
          y1(i1)=(y1(i1)-xmin)*scal
         enddo

       endif
      endif

      call outfile(fout,iunit,1)
      ntot=0

      if (kmin.eq.0) then
c     search all neighbours with distance < eps

        xperc=xperc/100.

        if (mdim.eq.1) then
          call crossrec(nmaxx,x1,nmaxy,y1,eps,
     .                   id,mmax,iunit,xperc,ntot,idown1,idown2)
         else
          call mcrossrec(nmaxx,x,nmaxy,y,eps,
     .              id,mmax,mdim,iunit,xperc,ntot,idown1,idown2)
        endif

c>-----------------------------------------------------------
      else 

       do i=(mmax-1)*id+1,nmaxx
        inot(i)=1
       enddo
       epsfac=1.1
       eps=eps/epsfac

       do 10 io=1,100
         eps=eps*epsfac
         if (mdim.eq.1) then
          call crossrec1(nmaxx,x1,nmaxy,y1,eps,
     .               id,mmax,kmin,iunit,inot,ntodo,ntot,idown1,idown2)
         else
          call mcrossrec1(nmaxx,x,nmaxy,y,eps,id,mmax,mdim,kmin,iunit,
     .                    inot,ntodo,ntot,idown1,idown2)
         endif
         if (iverb.eq.1) print*,eps,ntodo,ntot
         if (ntodo.eq.0) goto 7
 10    continue
c>--------------------------------------------------------
      endif

 7    close(iunit)
      print*,ntot,' points contained in the recurrence plot'
      if (kmin.gt.0) print*,'last epsilon:',eps
      if (kmin.gt.0) print*,'average number of neighbours:',
     .                      real(ntot)*idown1/nmaxx
      if (kmin.eq.0) print*,'using eps=',eps
      end
c========================================================
      subroutine usage()
c usage message

      call whatineed(
     ."[ -m#,# -d# -r# -k# -o outfile -l# -x# -L# -X# -c#[,#] -%# -V# 
     . -n -h] file1 file2")
      call popt("m","# of components, embedding dimension [1,2]")
      call popt("c","columns to be read [1,2,3,...,# of components]")
      call popt("d","delay [1]")
      call popt("r",
     ."diameter of the neighbourhood as absolute value [.001]")
      call popt("k",
     ."find the # closest points, starting with diameter r")
      call popt("%",
     ."print only percentage of dots [100], no effect if -k is set")
      call popt("l","length of time series 1 to be read [all data]")
      call popt("x","# of initial lines in 1 to be skipped [0]")
      call popt("s","use only every # delay vector of file 1 [1]")
      call popt("L","length of time series 2 to be read [all data]")
      call popt("X","# of initial lines in 2 to be skipped [0]")
      call popt("S","use only every # delay vector of file 2 [1]")
      call popt("n","if set: do NOT normalize data to unit interval")
      call pout("file1_file2_xrec")
      call pall()
      stop
      end
c>--------------------------------------------------------------------
      subroutine crossrec(nmaxx,x,nmaxy,y,eps,
     .                    id,m,iunit,xperc,ntot,idown1,idown2)
      parameter(im=100,ii=100000000,nx=1000000)
      dimension y(nx),x(nx)
      dimension jh(0:im*im),jpntr(nx),nlist(nx)
      nseed=13413241

      if(nmaxx.gt.nx.or.nmaxy.gt.nx) stop "crossrec: make nx larger."
      mb=min(m,2)

      call base(nmaxy,y,id,mb,jh,jpntr,eps)
      nnull=(m-1)*id+1
      do 20 n=nnull,nmaxx,idown1
         call neigh(nx,x,y,n,nx,id,mb,jh,jpntr,eps,nlist,nfound)
         do 30 nn=1,nfound                   ! all neighbours in two dimensions
            np=nlist(nn)
            if (np.lt.nnull) goto 30
            if (mod(np-nnull,idown2).ne.0) goto 30
            do 40 i=mb,m-1
               if(abs(x(n-i*id)-y(np-i*id)).ge.eps) goto 30
 40         continue
            call random(nseed,rr)
            if (rr.le.xperc) write(iunit,*)n,np
            if (rr.le.xperc) ntot=ntot+1
 30         continue
 20      continue
      end
c>--------------------------------------------------------------------
      subroutine mcrossrec(nmaxx,x,nmaxy,y,eps,
     .                id,m,mdim,iunit,xperc,ntot,idown1,idown2)
      parameter(im=100,nx=1000000,mx=10)
      dimension y(nx,mx),x(nx,mx)
      dimension jh(0:im*im),jpntr(nx),nlist(nx)
      dimension vx(mx)
      nseed=134512331

      if(nmaxx.gt.nx.or.nmaxy.gt.nx) stop "mcrossrec: make nx larger."

      call mbase(nmaxy,mdim,nx,y,id,1,jh,jpntr,eps)
      nnull=(m-1)*id+1
      do 20 n=nnull,nmaxx,idown1
         do ii=1,mdim
          vx(ii)=x(n,ii)
         enddo
         call mneigh2(nmaxy,mdim,y,nx,vx,jh,jpntr,eps,nlist,nfound)
         do 30 nn=1,nfound               ! all neighbours in mdim dimensions
            np=nlist(nn)
            if (np.lt.nnull) goto 30
            if (mod(np-nnull,idown2).ne.0) goto 30
            do 40 i=1,m-1
             do 41 iim=1,mdim
               if(abs(x(n-i*id,iim)-y(np-i*id,iim)).ge.eps) goto 30
 41          continue
 40         continue
             call random(nseed,rr)
             if (rr.le.xperc) write(iunit,*)n,np
             if (rr.le.xperc) ntot=ntot+1
 30      continue
 20   continue
      return
      end

      subroutine crossrec1(nmaxx,x,nmaxy,y,eps,id,m,kmin,iunit,
     .                     inot,ntodo,ntot,idown1,idown2)
      parameter(im=100,ii=100000000,nx=1000000)
      dimension y(nmaxy),x(nmaxx),inot(nmaxx)
      dimension jh(0:im*im),jpntr(nx),nlist(nx)
      ntodo=0

      if(nmaxx.gt.nx.or.nmaxy.gt.nx) stop "crossrec1: make nx larger."
      mb=min(m,2)

      call base(nmaxy,y,id,mb,jh,jpntr,eps)
      nnull=(m-1)*id+1
      do 20 n=nnull,nmaxx,idown1
         if (inot(n).eq.0) goto 20
         ntodo=ntodo+1
         call neigh(nx,x,y,n,nx,id,mb,jh,jpntr,eps,nlist,nfound)
         if (nfound.lt.kmin) goto 20
         nreal=0
         do 30 nn=1,nfound                   ! all neighbours in two dimensions
            np=nlist(nn)
            if(np.lt.nnull) goto 30
            if (mod(np-nnull,idown2).ne.0) goto 30
            do 40 i=mb,m-1
               if(abs(x(n-i*id)-y(np-i*id)).ge.eps) goto 30
 40         continue
            nreal=nreal+1
            nlist(nreal)=np
 30      continue
         if (nreal.lt.kmin) goto 20
           ntodo=ntodo-1
           inot(n)=0
           ntot=ntot+nreal
           do in=1,nreal
            write(iunit,*)n,nlist(in)
           enddo
 20      continue
      end
c>--------------------------------------------------------------------
      subroutine mcrossrec1(nmaxx,x,nmaxy,y,eps,
     .                      id,m,mdim,kmin,iunit,inot,
     .                      ntodo,ntot,idown1,idown2)
      parameter(im=100,nx=1000000,mx=10)
      dimension y(nx,mx),x(nx,mx),inot(nmaxx)
      dimension jh(0:im*im),jpntr(nx),nlist(nx)
      dimension vx(mx)
      ntodo=0

      if(nmaxx.gt.nx.or.nmaxy.gt.nx) stop "mcrossrec1: make nx larger."

      call mbase(nmaxy,mdim,nx,y,id,1,jh,jpntr,eps)
      nnull=(m-1)*id+1
      do 20 n=nnull,nmaxx,idown1
         if (inot(n).eq.0) goto 20
         ntodo=ntodo+1
         do ii=1,mdim
          vx(ii)=x(n,ii)
         enddo
         call mneigh2(nmaxy,mdim,y,nx,vx,jh,jpntr,eps,nlist,nfound)
         if (nfound.lt.kmin) goto 20
         nreal=0
         do 30 nn=1,nfound               ! all neighbours in mdim dimensions
            np=nlist(nn)
            if(np.lt.nnull) goto 30
            if (mod(np-nnull,idown2).ne.0) goto 30
            do 40 i=0,m-1
             do 41 iim=1,mdim
               if(abs(x(n-i*id,iim)-y(np-i*id,iim)).ge.eps) goto 30
 41          continue
 40         continue
            nreal=nreal+1
            nlist(nreal)=np 
 30      continue
         if (nreal.lt.kmin) goto 20
           inot(n)=0
           ntot=ntot+nreal
           do in=1,nreal
            write(iunit,*)n,nlist(in)
           enddo
           ntodo=ntodo-1
 20      continue
      return
      end

      subroutine random(iseed,s)
c
c     random number generator of Park & Miller
      integer*8 ifac,ibase,iargument
      ifac=7**5
      ibase=2**30-1
      im=im+2**30
      iargument=iseed
      iargument=mod(iargument*ifac,ibase)
      s=float(iargument)/float(ibase)
      iseed=iargument
      return
      end
c>------------------------------------
