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
c   lorenz.f
c   integrates the Lorenz system with Runge Kutta fourth order
c   author: H. Kantz (2007) based on earlier versions
c   with optional noise
c===========================================================================
c

      real*8 x(3),u(3,3),sliap(3),bb,ss,rr,r1,r2,dh,s
      character*72 fout
      data iverb/1/

      iverb=ican('V',iverb)
      call whatido("integration of the Lorenz system",iverb)
      irun=imust('l')
      itrans=ican('x',100)
      rr=fcan('R',28.0)
      ss=fcan('S',10.0)
      bb=fcan('B',2.666666667)
      isamp=ican('f',100)
      sn=fcan('r',0.)
c      ilyap=lopt('L',1)

      isout=igetout(fout,iverb)

      if(isout.eq.1) fout="lorenz.dat"
      call outfile(fout,iunit,iverb)

cc intermittency parameters
c       ss=10.d0
c       rr=166.11d0
c       bb=8.d0/3.d0

      iseed1=6456423
      iseed2=7243431

      xav=0.
      xsq=0.
      rsq=0.

c     step width of Runge Kutta integration dh:
      dh=.0005d0
c     time intervals between re-orthogonalization of tangent space
c            vectors: 0.01 units of time.
      ireno=.01d0/dh
c     length of transient in iteration steps:
      itrans=real(itrans)/dh
      totaltime=real(irun)/real(isamp)
      istep=1.d0/dh/isamp

      if (iverb.eq.1) 
     . write(istderr(),*)'Lorenz trajectory covering',totaltime,
     .                  ' time units'

c      x(1)=sqrt(s*(r+1.d0))+2.
c      x(2)=x(1)-1.d0
c      x(3)=r

      x(1)=5.
      x(2)=-10.
      x(3)=3.

      do 1 i=1,3
       sliap(i)=0.d0
       do j=1,3
        u(i,j)=0.d0
       enddo
       u(i,i)=1.d0
1     continue

      do 10 i2=1,itrans

        call RUKU(3,x,u,dh,bb,ss,rr)

        if (mod(i2,ireno).eq.0) then
          call norm(u,1,s)
          do i=2,3
            do j=1,i-1
              call proj(u,i,j)
            enddo
           call NORM(u,i,s)
          enddo
        endif

10    continue

      write(iunit,101)x(1),x(2),x(3)

 100  continue
      do 20 i2=1,irun*istep
c       add noise
        if (sn.gt.0.0) then
          call gauss(r1,r2,iseed1,iseed2)
          x(1)=x(1)+r1*sn
          x(2)=x(2)+r2*sn
          call gauss(r1,r2,iseed1,iseed2)
          x(3)=x(3)+r1*sn
          xav=xav+x(1)
          xsq=xsq+x(1)**2
          rsq=rsq+r1*r1
        endif
        call RUKU(3,x,u,dh,bb,ss,rr)
        if (mod(i2,istep).eq.0) write(iunit,101)x(1),x(2),x(3)
        if (mod(i2,ireno).eq.0) then 
c         Gram Schmidt Orthonormierung
          call norm(u,1,s)
          sliap(1)=sliap(1)+log(s)
          do i=2,3
            do j=1,i-1
             call proj(u,i,j)
            enddo
            call NORM(u,i,s)
            sliap(i)=sliap(i)+log(s)
          enddo
        endif

 20   continue

      if (sn.gt.0.0) then
        xav=xav/irun/istep
        xsq=xsq/irun/istep
        rsq=rsq/irun/istep
        rlevel=sqrt(rsq)/sqrt(xsq-xav*xav)*100.
        if (iverb.eq.1) 
     .   print*,'noise level in percent of x-coordinate',rlevel
      endif
      if (iverb.eq.1) then
       write(istderr(),*)
       write(istderr(),*)'Lyapunov exponents [1/unit time]'
       do i=1,3
        write(istderr(),*)real(sliap(i)/totaltime)
       enddo
      endif

 101  format(2x,3f10.3)

      stop
      end

      subroutine FORCE(x,ff,dh,bb,ss,rr)
      real*8 x(3),ff(3),dh,bb,ss,rr

        ff(1)=dh*ss*(x(2)-x(1))
        ff(2)=dh*(x(1)*(-x(3)+rr)-x(2))
        ff(3)=dh*(x(1)*x(2)-bb*x(3))

      return
      end

      subroutine LFORCE(x,u,fl,dh,bb,ss,rr)
      real*8 x(3),u(3,3),dh,fl(3,3),bb,ss,rr

       do j=1,3
         fl(j,1)=dh*ss*(u(j,2)-u(j,1))
         fl(j,2)=dh*(u(j,1)*(rr-x(3))-x(1)*u(j,3)-u(j,2))
         fl(j,3)=dh*(u(j,1)*x(2)+x(1)*u(j,2)-bb*u(j,3))
       enddo
      return
      end

      subroutine RUKU(n,x,u,dh,bb,ss,rr)
c     4th-order Runge Kutta
c     initial point x
c     final point y
c     stepsize dh
c     add subroutine force
      
      implicit real*8 (a-h,o-z)
      real*8 x(3),ff1(3),ff2(3),ff3(3),ff4(3),dummy(3)
      real*8 u(3,3),fl1(3,3),fl2(3,3),fl3(3,3),fl4(3,3)
      real*8 dl(3,3)

      call force(x,ff1,dh,bb,ss,rr)
      call LFORCE(x,u,fl1,dh,bb,ss,rr)

      do i=1,n
      dummy(i)=ff1(i)*.5d0+x(i)
        do j=1,3
        dl(i,j)=fl1(i,j)*.5d0+u(i,j)
        enddo
      enddo

      call force(dummy,ff2,dh,bb,ss,rr)
      call LFORCE(dummy,dl,fl2,dh,bb,ss,rr)

      do i=1,n
      dummy(i)=ff2(i)*.5d0+x(i)
        do j=1,3
        dl(i,j)=fl2(i,j)*.5d0+u(i,j)
        enddo
      enddo

      call force(dummy,ff3,dh,bb,ss,rr)
      call LFORCE(dummy,dl,fl3,dh,bb,ss,rr)

      do i=1,n
      dummy(i)=ff3(i)+x(i)
        do j=1,3
        dl(i,j)=fl3(i,j)+u(i,j)
        enddo
      enddo

      call force(dummy,ff4,dh,bb,ss,rr)
      call LFORCE(dummy,dl,fl4,dh,bb,ss,rr)

      do i=1,n
      x(i)=x(i)+ff1(i)/6.d0+ff2(i)/3.d0+ff3(i)/3.d0+ff4(i)/6.d0
        do j=1,3
        u(i,j)=u(i,j)+fl1(i,j)/6.d0+fl2(i,j)/3.d0+fl3(i,j)/3.d0
     +               +fl4(i,j)/6.d0
        enddo
      enddo

      return
      end

      subroutine NORM(u,i,s)
      real*8 u(3,3),s

      s=0.d0
      do 10 j=1,3
10    s=s+u(i,j)**2
      s=sqrt(s)
      si=1.d0/s
      do 20 j=1,3
20    u(i,j)=u(i,j)*si
      return
      end

      subroutine PROJ(u,i,j)
      real*8 u(3,3),s
      s=0.d0
      do 10 k=1,3
10      s=s+u(i,k)*u(j,k)
      do 20 k=1,3
20      u(i,k)=u(i,k)-s*u(j,k)
      return
      end

c>-------------------------------------------------------
      subroutine gauss(r1,r2,iseed1,iseed2)

      real*8 r1,r2,p,phi,r
      pii=8.d0*atan(1.d0)

      call RANDOM1(p,iseed1)
      call RANDOM1(phi,iseed2)
       
      phi=phi*pii
      r=sqrt(-log(1.d0-p)*2.d0)

      r1=r*sin(phi)
      r2=r*cos(phi)
      return
      end
c>-------------------------------------------------------
      subroutine RANDOM1(r,iseed)
c
c     random number generator of Park & Miller
c     random numbers in [0,1] !!!
      real*8 r
      integer*8 ia,im,ix
      ia=7**5
      im=2147483647
      ix=iseed
      ix=mod(ia*ix,im)
      r=dfloat(ix)/dfloat(im)
      iseed=ix
      return
      end
c>------------------------------------------------------------------
      subroutine usage()
c usage message

      call whatineed(
     .   "-l# [-f# -r# -R# -S# -B# -o outfile -x# -V# -h]")
      call popt("l","length of trajectory x,y,z")
      call popt("f","sample points per unit time [100]")
      call popt("r","absolute noise amplitute [0]")
      call popt("R","parameter r [28]")
      call popt("S","parameter sigma [10]")
      call popt("B","parameter b [8/3]")
      call popt("x","transient discarded [100 units of time]")
c      call popt("L","if present: compute Lyapunov exponents")
      call pout("lorenz.dat")
      call pall()
      stop
      end



