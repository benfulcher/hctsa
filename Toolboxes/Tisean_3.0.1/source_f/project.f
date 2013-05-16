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
c   nonlinear noise reduction
c   see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997,2004)
c   authors T. Schreiber, H. Kantz, R. Hegger (1998) based on earlier versions
c===========================================================================
      parameter(nx=100000)
      dimension x(nx), x0(nx), xc(nx)
      character*72 file, fout
      data imax/1/
      data iverb/1/

      call whatido("nonlinear noise reduction (see also: noise)",iverb)
      m=imust("m")
      nq=m-imust("q")
      eps=fmust("r",eps)
      kmin=imust("k")
      imax=ican("i",imax)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call readfile(nmax,x,nexcl,jcol,file,iverb)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_")

      do 10 n=1,nmax
 10      x0(n)=x(n)
      do 20 it=1,imax
         call clean(nmax,x,xc,m,kmin,nq,eps,iverb)
         if(fout.ne." ".or.isout.eq.1.or.it.eq.imax) then
            if(isout.eq.1) call suffix(fout,"c")
            call outfile(fout,iunit,iverb)
            do 30 n=1,nmax
 30            write(iunit,*) xc(n), x0(n)-xc(n)
            if(iunit.ne.istdout()) close(iunit)
            if(iv_io(iverb).eq.1) call writereport(nmax,fout)
         endif
         eps=0
         do 40 n=1,nmax
            eps=eps+(xc(n)-x(n))**2
 40         x(n)=xc(n)
         eps=sqrt(eps/nmax)
 20      if(iv_io(iverb).eq.1) 
     .      write(istderr(),*) 'New diameter of neighbourhoods is ', eps
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-m# -q# -r# -k# [-i# -o outfile -l# -x# -c# -V# -h] file")
      call popt("m","embedding dimension")
      call popt("q","dimension of manifold")
      call popt("r","radius of neighbourhoods")
      call popt("k","minimal number of neighbours")
      call popt("i","number of iterations (1)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_c, file_cc (etc.)")
      call pall()
      call ptext("Verbosity levels (add what you want):")
      call ptext("          1 = input/output" )
      call ptext("          2 = state of neighbour search")
      write(istderr(),'()') 
      stop
      end

      subroutine clean(nmax,y,yc,m,kmin,nq,d,iverb)
      parameter(im=100,ii=100000000,nx=100000,mm=15,small=0.0001) 
      dimension y(nmax),yc(nmax),r(mm),ju(nx),c(mm,mm),cm(mm),
     .  jh(0:im*im),jpntr(nx),nlist(nx), zcm(mm,nx)

      if(nmax.gt.nx.or.m.gt.mm) stop "clean: make mm/nx larger."
      sr=2*small+m-2                                        ! ${\rm tr}(1/r)=1$
      do 10 i=1,m
         r(i)=sr
 10      if(i.eq.m.or.i.eq.1) r(i)=sr/small
      do 20 i=1,nmax
 20      yc(i)=y(i)
      do 30 istep=1,2
         eps=d
         iu=nmax-m+1
         do 40 i=1,iu
 40         ju(i)=i+m-1
 1       call base(nmax,y,1,m,jh,jpntr,eps)
         iunp=0
         do 50 nn=1,iu                                        ! find neighbours
            n=ju(nn)
            call neigh(nmax,y,y,n,nmax,1,m,jh,jpntr,eps,nlist,nfound)
            if(nfound.lt.kmin) then               ! not enough neighbours found
               iunp=iunp+1                                ! mark for next sweep
               ju(iunp)=n
            else                                      ! fine: enough neighbours
               do 90 i=1,m                              ! centre of mass vector
                  s=0
                  do 100 np=1,nfound
 100                 s=s+y(nlist(np)-m+i)
 90               cm(i)=s/nfound
               if(istep.eq.1) then                  ! just store centre of mass
                  do 110 i=1,m
 110                 zcm(i,n)=cm(i)
               else
                  do 120 i=1,m                ! corrected centre of mass vector
                     s=0
                     do 130 np=1,nfound
 130                    s=s+zcm(i,nlist(np))
 120                 cm(i)=2*cm(i)-s/nfound
                  do 140 i=1,m                      ! compute covariance matrix
                     do 140 j=i,m
                        s=0
                        do 150 np=1,nfound
                           jm=nlist(np)-m
 150                       s=s+(y(jm+i)-cm(i))*(y(jm+j)-cm(j))
                        c(i,j)=r(i)*r(j)*s/nfound
 140                    c(j,i)=c(i,j)
                 call eigen(c,m)               ! find eigenvectors (decreasing)
                 do 160 i=1,m
                    s=0
                    do 170 iq=m-nq+1,m
                       do 170 j=1,m
 170                      s=s+(y(n-m+j)-cm(j))*c(i,iq)*c(j,iq)*r(j)
 160                yc(n-m+i)=yc(n-m+i)-s/r(i)/r(i)
               endif
            endif
 50         continue
         iu=iunp
         if(iv_uncorr(iverb).eq.1) 
     .      write(istderr(),*) "With ", eps, iunp, " uncorrected"
         eps=eps*sqrt(2.)
 30      if(iunp.ne.0) goto 1
      end

c driver for diagonalisation routines
c see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997)
c Copyright (C) T. Schreiber (1997)

      subroutine eigen(c,kk)
      parameter(md=15)
      dimension c(md,md),d(md),w1(md),w2(md),z(md,md)
      if(kk.gt.md) stop "eigen: make md larger."

      call rs(md,kk,c,d,1,z,w1,w2,ierr)
      do 10 i=1,kk
         do 10 j=1,kk
 10         c(i,j)=z(i,kk+1-j)
      end

