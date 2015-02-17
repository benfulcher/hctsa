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
c   Create multivariate surrogate data
c   author T. Schreiber (1999)
c===========================================================================
      parameter(nx=100000,mx=20)
      dimension xx(nx,mx), x(nx,mx), y(nx,mx), xamp(nx,mx), 
     .   xsort(nx,mx), list(nx), icol(mx), rwork(nx)
      character*72 file, fout
      data nsur/1/, imax/-1/
      external rand
      data iverb/15/

      call whatido("Create Multivariate Surrogate data",iverb)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      nsur=min(999,ican("n",nsur))
      imax=ican("i",imax)
      ispec=lopt("S",1)
      r=rand(sqrt(abs(fcan("I",0.0))))
      mcmax=ican("m",0)
      call columns(mc,mx,icol)
      if(mcmax.eq.0) mcmax=max(1,mc)
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call xreadfile(nmax,mcmax,nx,xx,nexcl,icol,file,iverb)
      nmaxp=nless(nmax)
      if(nmaxp.ne.nmax.and.iv_io(iverb).eq.1) 
     .   write(istderr(),*) "surrogates: using first ", nmaxp
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_surr")
      if(nsur.gt.1.and.isout.eq.1) call suffix(fout,"_000")

      do 10 isur=1,nsur
         if(nsur.gt.1.and.isout.eq.1) 
     .      write(fout(index(fout," ")-3:72),'(i3.3)') isur
         do 20 m=1,mcmax
            do 30 n=1,nmaxp
               x(n,m)=xx(n,m)
               y(n,m)=x(n,m)
               xamp(n,m)=x(n,m)
 30            xsort(n,m)=x(n,m)
            call store_spec(nmaxp,xamp(1,m),0)
            call sort(nmaxp,xsort(1,m),list)
            do 40 n=1,nmaxp
 40            rwork(n)=rand(0.0)
            call rank(nmaxp,rwork,list)
 20         call index2sort(nmaxp,x(1,m),list)
         it=-1
         dspec=r1mach(2)
 1       it=it+1
         do 50 m=1,mcmax
            do 50 n=1,nmaxp
 50            y(n,m)=x(n,m)
         ds0=dspec
         dspec=toxspec(nmaxp,mcmax,nx,xamp,y)
         if(imax.ge.0.and.it.ge.imax) goto 2
         do 60 m=1,mcmax
 60         call todist(nmaxp,xsort(1,m),y(1,m),x(1,m))
         if(dspec.lt.ds0) goto 1
 2       continue
         if(ispec.gt.0) then
            call xwritefile(nmaxp,mcmax,nx,y,fout,iverb)
         else
            call xwritefile(nmaxp,mcmax,nx,x,fout,iverb)
         endif
 10      if(iv_surr(iverb).eq.1) write(istderr(),*) 
     .      fout(1:index(fout," ")), ' (', it, 
     .      ' iterations, relative discrepancy ', dspec,   ')'
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-n# -i# -S -I# -o outfile -l# -x# -m# -c#[,#] -V# -h] file")
      call popt("n","number of surrogates (1)")
      call popt("i","number of iterations (until no change)")
      call popt("S","make spectrum exact rather than distribution")
      call popt("I","seed for random numbers")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("m","number of columns to be read (1)")
      call popt("c","columns to be read (1)")
      call pout("file_surr(_nnn)")
      call pall()
      call ptext("Verbosity levels (add what you want):")
      call ptext("          1 = input/output" )
      call ptext("          2 = iterations / discrepancy")
      stop
      end

      function toxspec(nmax,mmax,nxx,a,x)
      parameter(nx=100000,mx=20,tol=1e-5)
      dimension x(nxx,mmax), a(nxx,mmax), w(nx,mx), w1(nx), 
     .   w2(nx), iw(15), goal(mx)

      if(nmax.gt.nx.or.mmax.gt.mx) stop "toxspec: make nx/mx larger."
      call rffti1(nmax,w2,iw)  
      do 10 m=1,mmax
         do 20 n=1,nmax
 20         w(n,m)=x(n,m)
         call rfftf1(nmax,x(1,m),w1,w2,iw)
         do 30 n=1,nmax
 30         x(n,m)=x(n,m)/real(nmax)
         x(1,m)=sqrt(a(1,m))
         do 40 n=2,(nmax+1)/2
            pha=atan2(x(2*n-1,m),x(2*n-2,m))
            x(2*n-2,m)=sqrt(a(2*n-2,m))
 40         x(2*n-1,m)=pha
 10      if(mod(nmax,2).eq.0) x(nmax,m)=sqrt(a(nmax,m))
      if(mmax.gt.1) then
         do 50 n=2,(nmax+1)/2
            do 60 m=1,mmax
 60            goal(m)=x(2*n-1,m)-a(2*n-1,m)
            alpha=alp(mmax,goal)
            do 50 m=1,mmax
 50            x(2*n-1,m)=alpha+a(2*n-1,m)
      endif
      do 70 m=1,mmax
         do 80 n=2,(nmax+1)/2
            c=x(2*n-2,m)*cos(x(2*n-1,m))
            s=x(2*n-2,m)*sin(x(2*n-1,m))
            x(2*n-1,m)=s
 80         x(2*n-2,m)=c
 70      call rfftb1(nmax,x(1,m),w1,w2,iw)
      toxspec=0
      do 90 m=1,mmax
         do 90 n=1,nmax
 90         toxspec=toxspec+(x(n,m)-w(n,m))**2
      toxspec=sqrt((toxspec/nmax)/mmax)
      end

      function alp(mmax,goal)
      dimension goal(mmax)
      data pi/3.1415926/

      f1=0
      f2=0
      do 10 m=1,mmax
         f1=f1+cos(goal(m))
 10      f2=f2+sin(goal(m))
      alp=atan2(f2,f1)
      scos=0
      do 20 m=1,mmax
 20      scos=scos+cos(alp-goal(m))
      if(scos.lt.0) alp=alp+pi
      end

      subroutine todist(nmax,dist,x,y)
      parameter(nx=100000)
      dimension x(nmax), dist(nmax), y(nmax), list(nx)

      if(nmax.gt.nx) stop "todist: make nx larger."
      call rank(nmax,x,list)
      do 10 n=1,nmax
 10      y(n)=dist(list(n))
      end

