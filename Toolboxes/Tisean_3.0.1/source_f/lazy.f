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
c   lazy.f
c   simple nonlinear noise reduction
c   see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997,2004)
c   author T. Schreiber (1998)
c===========================================================================
      parameter(nx=1000000)
      dimension x(nx), x0(nx), xc(nx)
      character*72 file, fout
      data eps/0./, frac/0./, imax/1/
      data iverb/1/

      call whatido("simple nonlinear noise reduction",iverb)
      m=imust("m")
      eps=fcan("r",eps)
      frac=fcan("v",frac)
      imax=ican("i",imax)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)
      if(eps.eq.0.and.frac.eq.0.) call usage()

      call nthstring(1,file)
      call readfile(nmax,x,nexcl,jcol,file,iverb)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_l")
      call rms(nmax,x,sc,sd)
      if(frac.gt.0) eps=sd*frac
      do 10 n=1,nmax
 10      x0(n)=x(n)
      do 20 it=1,imax
         call nrlazy(nmax,x,xc,m,eps)
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
         if(eps.eq.0.) then
            if(iv_io(iverb).eq.1) write(istderr(),*) 
     .      'Zero correction, finished'
            stop
         endif
 20      if(iv_io(iverb).eq.1) write(istderr(),*) 
     .      'New diameter of neighbourhoods is ', eps
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-m# [-r# | -v#] [-i# -o outfile -l# -x# -c# -V# -h] file")
      call ptext("either -r or -v must be present")
      call popt("m","embedding dimension")
      call popt("r","absolut radius of neighbourhoods")
      call popt("v","same as fraction of standard deviation")
      call popt("i","number of iterations (1)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_lc,  file_lcc (etc.)")
      call pall()
      stop
      end

      subroutine nrlazy(nmax,y,yc,m,eps)
      parameter(im=100,ii=100000000,nx=1000000) 
      dimension y(nmax),yc(nmax),jh(0:im*im),jpntr(nx),nlist(nx)

      if(nmax.gt.nx) stop "nrlazy: make nx larger."
      call base(nmax,y,1,m,jh,jpntr,eps)
      do 10 n=1,nmax
 10      yc(n)=y(n)   
      do 20 n=m,nmax           
         call neigh(nmax,y,y,n,nmax,1,m,jh,jpntr,eps,nlist,nfound)
         av=0
         do 30 nn=1,nfound            
 30         av=av+y(nlist(nn)-(m-1)/2)              ! average middle coordinate
 20      yc(n-(m-1)/2)=av/nfound
      end
