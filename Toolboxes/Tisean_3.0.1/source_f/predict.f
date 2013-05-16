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
c   simple nonlinear prediction, fast neighbour search
c   see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997,2004)
c   author T. Schreiber (1998)
c===========================================================================
      parameter(nx=1000000)
      dimension x(nx), y(nx)
      character*72 file, fout
      data eps/0./, frac/0./, ifc/1/
      data iverb/1/

      call whatido("prediction with locally constant fits",iverb)
      id=imust("d")
      m=imust("m")
      eps=fcan("r",eps)
      frac=fcan("v",frac)
      ifc=ican("s",ifc)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)
      if(eps.eq.0.and.frac.eq.0.) call usage()

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iverb)
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_pred")
         call rms(nmax,x,sc,sd)
         if(frac.gt.0) eps=sd*frac
         iun=istdout()
         if(fout.eq." ") iun=istderr()
         write(iun,*) "err: ", fcerror(nmax,x,y,m,id,ifc,eps), 
     .      " "//file(1:index(file," ")-1)
 10      call writefile(nmax,y,fout,iverb)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# -m# [-r# | -v#]"//
     .   " [-s# -o outfile -l# -x# -c# -V# -h] file(s)")
      call ptext("either -r or -v must be present")
      call popt("d","delay")
      call popt("m","embedding dimension")
      call popt("r","absolute radius of neighbourhoods")
      call popt("v","same as fraction of standard deviation")
      call popt("s","time steps ahead forecast (one step)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_pred")
      call pall()
      stop
      end

      function fcerror(nmax,y,yp,m,id,ifc,eps)
      parameter(im=100,ii=100000000,nx=1000000) 
      dimension y(nmax),yp(nx),jh(0:im*im),jpntr(nx),nlist(nx)

      if(nmax.gt.nx) stop "fcerror: make nx larger."
      call base(nmax-ifc,y,id,m,jh,jpntr,eps)
      fcerror=0

      call rms(nmax,y,sx,sd)
      do 10 n=1,(m-1)*id+ifc
 10      yp(n)=sx
      do 20 n=(m-1)*id+1,nmax-ifc           
         call neigh(nmax,y,y,n,nmax,id,m,jh,jpntr,eps,nlist,nfound)
         av=0
         do 30 nn=1,nfound            
 30         if(nlist(nn).ne.n) av=av+y(nlist(nn)+ifc) 
         if(nfound.gt.1) then
            yp(n+ifc)=av/(nfound-1)
         else
            yp(n+ifc)=sx
         endif
 20      fcerror=fcerror+(y(n+ifc)-yp(n+ifc))**2
      fcerror=sqrt(fcerror/(nmax-ifc-(m-1)*id))
      end

