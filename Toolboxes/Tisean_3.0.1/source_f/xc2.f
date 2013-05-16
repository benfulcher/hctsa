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
c=========================================================================
c
c   cross-correlation integral xc2
c   see  H. Kantz, Phys.Rev.E49, 5091 (1994)
c
c   authors: T. Schreiber & H. Kantz (1998)
c   multivariate version: H. Kantz (Jan 2007)
c
c=========================================================================
c
      parameter(nx=1000000,me=30,meps=1000,mx=10)
      dimension x(nx,mx), c(me,meps), eps(meps), mdeps(meps), icol(mx)
      dimension y(nx,mx)
      integer mlist(2)
      character*72 file1, file2, fout
      data ipmin/1000/, res/2./, eps0/1e-30/, epsm/1e30/, id0/1/
      data iverb/1/

c=======================================================================
c assume two input files (file1, file2) 
c - with identical structure in terms of colums
c - with identical embedding spaces 
c - with individual length and exclusions: l, L, x, X
c - multivariate data: maximum length nx, maximum dimension mx
c
c   norm: max-norm 
c
c   no rescaling, since xc2 does not make sense if datasets are 
c                           not used in their original scalings
c=================================================================

      call whatido("cross correlation sum of two data sets",iverb)
      id=ican("d",id0)
      mmax=2
      mdim=1

      call imcan("M",2,mc,mlist)
      if (mc.ge.2) then
       mmax=mlist(2)
       mdim=mlist(1)
       if (mmax*mdim.lt.2) stop 'Increase embedding dimension'
       if (mc.gt.2) print*, 'extra arguments of -m ignored'
      endif

      ntmin=0
      ncmin=ican("n",1000)
      ipmin=ican("N",ipmin)
      res=fcan("#",res)
      feps=2**(1./res)
      eps0=fcan("r",eps0)
      epsm=fcan("R",epsm)
      nmaxx=ican("l",nx)
      nmaxy=ican("L",nx)
      nexcl1=ican("x",0)
      nexcl2=ican("X",0)

      call columns(mc,mx,icol)

      if (mc.gt.0.and.mc.ne.mdim) stop 'improper number of columns'
      isout=igetout(fout,0)
      if(fout.eq." ") isout=1

      call nthstring(1,file1)
      if(file1.eq."-") stop "first input file name missing"
      call xreadfile(nmaxx,mdim,nx,x,nexcl1,icol,file1,iverb)
      call nthstring(2,file2)
      if(file2.eq."-") stop "second input file name missing"
      call xreadfile(nmaxy,mdim,nx,y,nexcl2,icol,file2,iverb)

      if(isout.eq.1) then 
        call addsuff(fout,file1,"_")
        call addsuff(fout,fout,file2)
        call addsuff(fout,fout,"_xc2")
      endif

      epsmax=0.
      do imx=1,mmax
       call minmax(nmaxx,x(1,imx),xmin,xmax)
       epsmax=1.001*max(xmax-xmin,epsmax)
      enddo
      do imx=1,mmax
       call minmax(nmaxy,y(1,imx),xmin,xmax)
       epsmax=1.001*max(xmax-xmin,epsmax)
      enddo

      neps=0

      do 10 epsl=log(min(epsm,epsmax)),log(eps0),-log(feps)
         neps=neps+1
         if(neps.gt.meps) stop "xc2: make meps larger"
         eps(neps)=exp(epsl)
         do 20 m=1,mmax*mdim
 20         c(m,neps)=0
         if (mdim.eq.1) then
          call crosscor(nmaxx,x,nmaxy,y,eps(neps)
     .                 ,id,mmax,c(1,neps),ncmin,ipmin)
         else
          call mcrosscor(nmaxx,x,nmaxy,y,eps(neps),
     .      id,mmax,mdim,c(1,neps),ncmin,ipmin)
         endif
         mdd=mmax*mdim
         mdd1=max(2,mdim)
         do 30 m=mdd1,mdd
 30         if(c(m,neps).eq.0.) goto 1
         m=mdd+1
 1       mdd=m-1
         if(mdd.eq.mdim-1) stop
         mdeps(neps)=mdd
         call outfile(fout,iunit,iverb)
         do 40 m=mdd1,mdeps(1)
            write(iunit,'(4h#m= ,i5)') m
            do 50 nn=1,neps
               if(mdeps(nn).lt.m) goto 2
 50            write(iunit,*) eps(nn), c(m,nn) 
 2          write(iunit,'()')
 40         write(iunit,'()')
         close(iunit)
 10      write(istderr(),*) eps(neps), mdd, c(mdd,neps)
      stop
      end
c>--------------------------------------------------------------------
      subroutine usage()
c usage message

      call whatineed(
     .   "-M#,# [-d# -n# -N# -## -r# -R#"//
     .   " -o outfile -l# -x# -L# -X# -c#[,#] -V# -h] file1 file2")
      call popt("M",
     ."# of components, maximal embedding dimension [1,2]")
      call popt("d","delay [1]")
      call popt("n","minimal number of center points [1000]")
      call popt("N","maximal number of pairs [1000]")
      call popt("#","resolution, values per octave [2]")
      call popt("r",
     .   "minimal scale to be probed (as long as pairs found)")
      call popt("R","maximal scale to be probed [xmax-xmin]")
      call popt("l","length of time series 1 to be read [all data]")
      call popt("x","# of initial lines of 1 to be skipped [0]")
      call popt("L","length of time series 2 to be read [all data]")
      call popt("X","# of initial lines of 2 to be skipped [0]")
      call popt("c",
     ."columns to be read [1,2,3,.., # of components]")
      call pout("file1_file2_xc2")
      call pall()
      stop
      end
c>--------------------------------------------------------------------
      subroutine crosscor(nmaxx,x,nmaxy,y,eps,id,m,c,ncmin,ipmin)
      parameter(im=100,ii=100000000,nx=1000000,mm=30)
      dimension y(nmaxy),x(nmaxx)
      dimension jh(0:im*im),ipairs(mm),c(m),jpntr(nx),nlist(nx)

      if(nmaxx.gt.nx.or.m.gt.mm) stop "crosscor: make mm/nx larger."
      if(nmaxy.gt.nx.or.m.gt.mm) stop "crosscor: make mm/nx larger."

      do 10 i=1,m-1
 10      ipairs(i)=0
      mb=min(m,2)
      call base(nmaxx,x,id,mb,jh,jpntr,eps)
      do 20 n=(m-1)*id+1,nmaxy
         call neigh(nx,y,x,n,nmaxx,id,mb,jh,jpntr,eps,nlist,nfound)
         do 30 nn=1,nfound                   ! all neighbours in two dimensions
            np=nlist(nn)
            if(np.lt.(m-1)*id+1) goto 30
            ipairs(1)=ipairs(1)+1
            do 40 i=mb,m-1
               if(abs(y(n-i*id)-x(np-i*id)).ge.eps) goto 30
 40            ipairs(i)=ipairs(i)+1            ! neighbours in 3..m dimensions
 30         continue
 20      if(n-(m-1)*id.ge.ncmin.and.ipairs(m-1).ge.ipmin) goto 1
      n=n-1
 1    s=real(n-(m-1)*id)*real(nmaxx-(m-1)*id) ! normalisation
      do 50 i=1,m-1
 50      if(s.gt.0.) c(i+1)=ipairs(i)/s
      end
c>--------------------------------------------------------------------
      subroutine mcrosscor(nmaxx,x,nmaxy,y,eps,id,m,mdim,c,ncmin,ipmin)

      parameter(im=100,nx=1000000,mm=30,mx=10)

      dimension y(nx,mx),x(nx,mx)
      dimension jh(0:im*im),ipairs(mm),c(m*mdim),jpntr(nx),nlist(nx)
      dimension vx(mm)

      if(nmaxx.gt.nx.or.m.gt.mm) stop "mcrosscor: make mm/nx larger."
      if(nmaxy.gt.nx.or.m.gt.mm) stop "mcrosscor: make mm/nx larger."

      if (m*mdim.gt.mm) stop 'embedding x spatial dimension < 30 !'
      do 10 i=1,m*mdim
 10      ipairs(i)=0

      call mbase(nmaxy,mdim,nx,y,id,1,jh,jpntr,eps)
      do 20 n=(m-1)*id+1,nmaxx
         do ii=1,mdim
          vx(ii)=x(n,ii)
         enddo
         call mneigh2(nmaxy,mdim,y,nx,vx,jh,jpntr,eps,nlist,nfound)
         do 30 nn=1,nfound               ! all neighbours in mdim dimensions
            np=nlist(nn)
            if(np.lt.(m-1)*id+1) goto 30
            ipairs(mdim)=ipairs(mdim)+1
            do 40 i=1,m-1
             do 41 iim=1,mdim
               if(abs(x(n-i*id,iim)-y(np-i*id,iim)).ge.eps) goto 30
               idim=i*mdim+iim
               ipairs(idim)=ipairs(idim)+1  ! neighbours mdim+1..m dimensions
 41           continue
 40          continue
 30         continue
 20      if(n-(m-1)*id.ge.ncmin.and.ipairs(m*mdim).ge.ipmin) goto 1
      n=n-1
 1    s=real(n-(m-1)*id)*real(nmaxy-(m-1)*id) ! normalisation
      do 50 i=mdim,mdim*m
 50      if(s.gt.0.) c(i)=ipairs(i)/s

      return
      end
c>---------------------------------------------------------------------






