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
c   
c   d1 with finite sample correction following Grassberger
c   subroutine for c1
c
c===========================================================================
      subroutine d1(nmax,mmax,nxx,y,id,m,ncmin,pr,pln,eln,nmin,kmax)
      parameter(im=100,ii=100000000,nx=100000,tiny=1e-20) 
      dimension y(nxx,mmax),jh(0:im*im),ju(nx),d(nx),jpntr(nx),
     .   nlist(nx),nwork(nx)
      external rand

      if(nmax.gt.nx) stop "d1: make nx larger."
      mt=(m-1)/mmax+1
      ncomp=nmax-(mt-1)*id
      kpr=int(exp(pr)*(ncomp-2*nmin-1))+1
      k=int(exp(pln)*(ncomp-2*nmin-1))+1
      if(k.gt.kmax) then
         ncomp=real(ncomp-2*nmin-1)*real(kmax)/k+2*nmin+1
         k=kmax
      endif         
      pln=psi(k)-log(real(ncomp-2*nmin-1))
      if(k.eq.kpr) return
      write(istderr(),*) 'Mass ', exp(pln),': k=', k, ', N=', ncomp 
      call rms(nmax,y,sc,sd)
      eps=exp(pln/m)*sd
      do 10 i=1,nmax-(mt-1)*id
 10      ju(i)=i+(mt-1)*id
      do 20 i=1,nmax-(mt-1)*id
         iperm=min(int(rand(0.0)*nmax-(mt-1)*id)+1,nmax-(mt-1)*id)
         ih=ju(i)
         ju(i)=ju(iperm)
 20      ju(iperm)=ih
      iu=ncmin
      eln=0
 1    call mbase(ncomp+(mt-1)*id,mmax,nxx,y,id,m,jh,jpntr,eps)
      iunp=0
      do 30 nn=1,iu                                           ! find neighbours
         n=ju(nn)
         call mneigh(nmax,mmax,nxx,y,n,nmax,id,m,jh,jpntr,eps,
     .      nlist,nfound)
         nf=0
         do 40 ip=1,nfound
            np=nlist(ip)
            nmd=mod(abs(np-n),ncomp)
            if(nmd.le.nmin.or.nmd.ge.ncomp-nmin) goto 40  ! temporal neighbours
            nf=nf+1
            dis=0
            mcount=0
            do 50 i=mt-1,0,-1
               do 50 is=1,mmax
                  mcount=mcount+1
                  if(mcount.gt.m) goto 2
 50               dis=max(dis,abs(y(n-i*id,is)-y(np-i*id,is)))
 2          d(nf)=dis
 40         continue
         if(nf.lt.k) then
            iunp=iunp+1                                   ! mark for next sweep
            ju(iunp)=n
         else
            e=which(nf,d,k,nwork)
            eln=eln+log(max(e,tiny))
         endif
 30      continue
      iu=iunp
      eps=eps*sqrt(2.)
      if(iunp.ne.0) goto 1
      eln=eln/(ncmin-(mt-1)*id)
      end

c digamma function
c Copyright (C) T. Schreiber (1998)

      function psi(i)
      dimension p(0:20)
      data p/0., 
     .  -0.57721566490,  0.42278433509,  0.92278433509,  1.25611766843,
     .   1.50611766843,  1.70611766843,  1.87278433509,  2.01564147795,
     .   2.14064147795,  2.25175258906,  2.35175258906,  2.44266167997,
     .   2.52599501330,  2.60291809023,  2.67434666166,  2.74101332832,
     .   2.80351332832,  2.86233685773,  2.91789241329,  2.97052399224/

      if(i.le.20) then
         psi=p(i)
      else
         psi=log(real(i))-1/(2.*i)
      endif
      end








