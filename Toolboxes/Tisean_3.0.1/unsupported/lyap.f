c Lyapunov exponent
c see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997)
c Copyright (C) T. Schreiber, H. Kantz (1997)

      subroutine lyap(nmax,y,id,m,eps,ifu,s,nmin,nfmin,ncmin)
      parameter (im=100,ii=100000000,nx=1000000,ifum=10000,tiny=1e-20)
      dimension y(nmax),s(0:ifu),sh(0:ifum),jh(0:im*im),
     .  jpntr(nx),nlist(nx)

      if(nmax.gt.nx.or.ifu.gt.ifum) stop "lyap: make nx/ifum larger."
      call base(nmax-ifu,y,id,m,jh,jpntr,eps)
      do 10 i=0,ifu
 10     s(i)=0
      nc=0
      do 20 n=(m-1)*id+1, nmax-ifu                           ! reference points
        call neigh(nmax-ifu,y,n,nmax,id,m,jh,jpntr,eps,nlist,nfound)
        do 30 i=0,ifu
 30       sh(i)=0
        nf=0
        do 40 nn=1,nfound
          np=nlist(nn)
          if(abs(n-np).le.nmin) goto 40
          nf=nf+1
          do 50 i=0,ifu                                     ! average distances
 50         sh(i)=sh(i)+abs(y(n+i)-y(np+i))                  
 40       continue
        if(nf.ge.nfmin) then              ! enough neighbours closer $\epsilon$
          nc=nc+1
          do 60 i=0,ifu                     ! average log of averaged distances
 60         s(i)=s(i)+log(max(sh(i)/nf,tiny))  
        endif
 20     if(nc.eq.ncmin) goto 1
 1    if(nc.eq.0) return                     ! no points with enough neighbours
      do 70 i=0,ifu
 70     s(i)=s(i)/nc
      write(istderr(),*) 'tried ', n,' reference points, found ', nc
      end
