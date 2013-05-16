c correlation integral c2 
c see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997)
c Copyright (C) T. Schreiber (1997)

      parameter(nx=1000000,me=30,meps=1000)
      dimension x(nx), c(me,meps), eps(meps), mdeps(meps)
      character*72 file, fout
      data ipmin/1000/, res/2./, eps0/1e-30/, epsm/1e30/

      call whatido("correlation sum (see also: c2naive)")
      id=imust("d")
      mmax=max(imust("M"),2)
      ntmin=imust("t")
      ncmin=imust("n")
      ipmin=ican("N",ipmin)
      res=fcan("#",res)
      feps=2**(1./res)
      eps0=fcan("r",eps0)
      epsm=fcan("R",epsm)
      nmaxx=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      iverb=igetv(1)
      isout=igetout(fout,iexv(iverb,1))
      if(fout.eq." ") isout=1

      do 10 ifi=1,nstrings()
         call nthstring(ifi,file)
         nmax=nmaxx
         call readfile(nmax,x,nexcl,jcol,file,iexv(iverb,1))
         if(file.eq."-") file="stdin"
         if(isout.eq.1) call addsuff(fout,file,"_c2")
         call minmax(nmax,x,xmin,xmax)
         md=mmax
         neps=0
         do 20 epsl=log(min(epsm,xmax-xmin)),log(eps0),-log(feps)
            neps=neps+1
            if(neps.gt.meps) goto 10    ! done, next file
            eps(neps)=exp(epsl)
            do 30 m=2,md
 30            c(m,neps)=0
            call correl(nmax,x,eps(neps),id,md,c(1,neps),
     .         ntmin,ncmin,ipmin)
            do 40 m=2,md
 40            if(c(m,neps).eq.0.) goto 1
 1             md=m-1
            if(md.eq.1) goto 10     ! done, next file
            mdeps(neps)=md
            call outfile(fout,iunit,iexv(iverb,1))
            do 50 m=2,mdeps(1)
               write(iunit,'(4h#m= ,i5)') m
               do 60 nn=1,neps
                  if(mdeps(nn).lt.m) goto 2
 60               write(iunit,*) eps(nn), c(m,nn) 
 2             write(iunit,'()')
 50            write(iunit,'()')
            close(iunit)
 20         write(istderr(),*) eps(neps), md, c(md,neps)
 10      continue
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# -M# -t# -n# [-N# -## -r# -R#"//
     .   " -o outfile -l# -x# -c# -V# -h] file(s)")
      call popt("d","delay")
      call popt("M","maximal embedding dimension (at least 2)")
      call popt("t","minimal time separation")
      call popt("n","minimal number of center points")
      call popt("N","maximal number of pairs (1000)")
      call popt("#","resolution, values per octave (2)")
      call popt("r",
     .   "minimal length to be probed (as long as pairs found)")
      call popt("R","maximal length to be probed (xmax-xmin)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1, or file,#)")
      call pout("file_c2")
      call pall()
      stop
      end

      subroutine correl(nmax,y,eps,id,m,c,nmin,ncmin,ipmin)
      parameter(im=100,ii=100000000,nx=1000000,mm=30)
      dimension y(nmax),jh(0:im*im),ipairs(mm),c(m),jpntr(nx),nlist(nx)

      if(nmax.gt.nx.or.m.gt.mm) stop "correl: make mm/nx larger."
      do 10 i=1,m-1
 10      ipairs(i)=0
      call base(nmax,y,id,2,jh,jpntr,eps)
      do 20 n=nmin+(m-1)*id+1,nmax
         call neigh(nmax,y,n,n-nmin,id,2,jh,jpntr,eps,nlist,nfound)
         ipairs(1)=ipairs(1)+nfound
         do 30 nn=1,nfound                   ! all neighbours in two dimensions
            np=nlist(nn)
            if(np.lt.(m-1)*id+1) goto 30
            do 40 i=2,m-1
               if(abs(y(n-i*id)-y(np-i*id)).ge.eps) goto 30
 40            ipairs(i)=ipairs(i)+1            ! neighbours in 3..m dimensions
 30         continue
 20      if(n-nmin-(m-1)*id.ge.ncmin.and.ipairs(m-1).ge.ipmin) goto 1
 1    s=real(n-nmin-(m-1)*id+1)*real(n-nmin-(m-1)*id)/2.        ! normalisation
      do 50 i=1,m-1
 50      c(i+1)=ipairs(i)/s
      end
