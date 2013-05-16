c Create surrogate data
c Copyright (C) T. Schreiber (1999)

      parameter(nx=1000000)
      dimension xx(nx), x(nx), y(nx), xamp(nx), xsort(nx), list(nx)
      character*72 file, fout
      data nsur/1/, imax/-1/
      external rand

      call whatido("Create Surrogate data")
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      nsur=min(999,ican("n",nsur))
      imax=ican("i",imax)
      ispec=ican("S",0)
      r=rand(sqrt(fcan("I",0.0)))
      iverb=igetv(15)
      isout=igetout(fout,iexv(iverb,1))

      call nthstring(1,file)
      call readfile(nmax,xx,nexcl,jcol,file,iexv(iverb,1))
      nmaxp=nless(nmax)
      if(nmaxp.ne.nmax.and.iexv(iverb,1).ne.0) 
     .   write(istderr(),*) "surrogates: using first ", nmaxp
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_surr")
      if(nsur.gt.1.and.isout.eq.1) call suffix(fout,"_000")

      do 10 isur=1,nsur
         if(nsur.gt.1.and.isout.eq.1) 
     .      write(fout(index(fout," ")-3:72),'(i3.3)') isur
         do 20 n=1,nmaxp
            x(n)=xx(n)
            y(n)=x(n)
            xamp(n)=x(n)
 20         xsort(n)=x(n)
         do 30 n=1,nmaxp
            nper=min(int(rand(0.0)*nmaxp)+1,nmaxp)
            h=x(n)
            x(n)=x(nper)
 30         x(nper)=h
         call store_spec(nmaxp,xamp,0)
         call sort(nmaxp,xsort,list)
         it=-1
         dspec=r1mach(2)
 1       it=it+1
         do 40 n=1,nmaxp
 40         y(n)=x(n)
         ds0=dspec
         dspec=totospec(nmaxp,xamp,y)
         if(imax.ge.0.and.it.ge.imax) goto 2
         call todist(nmaxp,xsort,y,x)
         if(dspec.lt.ds0) goto 1
 2       continue
         if(ispec.gt.0) then
            call writefile(nmaxp,y,fout,iexv(iverb,1))
         else
            call writefile(nmaxp,x,fout,iexv(iverb,1))
         endif
 10      if(iexv(iverb,2).ne.0) write(istderr(),*) 
     .      fout(1:index(fout," ")), ' (', it, 
     .      ' iterations, relative discrepancy ', dspec,   ')'
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-n# -i# -S -I# -o outfile -l# -x# -c# -V# -h] file")
      call popt("n","number of surrogates (1)")
      call popt("i","number of iterations (until no change)")
      call popt("S","make spectrum exact rather than distribution")
      call popt("I","seed for random numbers")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_surr(_nnn)")
      call pall()
      stop
      end

      function totospec(nmax,a,x)
      parameter(nx=1000000)
      dimension x(nmax), a(nmax), w(nx), w1(nx), w2(nx), iw(15)

      if(nmax.gt.nx) stop "totospec: make nx larger."
      do 10 n=1,nmax
 10      w(n)=x(n)
      call rffti1(nmax,w2,iw)  
      call rfftf1(nmax,x,w1,w2,iw)
      do 20 n=1,nmax
 20      x(n)=x(n)/real(nmax)
      x(1)=sqrt(a(1))
      do 30 n=2,(nmax+1)/2
         ab=a(2*n-2)/(x(2*n-2)**2+x(2*n-1)**2)
         x(2*n-2)=x(2*n-2)*sqrt(ab)
 30      x(2*n-1)=x(2*n-1)*sqrt(ab)
      if(mod(nmax,2).eq.0) x(nmax)=sqrt(a(nmax))
      call rfftb1(nmax,x,w1,w2,iw)
      totospec=0
      do 40 n=1,nmax
 40      totospec=totospec+(x(n)-w(n))**2
      totospec=sqrt(totospec/nmax)
      end


