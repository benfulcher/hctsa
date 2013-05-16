c Lyapunov exponent
c see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997)
c Copyright (C) T. Schreiber, H. Kantz (1997)

      parameter(nx=1000000,ifum=10000)
      dimension x(nx), s(0:ifum)
      character*72 file, fout
      data ncmin/nx/, ifu/100/, kmin/10/

      call whatido("maximal Lyapunov exponent")
      id=imust("d")
      ntmin=imust("t")
      mmin=max(imust("m"),2)
      mmax=imust("M")
      eps=fmust("r")
      ifu=min(ifum,ican("s",ifu))
      ncmin=ican("n",ncmin)
      kmin=ican("k",kmin)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      iverb=igetv(1)
      isout=igetout(fout,iexv(iverb,1))

      call nthstring(1,file)
      call readfile(nmax,x,nexcl,jcol,file,iexv(iverb,1))
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_lyap")
      call outfile(fout,iunit,iexv(iverb,1))

      do 10 m=mmin,mmax
         call lyap(nmax,x,id,m,eps,ifu,s,ntmin,kmin,ncmin)
         do 20 j=0,ifu
 20         write(iunit,*) j, exp(s(j))
         write(iunit,'()')
         write(iunit,'()')
 10      if(iexv(iverb,1).ne.0) write(istderr(),*) "Finished dim ", m
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# -m# -M# -t# -r# "//
     .   "[-s# -n# -k# -o outfile -l# -x# -c# -V# -h] file")
      call popt("d","delay")
      call popt("m","minimal embedding dimension (at least 2)")
      call popt("M","maximal embedding dimension (at least 2)")
      call popt("t","minimal time separation")
      call popt("r","diameter of initial neighbourhoods")
      call popt("s","time steps over which expansion is followed (100)")
      call popt("n","minimal number of reference points (all)")
      call popt("k",
     .   "minimal number of neighbours for reference points (10)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_lyap")
      call pall()
      stop
      end
