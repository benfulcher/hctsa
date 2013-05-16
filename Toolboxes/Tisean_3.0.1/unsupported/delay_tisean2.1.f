c delay coordinates
C Copyright (C) Thomas Schreiber (1999)

      parameter(nx=1000000)
      dimension x(nx)
      character*72 file, fout
      data id/1/, m/2/
      data iverb/1/

      call whatido("embed using delay coordinates",iverb)
      id=ican("d",id)
      m=ican("m",m)
      if(m.gt.1000) then
         write(istderr(),'(a,a)') 
     .     "delay: cannot handle embedding dimension "//
     .     "larger than than 1000"
         m=1000
      endif
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call readfile(nmax,x,nexcl,jcol,file,iverb)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_delay")
      call outfile(fout,iunit,iverb)
      do 10 n=(m-1)*id+1,nmax
 10      write(iunit,'(1000g16.7)') (x(n-(j-1)*id), j=m,1,-1)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-d# -m# -o outfile -l# -x# -c# -V# -h] file")
      call popt("d","delay (1)")
      call popt("m","embedding dimension (2)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_delay")
      call pall()
      stop
      end


