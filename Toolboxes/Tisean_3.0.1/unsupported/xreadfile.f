c multivariate i/o utilities
c Copyright (C) T. Schreiber (1999)

      subroutine xreadfile(nmax,mmax,nx,x,nexcl,icol,file,iverb)
c read columns as seperate time series
      parameter(mline=100)
      dimension x(nx,mmax), icol(mmax), dum(mline)
      character*(*) file
      character*10000 aline

      iv=iv_io(iverb)
      if(iv.ne.0) write(istderr(),*) 
     .   'reading from columns', (icol(i),i=1,mmax)
      call infile(file,iunit,iverb)
      mlast=0
      do 10 i=1,mmax
 10      mlast=max(mlast,icol(i))
      if(mlast.gt.mline) stop "xreadfile: make mline larger."
      lc=0
      do 20 n=1,nexcl
         lc=lc+1
 20      read(iunit,*,end=999)
      do 30 n=1,nmax
 1       lc=lc+1
         naline=(mlast+1)*100
         read(iunit,'(a)',err=2,end=999) aline(1:naline)
         if(aline(naline:naline).ne." ") then
            write(istderr(),*) "xreadfile: input line ", lc, " too long"
            stop
         endif
         read(aline,*,err=2,end=999) (dum(i),i=1,mlast)
         do 40 i=1,mmax
 40         x(n,i)=dum(icol(i))
         if(n.eq.1.and.iv_echo(iverb).eq.1) then
            write(istderr(),*) 
     .         "xreadfile: first data item read from line ", lc
            write(istderr(),*) (x(n,i), i=1,mmax)
         endif
         goto 30
 2       if(iv.ne.0) write(istderr(),*) "data in line ", lc, " ignored"
         goto 1
 30      continue
      if(iv.ne.0) write(istderr(),*) '*** readfile: warning:'//
     .   ' maybe not the whole file has been used'
 999  nmax=n-1
      if(iunit.ne.istdin()) close(iunit)
      if(iv.ne.0) call readreport(nmax,file)
      end

      subroutine xwritecfile(nmax,mmax,nx,x,file,iverb,comm)
c write comment and nmax points
      dimension x(nx,mmax)
      character*(*) file,comm

      if(mmax.gt.1000) then
         write(istderr(),*) "xwritecfile: "//
     .      "cannot write more than 1000 columns"
         stop
      endif
      call outfile(file,iunit,iverb)
      if(comm.ne." ") write(iunit,'(a)') comm
      do 10 n=1,nmax
 10      write(iunit,'(1000g16.7)') (x(n,i),i=1,mmax)
      if(iunit.eq.istdout()) then
         write(iunit,*)
         write(iunit,*)
      else
         close(iunit)
      endif
      if(iv_io(iverb).eq.1) call writereport(nmax,file)
      end

      subroutine xwritefile(nmax,mmax,nx,x,file,iverb)
c write nmax points
      dimension x(nx,mmax)
      character*(*) file

      call xwritecfile(nmax,mmax,nx,x,file,iverb," ")
      end

      subroutine columns(mc,mmax,icol)
      dimension icol(*)

      call imcan("c",mmax,mc,icol)
      icmax=0
      do 10 m=1,mc
 10      icmax=max(icmax,icol(m))
      do 20 m=mc+1,mmax
         icmax=icmax+1
 20      icol(m)=icmax
      end



