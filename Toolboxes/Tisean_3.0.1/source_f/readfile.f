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
c   i/o utilities for TISEAN f-sources
c   author T. Schreiber (1998) based on earlier versions
c===========================================================================
      subroutine readfile(nmax,x,nexcl,icol,file,iverb)
c read at most nmax points, return nmax
      dimension x(nmax)
      character*(*) file

      iv=iv_io(iverb)
      if(icol.eq.0) icol=igetcol(file)
      if(icol.gt.0.and.iv.ne.0) 
     .   write(istderr(),*) 'reading from column', icol
      call infile(file,iunit,iverb)
      lc=0
      do 10 n=1,nexcl
         lc=lc+1
 10      read(iunit,*,end=999)
      do 20 n=1,nmax
 1       lc=lc+1
         read(iunit,*,err=2,end=999) (dum,i=1,icol-1), x(n)
         goto 20
 2       if(iv.ne.0) write(istderr(),*) "data in line ", lc, " ignored"
         goto 1
 20      continue
      if(iv.ne.0) write(istderr(),*) '*** readfile: warning:'//
     .   ' maybe not the whole file has been used'
 999  nmax=n-1
      if(iunit.ne.istdin()) close(iunit)
      if(iv.ne.0) call readreport(nmax,file)
      if(icol.gt.0.and.file.ne."-") call putcol(file,icol)
      end

      function igetcol(file)
      character*(*) file

      igetcol=0
      do 10 i=len(file),1,-1
 10      if(file(i:i).eq.",") goto 1
 1    if(i.eq.0) return
      read(file(i+1:len(file)),'(i10)',err=999) igetcol
      file(i:len(file))=" "
 999  continue
      end

      subroutine putcol(file,icol)
      character*(*) file

      if(icol.le.9) then
         write(file(index(file," "):index(file," ")+1),'(1h,,i1)') icol
      else 
         write(file(index(file," "):index(file," ")+2),'(1h,,i2)') icol
      endif
      end

      subroutine writecfile(nmax,x,file,iverb,comm)
c write comment and nmax points
      dimension x(nmax)
      character*(*) file,comm

      call outfile(file,iunit,iverb)
      if(comm.ne." ") write(iunit,'(a)') comm
      do 10 n=1,nmax
 10      write(iunit,*) x(n)
      if(iunit.eq.istdout()) then
         write(iunit,*)
         write(iunit,*)
      else
         close(iunit)
      endif
      if(iv_io(iverb).eq.1) call writereport(nmax,file)
      end

      subroutine writefile(nmax,x,file,iverb)
c write nmax points
      dimension x(nmax)
      character*(*) file

      call writecfile(nmax,x,file,iverb," ")
      end

      subroutine infile(file,iunit,iverb)
c open file for read on iunit=ifile(), or iunit=istdin() if "-"      
      character*(*) file

      if(file.eq."-") then
         iunit=istdin()
         if(iv_io(iverb).eq.1) write(istderr(),*) "reading from stdin"
         return
      endif
      iunit=ifilein()
      open(iunit,file=file,status="old",err=999)
      if(iv_io(iverb).eq.1) write(istderr(),'(a,a,a)') 
     .   "opened ",file(1:index(file," ")-1), " for input"
      return
 999  write(istderr(),'(a,a)') "Cannot open input file ",
     .   file(1:index(file," ")-1)
      stop
      end

      subroutine outfile(file,iunit,iverb)
c open file for write on iunit=ifileout(), or iunit=istdout() if file=" "      
      character*(*) file

      if(file.eq." ") then
         iunit=istdout()
         if(iv_io(iverb).eq.1) write(istderr(),*) "writing to stdout"
         return
      endif
      iunit=ifileout()
      open(iunit,file=file,status='unknown',err=999)
      if(iv_io(iverb).eq.1) write(istderr(),'(a,a,a)') 
     .   "opened ",file(1:index(file," ")-1), " for output"
      return
 999  write(istderr(),'(a,a)') "Cannot open output file ",
     .   file(1:index(file," ")-1)
      stop
      end

      subroutine suffix(base,suff)
c append stuff after last nonblank character in base
      character*(*) base, suff

      base=base(1:index(base," ")-1)//suff
      end

      subroutine addsuff(target,base,suff)
c append stuff after last nonblank character in base
      character*(*) target,base, suff

      target=base(1:index(base," ")-1)//suff
      end

      subroutine readreport(nmax,file)
c report on numbers read
      character*(*) file

      if(file.eq."-") then
         write(istderr(),'(i10,a)') nmax, ' values read from stdin'
      else
         write(istderr(),'(i10,a,a)') nmax, ' values read from file: ', 
     .      file(1:index(file," ")-1)
      endif
      if(nmax.ne.0) return
      if(file.eq."-") then
         write(istderr(),'(a)') "No input given - aborting."
      else
         write(istderr(),'(a,a,a)') "Input file ",
     .      file(1:index(file," ")-1), " empty - aborting."
      endif
      call usage()
      end

      subroutine writereport(nmax,file)
c report on numbers written
      character*(*) file

      if(file.eq." ") then
         write(istderr(),'(i10,a)') nmax, ' values written to stdout'
      else
         write(istderr(),'(i10,a,a)') nmax, ' values written to file: ', 
     .      file(1:index(file," ")-1)
      endif
      end

