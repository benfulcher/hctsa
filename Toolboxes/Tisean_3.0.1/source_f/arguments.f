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
c   arguments.f
c   get command line arguments 
c   author T. Schreiber (1998)
c===========================================================================

      subroutine argdel(i)
      parameter(margs=1000)
      dimension largs(margs)
      common /args/ nargs, largs

      if(i.eq.0) then
         nargs=min(margs,iargc())
         do 10 n=1,nargs
 10         largs(n)=1
      else
         if(i.gt.iargc()) return
         if(largs(i).eq.0) return
         largs(i)=0
         nargs=nargs-1
      endif
      end

      function nstrings()
      parameter(margs=1000)
      dimension largs(margs)
      common /args/ nargs, largs

      nstrings=max(nargs,1)
      end

      subroutine nthstring(n,string)
      parameter(margs=1000)
      dimension largs(margs)
      common /args/ nargs, largs
      character*(*) string

      iv=0
      do 10 i=1,iargc()
         if(largs(i).eq.1) iv=iv+1
 10      if(iv.eq.n) goto 1
      string="-"
      return
 1    call getarg(i,string)
      end

      function imust(c)
c get mandatory integer argument, call usage statement if missing
      character c

      imust=iopt(c,1,ierr)
      if(ierr.ne.0) call usage()
      end

      function fmust(c)
c get mandatory real argument, call usage statement if missing
      character c

      fmust=fopt(c,1,ierr)
      if(ierr.ne.0) call usage()
      end

      subroutine smust(c,string)
c get mandatory string argument, call usage statement if missing
      character c
      character*(*) string

      call sopt(c,1,string,ierr)
      if(ierr.ne.0) call usage()
      end

      function ican(c,idef)
c get optional integer argument, provide default if missing
      character c

      ican=iopt(c,1,ierr)
      if(ierr.ne.0) ican=idef
      end
      
      function fcan(c,fdef)
c get optional real argument, provide default if missing
      character c

      fcan=fopt(c,1,ierr)
      if(ierr.ne.0) fcan=fdef
      end

      subroutine stcan(c,string,dstring)
c get optional string argument, provide default if missing
      character c
      character*(*) string, dstring

      call sopt(c,1,string,ierr)
      if(ierr.ne.0) string=dstring
      end

      function igetout(fout,iverb)
c gets alternate output file name, default " "
c return 1 if fout must be determined from input file name
      character*(*) fout

      igetout=0
      call stcan("o",fout," ")
      if(fout.ne." ".and.nstrings().gt.1.and.iv_io(iverb).ne.0) 
     .   write(istderr(),*) '*** single output file for multiple'//
     .   ' input files - results may be overwritten'
      if(fout.ne." ") return
      igetout=lopt("o",1)
      end

      subroutine imcan(c,mmax,mc,ilist)
c get optional integer argument with multiple comma separated values
      character c
      character*72 string
      dimension ilist(*)

      call stcan(c,string," ")
      string(index(string," "):index(string," "))=","
      do 10 m=1,mmax
         if(index(string,",").le.1) goto 1
         read(string(1:index(string,",")-1),*,err=1,end=1) ilist(m)
 10      string=string(index(string,",")+1:72)
 1    mc=m-1
      end

      subroutine fmcan(c,mmax,mc,flist)
c get optional real argument with multiple comma separated values
      character c
      character*72 string
      dimension flist(*)

      call stcan(c,string," ")
      string(index(string," "):index(string," "))=","
      do 10 m=1,mmax
         if(index(string,",").le.1) goto 1
         read(string(1:index(string,",")-1),*,err=1,end=1) flist(m)
 10      string=string(index(string,",")+1:72)
 1    mc=m-1
      end
