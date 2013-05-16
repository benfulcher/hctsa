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
c   multivariate i/o utilities for TISEAN f-sources
c   author T. Schreiber (1999)
c===========================================================================
      subroutine xreadfile(nmax,mmax,nx,x,nexcl,icol,file,iverb)
c read columns as seperate time series
      parameter(mline=1000)
      dimension x(nx,mmax), icol(mmax), dum(mline)
      character*(*) file

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
         read(iunit,*,err=2,end=999)  (dum(i),i=1,mlast)
         do 40 i=1,mmax
 40         x(n,i)=dum(icol(i))
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



