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
c   commandline.f
c   get command line options
c   author T. Schreiber (1998)
c===========================================================================

      function iopt(c,ith,ierr)
c get ith occurence of switch -c as integer
      character*72 argv
      character c

      iopt=0
      ifound=0
      do 10 i=1,iargc()
         call getarg(i,argv)
         if(argv(1:2).eq.'-'//c) then
            ifound=ifound+1
            if(ifound.eq.ith) then
               call argdel(i)
               if(argv(3:72).ne.' ') then
                  iopt=i_s(argv(3:72),ierr)
               else if(i+1.le.iargc()) then
                  call getarg(i+1,argv)
                  iopt=i_s(argv,ierr)
                  if(ierr.eq.0) call argdel(i+1)
               else
                  ierr=1
               endif
               return
            endif
         endif
 10      continue
      ierr=1
      end

      function fopt(c,ith,ierr)
c get ith occurence of switch -c as real
      character*72 argv
      character c

      fopt=0
      ifound=0
      do 10 i=1,iargc()
         call getarg(i,argv)
         if(argv(1:2).eq.'-'//c) then
            ifound=ifound+1
            if(ifound.eq.ith) then
               call argdel(i)
               if(argv(3:72).ne.' ') then
                  fopt=f_s(argv(3:72),ierr)
               else if(i+1.le.iargc()) then
                  call getarg(i+1,argv)
                  fopt=f_s(argv,ierr)
                  if(ierr.eq.0) call argdel(i+1)
               else
                  ierr=1
               endif
               return
            endif
         endif
 10      continue
      ierr=1
      end

      subroutine sopt(c,ith,string,ierr)
c get ith occurence of switch -c as string
      character*(*) string
      character c

      ifound=0
      do 10 i=1,iargc()
         call getarg(i,string)
         if(string(1:2).eq.'-'//c) then
            ifound=ifound+1
            if(ifound.eq.ith) then
               call argdel(i)
               if(string(3:).ne.' ') then
                  string=string(3:)
                  ierr=0
               else if(i+1.le.iargc()) then
                  call getarg(i+1,string)
                  if(string(1:1).eq."-") then
                     ierr=1
                     return
                  endif
                  call argdel(i+1)
                  ierr=0
               else
                  ierr=1
               endif
               return
            endif
         endif
 10      continue
      ierr=1
      end

      function lopt(c,ith)
c test if ith occurence of switch -c is present
      character*72 argv
      character c

      lopt=0
      ifound=0
      do 10 i=1,iargc()
         call getarg(i,argv)
         if(argv(1:2).eq.'-'//c) then
            ifound=ifound+1
            if(ifound.eq.ith) then
               lopt=1
               call argdel(i)
               return
            endif
         endif
 10      continue
      end

      function iget(inum)
c get inum'th argument as integer
      character*72 argv
      
      iget=0
      call getarg(inum,argv)
      if(argv.eq.' ') 
     .write(istderr(),'(a,i10)') "iget: missing integer argument",inum
      iget=i_s(argv,ierr)
      if(ierr.ne.0) 
     .write(istderr(),'(a,i10)') "iget: integer argument expected:",inum
      end

      function fget(inum)
c get inum'th argument as real
      character*72 argv
      
      fget=0
      call getarg(inum,argv)
      if(argv.eq.' ') 
     .   write(istderr(),'(a)') "fget: missing real argument",inum
      fget=f_s(argv,ierr)
      if(ierr.ne.0) 
     .   write(istderr(),'(a)') "fget: real argument expected:;",inum
      end
