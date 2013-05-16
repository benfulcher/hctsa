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
c   embed using principal components
C   author Thomas Schreiber (1998)
c===========================================================================
      parameter(nx=1000000, me=500)
      dimension x(nx), c(me,me), d(me), xc(me), z(me,me)
      character*72 file, fout
      data id/1/, isvd/2/
      data iverb/1/

      call whatido("embed using principal components",iverb)
      m=imust("m")
      if(m.gt.me) stop "svd: make me larger."
      id=ican("d",id)
      isvd=min(ican("q",isvd),m)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      jcol=ican("c",0)
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call readfile(nmax,x,nexcl,jcol,file,iverb)
      call normal(nmax,x,sc,sd)
      call svd_vectors(nmax,m,id,x,c,z,d)
      if(iv_io(iverb).eq.1) write(istderr(),*) 
     .   "#, fraction of variance, accumulative fraction"
      ctot=0.
      do 10 i=1,m
 10      ctot=ctot+d(m+1-i)
      cacc=0.
      do 20 i=1,m
         cacc=cacc+d(m+1-i)
 20      if(iv_io(iverb).eq.1) 
     .      write(istderr(),*) i, d(m+1-i)/ctot, cacc/ctot
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_pc")
      call outfile(fout,iunit,iverb)
      do 30 n=(m-1)*id+1,nmax
         do 40 i=1,isvd
            s=0
            do 50 j=1,m
 50            s=s+z(j,m+1-i)*x(n-(j-1)*id)
 40         xc(i)=s
 30      write(iunit,*) (xc(i),i=1,isvd)
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-m# [-d# -q# -o outfile -l# -x# -c# -V# -h] file")
      call popt("m","initial embedding dimension")
      call popt("d","delay for initial embedding (1)")
      call popt("q","number of principal components (2)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_pc")
      call pall()
      stop
      end

      subroutine svd_vectors(nmax,m,id,x,c,z,d)
      parameter(me=500)
      dimension x(nmax), c(me,*), d(m), w1(me), w2(me), z(me,*)
  
      if(m.gt.me) stop "svd_vectors: make me larger."
      do 10 i=1,m
         do 10 j=i,m
            s=0.
            do 20 n=(m-1)*id+1,nmax
 20            s=s+x(n-(i-1)*id)*x(n-(j-1)*id)
            c(i,j)=s/(nmax-(m-1)*id)
 10         c(j,i)=c(i,j)
      call rs(me,m,c,d,1,z,w1,w2,ierr)
      end
