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
c   c2t.f
c   Takens' estimator from c2
c   author T. Schreiber (1998)
c
c   Modified for hctsa: all real variables are double precision, and the
c   points are sorted by an insertion sort in double precision.
c===========================================================================

      implicit double precision (a-h,o-z)
      parameter(meps=1000)
      dimension e(meps), c(meps)
      character*72 file, fout, aline
      data iverb/1/

      call whatido("Takens' estimator from correlation sum data",iverb)
      isout=igetout(fout,iverb)
      if(nstrings().eq.0) call usage()
      call nthstring(1,file)
      call infile(file,iunit,iverb)
      if(isout.eq.1) call addsuff(fout,file,"_t")
      call outfile(fout,iunit2,iverb)
 1    read(iunit,'(a)',end=999) aline
 4    if(aline(1:1).ne."#") goto 1
      if(aline(1:1).eq."#") 
     .   read(aline(index(aline,"m=")+2:72),'(i20)',err=1) m
      me=0
 2    read(iunit,'(a)') aline
      if(aline(1:72).eq." ") goto 3
      read(aline,*,err=999,end=999) ee, cc
      if(cc.le.0.) goto 3
      me=me+1
      e(me)=log(ee)
      c(me)=log(cc)
      goto 2
 3    write(iunit2,'("#m= ",i5)') m
c     sort the points by increasing e (insertion sort, in double precision)
      do 5 k=2,me
         et=e(k)
         ct=c(k)
         j=k-1
 6       if(j.lt.1) goto 7
         if(e(j).le.et) goto 7
         e(j+1)=e(j)
         c(j+1)=c(j)
         j=j-1
         goto 6
 7       e(j+1)=et
         c(j+1)=ct
 5    continue
      cint=0
      do 10 i=2,me
         b=(e(i)*c(i-1)-e(i-1)*c(i))/(e(i)-e(i-1))
         a=(c(i)-c(i-1))/(e(i)-e(i-1))
         if(a.ne.0) then
            cint=cint+(exp(b)/a)*(exp(a*e(i))-exp(a*e(i-1)))
         else
            cint=cint+exp(b)*(e(i)-e(i-1))
         endif
 10      write(iunit2,*) exp(e(i)), exp(c(i))/cint
      write(iunit2,'()') 
      write(iunit2,'()') 
      goto 4
 999  stop
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-o outfile -V# -h] file")
      call pout("file_t")
      call pall()
      stop
      end
