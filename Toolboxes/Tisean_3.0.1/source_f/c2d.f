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
c   c2d.f
c   local slopes from c2
c   author T. Schreiber (1998)
c===========================================================================

      parameter(meps=1000)
      dimension e(meps), c(meps)
      character*72 file, fout, aline
      data iav/1/
      data iverb/1/

      call whatido("local slopes from c1/c2 correlation sum data",iverb)
      iav=ican('a',iav)
      isout=igetout(fout,iverb)
      if(nstrings().eq.0) call usage()
      call nthstring(1,file)
      call infile(file,iunit,iverb)
      if(isout.eq.1) call addsuff(fout,file,"_d")
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
 3    write(iunit2,'(4h#m= ,i5)') m
      do 30 j=iav+1,me-iav
         call slope(e(j-iav),c(j-iav),2*iav+1,s)
 30      if(s.gt.0.) write(iunit2,*) exp(0.5*(e(j+iav)+e(j-iav))),  s
      write(iunit2,'()') 
      write(iunit2,'()') 
      goto 4
 999  stop
      end

      subroutine slope(x,y,n,a)
      dimension x(n),y(n)

      sx=0.
      sa=0
      a=0.
      do 10 i=1,n
 10      sx=sx+x(i)
      do 20 i=1,n
         sa=sa+(x(i)-sx/n)**2
 20      a=a+y(i)*(x(i)-sx/n)
      a=a/sa
      end


      subroutine usage()
c usage message

      call whatineed(
     .   "[-a# -o outfile -V# -h] file")
      call popt("a","average using -#,...,+# [1]")
      call pout("file_d")
      call pall()
      stop
      end


