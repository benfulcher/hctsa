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
c   delay coordinates for periodic orbits
C   Copyright (C) Thomas Schreiber (1998)
c===========================================================================
      parameter(nx=1000)
      dimension x(nx)
      character*72 file, fout
      data m/2/
      data iverb/1/

      call whatido("embed using delay coordinates",iverb)
      id=imust("d")
      m=ican("m",m)
      ipp=ican("p",0)
      isout=igetout(fout,iverb)
      call nthstring(1,file)
      call infile(file,iunit,iverb)
      if(isout.eq.1) call addsuff(fout,file,"_delay")
      call outfile(fout,iunit2,iverb)

 1    read(iunit,*,err=1,end=999) ipor, dum1, dum2
      do 10 ip=1,ipor
 10      read(iunit,*,end=999) idum, x(ip)
      if(ipp.ne.0.and.ipor.ne.ipp) goto 1
      do 20 ip=1,ipor+1
 20      write(iunit2,*) (x(mod(ip-(j-1)*id-1+m*ipor,ipor)+1), j=m,1,-1)
      write(iunit2,'()')
      write(iunit2,'()')
      goto 1
 999  continue
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# [-m# -p# -o outfile -l# -x# -c# -V# -h] file")
      call popt("d","delay")
      call popt("m","embedding dimension (2)")
      call popt("p","period of orbit (1)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_delay")
      call pall()
      stop
      end


