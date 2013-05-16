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
c   endtoend.f
c   Determine end-to-end mismatch before making surrogate data
c   author T. Schreiber (1999)
c===========================================================================

      parameter(nx=100000,mx=20)
      dimension x(nx,mx), icol(mx)
      character*72 file, fout
      data iverb/15/

      call whatido("Determine end-to-end mismatch",iverb)
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      wjump=fcan("j",0.5)
      mcmax=ican("m",0)
      call columns(mc,mx,icol)
      if(mcmax.eq.0) mcmax=max(1,mc)
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call xreadfile(nmax,mcmax,nx,x,nexcl,icol,file,iverb)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_end")
      call outfile(fout,iunit,iverb)
      nmaxp=nmax
      etotm=mcmax
 1    nmaxp=nless(nmaxp)
      call jump(nmax,nmaxp,nx,x,mcmax,wjump,etot,ejump,eslip,njump)
      if(etot.lt.etotm) then
         etotm=etot
         write(iunit,'(a,i7,a,i7,a,f5.1,a)')
     .      "length:", nmaxp, 
     .      "  offset: ", nexcl+njump,
     .      "  lost: ", real(nmax-nmaxp)/real(nmax)*100, " %"
         write(iunit,*) "      jump: ", ejump*100, " %"
         write(iunit,*) "      slip: ", eslip*100, " %" 
         write(iunit,*) "  weighted: ", etot*100, " %"
         write(iunit,'()')
      endif
      if(etot.lt.1e-5) stop
      nmaxp=nmaxp-1
      if(nmaxp.gt.2) goto 1
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-j# -o outfile -l# -x# -c# -V# -h] file")
      call popt("j","weight given to jump relative to slip (0.5)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("m","number of columns to be read (1)")
      call popt("c","columns to be read (1)")
      call pout("file_end")
      call pall()
      stop
      end

      subroutine jump(nmax,nmaxp,nx,x,mcmax,wjump,etot,ejump,eslip,
     .   njump)
c loop through time ofsets to minimize jump effect
      dimension x(nx,*)

      etot=mcmax
      do 10 nj=0,nmax-nmaxp
         xj=0
         sj=0
         do 20 m=1,mcmax
            xj=xj+xjump(nmaxp,x(1+nj,m))
 20         sj=sj+sjump(nmaxp,x(1+nj,m))
         if(wjump*xj+(1-wjump)*sj.ge.etot) goto 10
         etot=wjump*xj+(1-wjump)*sj
         ejump=xj
         eslip=sj
         njump=nj
 10      continue
      end

      function xjump(nmax,x)
c contribution of end effect to 1st derivative
      dimension x(*)

      call rms(nmax,x,sc,sd)
      xjump=0
      if(sd.eq.0.) return
      xjump=(x(1)-x(nmax))**2/(nmax*sd**2)
      end

      function sjump(nmax,x)
c contribution of end effect to 2nd derivative
      dimension x(*)

      call rms(nmax,x,sc,sd)
      sjump=0
      if(sd.eq.0.) return
      sjump=((x(nmax)-x(nmax-1))-(x(2)-x(1)))**2 / (nmax*sd**2)
      end

