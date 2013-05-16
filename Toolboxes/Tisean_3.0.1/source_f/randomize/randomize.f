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
c   constrained randomization
c   author T. Schreiber (1999)
c===========================================================================
      parameter(nx=100000,mx=20) 
      double precision time
      dimension x(nx,mx), y(nx,mx), xx(nx,mx), icol(mx)
      character*72 file, fout, comment
      common nmax,cost,temp,cmin,rate,x
      external rand
      data wr/0.9/, nsur/1/
      data iverb/15/

      call whatido("constrained randomization",iverb)
      call what_cost()
      call what_cool()
      call what_permute()
      rr=rand(ican("I",0)/real(2**22))
      nsur=min(999,ican("n",nsur))
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      mcmax=ican("m",0)
      call columns(mc,mx,icol)
      if(mcmax.eq.0) mcmax=max(1,mc)
      wr=fcan("u",wr)
      call opts_cost(mcmax)
      call opts_cool()
      call opts_permute()
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call xreadfile(nmax,mcmax,nx,xx,nexcl,icol,file,iverb)
      call cost_transform(nmax,mcmax,nx,xx)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_rnd")
      if(nsur.gt.1) call suffix(fout,"_000")

      do 10 isur=1,nsur
         do 20 n=1,nmax
            do 20 m=1,mcmax
 20            x(n,m)=xx(n,m)
         rate=1
         cost=r1mach(2)
         cmin=cost
         if(nsur.gt.1) write(fout(index(fout," ")-3:72),'(i3.3)') isur
         temp=cool_init()
         call cost_init()
         call permute_init()
         call cost_full(iv_vcost(iverb))
         cmin=cost
         time=0
 1       time=time+1.
         call permute(n1,n2)
         cmax=cost-temp*log(rand(0.0))    ! maximal acceptable cost
         call cost_update(n1,n2,cmax,iaccept,iv_vmatch(iverb))
         tnew=cool(iaccept,iend,iv_cool(iverb))
         if(tnew.ne.temp.or.cost.lt.cmin*wr) then
            cc=cost
            call cost_full(iv_vcost(iverb))
            if(iv_match(iverb).eq.1) write(istderr(),*) 
     .         "cost function mismatch: ", abs((cc-cost)/cost)
         endif
         temp=tnew
         if(cost.lt.cmin*wr) then
            if(iv_cost(iverb).eq.1) write(istderr(),*) 
     .         "after ",real(time)," steps at T=",temp," cost: ",cost
            cmin=cost
            call cost_inverse(nmax,mcmax,nx,x,y)
            write(comment,'(8h# cost: ,g15.5)') cost
            call xwritecfile(nmax,mcmax,nx,y,fout,iverb,comment)
         endif
         if(iend.ne.1) goto 1
         write(comment,'(8h# cost: ,g15.5)') cost
         call writecfile(nmax,mcmax,nx,y,fout,iverb,comment)
 10      continue
      end

      subroutine usage()
c usage message

      call whatineed("[-n# -u# -I# -o outfile -l# -x# -c# -V# -h]"//
     .   " [cost opt.] [cooling opt.] [permutation opt.] file")
      call popt("n","number of surrogates (1)")
      call popt("u","improvement factor before write (0.9)")
      call popt("I","seed for random numbers (0)")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column to be read (1 or file,#)")
      call pout("file_rnd(_nnn)")
      call pall()
      call ptext("Verbosity levels (add what you want):")
      call ptext("          1 = input/output" )
      call ptext("          2 = current cost if improved")
      call ptext("          4 = cost mismatch")
      call ptext("          8 = temperature etc. at cooling")
      call ptext("         16 = verbose cost if improved")
      call ptext("         32 = verbose cost mismatch")
      write(istderr(),'()') 
      call usage_cost()
      write(istderr(),'()') 
      call usage_cool()
      write(istderr(),'()') 
      call usage_permute()
      write(istderr(),'()') 
      stop
      end

