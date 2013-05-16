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
c   clustering a dissimilarity matrix
c   see Schreiber and Schmitz, Phys. Rev. Lett. 79 (1997) 1475
c   author T. Schreiber (1998)
c===========================================================================

      parameter(npmax=1000)
      dimension d(npmax,npmax), iu(npmax), ifix(npmax)
      character*72 file, fout, filex
      data iverb/3/

      call whatido("clustering a dissimilarity matrix",iverb)
      ncl=imust("#")
      iflag=lopt("=",1)
      call stcan('X',filex,' ')
      isout=igetout(fout,iverb)

      call nthstring(1,file)
      call infile(file,iunit,iverb)
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_clust")
      do 10 i=1,npmax
         do 10 j=1,npmax
 10         d(i,j)=-1e20
      np=0
 1    read(iunit,*,end=999) i,j,dij
      d(i,j)=dij
      np=max(i,j,np)
      goto 1
 999  if(iv_io(iverb).eq.1) write(0,'(a,i)') "matrix size ", np
      dmean=0
      nd=0
      do 20 i=1,np
         do 20 j=1,np
            if(d(i,j).ne.-1e20) then
               nd=nd+1
               dmean=dmean+d(i,j)
            endif
 20         continue
      do 30 i=1,np
         do 30 j=1,np
 30         if(d(i,j).eq.-1e20) d(i,j)=dmean/nd
      do 40 i=1,np
 40      ifix(i)=0
      if(filex.ne." ") then
         open(10,file=filex,status='old',err=998)
         nfix=0
 2       read(10,*,end=998,err=2) i, iff
         if(i.lt.1.or.i.gt.np.or.iff.gt.ncl.or.iff.lt.1) goto 1
         ifix(i)=iff
         nfix=nfix+1
      endif
 998  if(nfix.eq.np) stop "all fixed."
      call clustering(np,d,npmax,ncl,nfix,ifix,iu,iverb,iflag)
      call outfile(fout,iunit,iverb)
      do 50 n=1,np
 50      write(iunit,*) iu(n), (costi(np,iu,d,n,ic,iflag),ic=1,ncl)
      end

      subroutine usage()
c usage message

      call whatineed("-## [-= -X xfile] file")
      call popt("#","number of clusters")
      call popt("=","if set, bias towards similar size clusters")
      call popt("X","list of indices with fixed cluster assignments")
      call pout("file_clust")
      call pall()
      call ptext("Verbosity levels (add what you want):")
      call ptext("          1 = input/output" )
      call ptext("          2 = state of clustering")
      call ptext("          8 = temperature / cost at cooling")
      stop
      end

      subroutine clustering(np,d,npmax,ncl,nfix,ifix,iu,iverb,iflag)
      parameter(nt0=20,tfac=10.,tstep=0.99,ntotmaxf=20,nsuccmaxf=2)
      external rand
      character*1 c
      dimension d(npmax,npmax), iu(*), ifix(*)
      equivalence (c,ic)
      data c/'A'/

      ntotmax=(np-nfix)*ntotmaxf
      nsuccmax=(np-nfix)*nsuccmaxf
      se=0.
      se2=0.
      do 10 nt=1,nt0
         call ranconf(np,iu,ncl,ifix)
         e=cost(np,iu,d,ncl,iflag)
         se=se+e
 10      se2=se2+e**2
      t=tfac*sqrt(se2/nt0-(se/nt0)**2)

      ntot=0
      nsucc=0
 1    call cconf(np,iu,ncl,nch,iuold,ifix)
      ec=cost(np,iu,d,ncl,iflag)
      ntot=ntot+1
      if(ec.lt.e.or.(rand(0.0).lt.exp(-(ec-e)/t))) then
         e=ec
         nsucc=nsucc+1
      else
         iu(nch)=iuold
      endif
      if(ntot.eq.ntotmax .or. nsucc.eq.nsuccmax) then
         if(nsucc.eq.0) return
         ntot=0
         nsucc=0
         if(iv_clust(iverb).eq.1) write(istderr(),'(80a1)') 
     .      (ic+iu(n)-1,n=1,np)
         if(iv_cool(iverb).eq.1) write(istderr(),*) t, e
         t=t*tstep
      endif
      goto 1
      end

      function cost(np,iu,d,ncl,iflag)
      parameter(npmax=1000)
      dimension d(npmax,npmax), iu(*), ictab(npmax)
      
      cost=0
      do 10 ic=1,ncl
         nic=0
         do 20 n=1,np
            if(iu(n).ne.ic) goto 20
            nic=nic+1
            ictab(nic)=n
 20         continue
         cc=0
         do 30 ii=1,nic
            i=ictab(ii)
            do 30 jj=1,nic
               j=ictab(jj)
 30            cc=cc+d(i,j)
 10      if(nic.gt.0) cost=cost+cc/(1+(1-iflag)*(nic-1))
      end

      function costi(np,iu,d,nn,ic,iflag)
      parameter(npmax=1000)
      dimension d(npmax,npmax), iu(*), ictab(npmax)
      
      costi=0
      nic=0
      do 20 n=1,np
         if(iu(n).ne.ic) goto 20
         nic=nic+1
         ictab(nic)=n
 20      continue
      cc=0
      do 30 jj=1,nic
         j=ictab(jj)
 30      cc=cc+d(nn,j)+d(j,nn)
      if(nic.gt.0) costi=0.5*cc/(1+(1-iflag)*(nic-1))
      end

      subroutine ranconf(np,iu,ncl,ifix)
      external rand
      dimension iu(*), ifix(*)

      do 10 n=1,np
         iu(n)=ifix(n)
 10      if(ifix(n).eq.0) iu(n)=min(int(rand(0.0)*ncl)+1,ncl)
      end

      subroutine cconf(np,iu,ncl,nch,iuold,ifix)
      external rand
      dimension iu(*), ifix(*)

 1    nch=min(int(rand(0.0)*np)+1,np)
      if(ifix(nch).ne.0) goto 1
      iuold=iu(nch)
      iu(nch)=iuold+int(rand(0.0)*(ncl-1))+1
      if(iu(nch).gt.ncl) iu(nch)=iu(nch)-ncl
      end
