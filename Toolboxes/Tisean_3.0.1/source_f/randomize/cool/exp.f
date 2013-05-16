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
c   part of the randomize-package for constraint surrogates
c   exponential cooling scheme
c   author T. Schreiber (1999)
c
c-------------------------------------------------------------------
c get options specific for cooling scheme
C
      subroutine opts_cool()
      common /coolcom/ 
     .   itini,tini,iafac,afac,cgoal,mtot,msucc,ntot,nsucc,mstop

      tini=fcan("T",0.)
      afac=fcan("a",0.)
      mtot=ican("S",20000)
      msucc=ican("s",2000)
      mstop=ican("z",200)
      cgoal=fcan("C",0.)
      end

c-------------------------------------------------------------------
c print version information on cooling scheme
C
      subroutine what_cool()
      call ptext("Cooling scheme: exponential")
      end

c-------------------------------------------------------------------
c print usage message specific for cooling scheme
C
      subroutine usage_cool()
      call ptext("Cooling options: [-T# -a# -S# -s# -z# -C#]")
      call popt("T","initial temperature (auto)")
      call popt("a","cooling factor (auto)")
      call popt("S","total steps before cooling (20000)")
      call popt("s","successful steps before cooling (2000)")
      call popt("z","minimal successful steps before cooling (200)")
      call popt("C","goal value of cost function (0.0)")
      end

c-------------------------------------------------------------------
c initialise all that is needed for cooling scheme
C
      function cool_init()
      common /coolcom/ 
     .   itini,tini,iafac,afac,cgoal,mtot,msucc,ntot,nsucc,mstop

      ntot=0
      nsucc=0
      itini=1
      if(tini.eq.0.) then
         tini=1e-4
         itini=0
      endif
      iafac=1
      if(afac.eq.0.) then
         afac=0.5
         iafac=0
      endif
      temp=tini
      cool_init=temp
      end
      
c-------------------------------------------------------------------
c determine new temperature depending on current cost function,
c acceptance status and history
c par can be used to pass information to the permutation scheme
c
      function cool(iaccept,iend,iv)
      common /coolcom/ 
     .   itini,tini,iafac,afac,cgoal,mtot,msucc,ntot,nsucc,mstop
      common nmax,cost,temp,cmin,rate

      iend=0
      cool=temp
      nsucc=nsucc+iaccept
      ntot=ntot+1
      if(ntot.lt.mtot.and.nsucc.lt.msucc) return
      rate=real(nsucc)/real(ntot)
      iend=1
      if(cost.le.cgoal) return
      if(itini.eq.0.and.temp.eq.tini.and.ntot.gt.1.5*nsucc) then
         tini=10*temp
         if(iv.ne.0) write(istderr(),*) 
     .      "increased initial temperature from ",
     .      temp, " to ", tini, " for melting"
         temp=tini
      else if(nsucc.le.mstop) then
         if(iafac.eq.1) return
         afac=sqrt(afac)
         mtot=mtot*sqrt(2.)
         temp=tini
         if(iv.ne.0) write(istderr(),*) "starting over: "
         if(iv.ne.0) write(istderr(),*) "   Cooling rate: ", afac, 
     .      " S:", mtot, " s: ", msucc
      else         
         temp=temp*afac
         if(iv.ne.0) write(istderr(),
     .      '(3hT: ,g15.6,4h S: ,i15,4h s: , i15,8h  cost: ,g15.6)') 
     .      temp, ntot, nsucc, cost
      endif
      iend=0
      ntot=0
      nsucc=0
      cool=temp
      end

