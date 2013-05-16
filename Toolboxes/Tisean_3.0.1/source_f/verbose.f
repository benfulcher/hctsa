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
c   utilities for TISEAN f-sources
c===========================================================================
      function igetv(idef)
c get verbosity level

      igetv=ican("V",-1)
      if(igetv.eq.-1.and.lopt("V",1).eq.1) igetv=2**15-1
      if(igetv.eq.-1) igetv=idef
      end

      function iexv(iverb,item)
c 1 if verbosity level includes item

      iexv=iand(iverb,item)/item
      end

c the following functions test for specific numerical verbosity values

      function iv_io(iverb)           ! report i/o activity
      iv_io=iexv(iverb,1)
      end

      function iv_echo(iverb)         ! echo first line of data read
      iv_echo=iexv(iverb,128)
      end

      function iv_cost(iverb)         ! current value of cost function 
      iv_cost=iexv(iverb,2)
      end

      function iv_match(iverb)        ! cost mismatch 
      iv_match=iexv(iverb,4)
      end

      function iv_cool(iverb)         ! temperature etc. at cooling 
      iv_cool=iexv(iverb,8)
      end

      function iv_vcost(iverb)        ! verbose cost if improved 
      iv_vcost=iexv(iverb,16)
      end

      function iv_vmatch(iverb)       ! verbose cost mismatch
      iv_vmatch=iexv(iverb,32)
      end

      function iv_10(iverb)           ! upo status after 10 points
      iv_10=iexv(iverb,16)
      end

      function iv_100(iverb)          ! upo status after 100 points
      iv_100=iexv(iverb,8)
      end

      function iv_1000(iverb)         ! upo status after 1000 points
      iv_1000=iexv(iverb,4)
      end

      function iv_upo(iverb)          ! print orbits found
      iv_upo=iexv(iverb,2)
      end

      function iv_surr(iverb)         ! print iterations / discrepancy
      iv_surr=iexv(iverb,2)           
      end

      function iv_uncorr(iverb)       ! neighbour search status
      iv_uncorr=iexv(iverb,2)           
      end

      function iv_clust(iverb)        ! clustering status
      iv_clust=iexv(iverb,2)           
      end




