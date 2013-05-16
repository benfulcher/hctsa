c template cost function
c Copyright (C) T. Schreiber (1999)

c-------------------------------------------------------------------
c get cost function specific options
c
      subroutine opts_cost()
      end

c-------------------------------------------------------------------
c print version information on cost function
c
      subroutine what_cost()
      call ptext("Cost function: generic, refer to source")
      end

c-------------------------------------------------------------------
c print cost function specific usage message
c
      subroutine usage_cost()
      call ptext("Cost function options: (none)")
      end

c-------------------------------------------------------------------
c initialise all that is needed for cost function
c
      subroutine cost_init()
      parameter(mitem=100)
      dimension c0(mitem), c(mitem)
      common /costcom/ c0, c, sd, sc

      do 10 n=1,mitem
 10      if(igeneric(n,c0).eq.0) return
      end

c-------------------------------------------------------------------
c initial transformation on time series and its inverse
c
      subroutine cost_transform(nmax,x)
      dimension x(nmax)
      parameter(mitem=100)
      dimension c0(mitem), c(mitem)
      common /costcom/ c0, c, sd, sc

      call normal1(nmax,x,sc,sd)
      end

      subroutine cost_inverse(nmax,x,y)
      dimension x(nmax), y(nmax)
      parameter(mitem=100)
      dimension c0(mitem), c(mitem)
      common /costcom/ c0, c, sd, sc
      
      do 10 n=1,nmax
 10      y(n)=x(n)*sd+sc
      end

c-------------------------------------------------------------------
c compute full cost function from scratch
c
      subroutine cost_full(iv)
      parameter(mitem=100)
      dimension c0(mitem), c(mitem)
      common /costcom/ c0, c, sd, sc
      common nmax,cost

      cost=0  
      if(iv.ne.0) write(istderr(),*) "goal / is  / miss"
      do 10 n=1,mitem
         if(igeneric(n,c).eq.0) return
         if(iv.ne.0) write(istderr(),*) c0(n), c(n), c0(n)-c(n)
 10      call aver(cost,c(n)-c0(n))
      end

c-------------------------------------------------------------------
c compute changed cost function on exchange of n1 and n2 
c
      subroutine cost_update(n1,n2,cmax,iaccept,iv)
      parameter(mitem=100,nx=100000)
      dimension c0(mitem), c(mitem), ccop(mitem), x(nx)
      common /costcom/ c0, c, sd, sc
      common nmax,cost,temp,cmin,rate,x

      comp=0
      iaccept=0
      do 10 n=1,mitem
         if(iugeneric(n,c,ccop,n1,n2).eq.0) goto 1
         call aver(comp,ccop(n)-c0(n))
 10      if(comp.ge.cmax) return
 1    call exch(n1,n2)      ! if got here: accept
      cost=comp  
      iaccept=1
      if(iv.ne.0) call panic(ccop)
      do 20 i=1,n-1
 20      c(i)=ccop(i)
      end

      subroutine panic(ccop)
      parameter(mitem=100)
      dimension c0(mitem), c(mitem), ccop(*)
      common /costcom/ c0, c, sd, sc

      write(istderr(),*) "update / fresh / discrepancy"
      do 10 n=1,mitem
         if(igeneric(n,c).eq.0) return
 10      write(istderr(),*) ccop(n), c(n), ccop(n)-c(n)
      end

c-------------------------------------------------------------------
c generic list of statistics
c here the number and order of quantities has to be coded at compile time
c add or delete items as you please
c igeneric 
c     will compute them one at a time and return value 0 when it is done
c iugeneric 
c     will cdo the same as an update when n1 and n2 are changed
c     it is easiest but inefficient just to call igeneric
c aver 
c     will do the averaging as desired
c
      function igeneric(n,c)
      parameter(nx=100000)
      dimension c(*), x(nx)
      common nmax,cost,temp,cmin,rate,x

      igeneric=1
      cc=0
      if(n.le.2) then      ! lag 1 and 2 autocorrelations
         do 10 i=n+1,nmax
 10         cc=cc+x(i-n)*x(i)
      else if(n.eq.3) then ! lag 1,1  3-point-autocorrelation
         do 20 i=2,nmax
 20         cc=cc+x(i-1)**2*x(i)  
      else if(n.eq.4) then ! lag 1,0  3-point-autocorrelation
         do 30 i=2,nmax
 30         cc=cc+x(i-1)*x(i)**2
      else if(n.eq.5) then ! lag 2,0  3-point-autocorrelation
         do 40 i=3,nmax
 40         cc=cc+x(i-2)**2*x(i)
      else if(n.eq.6) then ! lag 2,1  3-point-autocorrelation
         do 50 i=3,nmax
 50         cc=cc+x(i-2)*x(i-1)*x(i)
      else if(n.eq.7) then ! lag 1,1,0  4-point-autocorrelation
         do 60 i=2,nmax
 60         cc=cc+x(i-1)**2*x(i)**2
      else if(n.eq.8) then ! lag 1,0,0  4-point-autocorrelation
         do 70 i=2,nmax
 70         cc=cc+x(i-1)*x(i)**3
      else if(n.eq.9) then ! lag 1,1,1  4-point-autocorrelation
         do 80 i=2,nmax
 80         cc=cc+x(i-1)**3*x(i)
      else                 ! no more items
         igeneric=0
      endif
      c(n)=cc
      end

      function iugeneric(n,c,ccop,n1,n2)
      parameter(nx=100000)
      dimension c(*), ccop(*), x(nx)
      common nmax,cost,temp,cmin,rate,x

c if for any value of n an efficient method is available to compute
c the change of statistic under exchange of n1 and n2, the call
c to igeneric should be replaced for that value of n.
c Generically, igeneric is at least O(nmax) while updates are often O(1).
c
c simple but SLOW:
c      call exch(n1,n2)
c      iugeneric=igeneric(n,ccop,nmax,x)
c      call exch(n1,n2)         ! must leave function unchanged
c      end
c
c better:

      cc=c(n)
      iugeneric=1
      if(n.eq.1.or.n.eq.2) then      ! lag 1 and 2 autocorrelations
         if(n1-n.ge.1) cc=cc-x(n1-n)*x(n1)
         if(n2-n.ge.1) cc=cc-x(n2-n)*x(n2)
         if(n1+n.le.nmax) cc=cc-x(n1)*x(n1+n)
         if(n2+n.le.nmax) cc=cc-x(n2)*x(n2+n)
         call exch(n1,n2)
         if(n1-n.ge.1) cc=cc+x(n1-n)*x(n1)
         if(n2-n.ge.1) cc=cc+x(n2-n)*x(n2)
         if(n1+n.le.nmax) cc=cc+x(n1)*x(n1+n)
         if(n2+n.le.nmax) cc=cc+x(n2)*x(n2+n)
      else if(n.eq.3) then ! lag 1,1  3-point-autocorrelation
         if(n1-1.ge.1) cc=cc-x(n1-1)**2*x(n1)
         if(n2-1.ge.1) cc=cc-x(n2-1)**2*x(n2)
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc-x(n1)**2*x(n1+1)
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc-x(n2)**2*x(n2+1)
         call exch(n1,n2)
         if(n1-1.ge.1) cc=cc+x(n1-1)**2*x(n1)
         if(n2-1.ge.1) cc=cc+x(n2-1)**2*x(n2)
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc+x(n1)**2*x(n1+1)
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc+x(n2)**2*x(n2+1)
      else if(n.eq.4) then ! lag 1,0  3-point-autocorrelation
         if(n1-1.ge.1) cc=cc-x(n1-1)*x(n1)**2
         if(n2-1.ge.1) cc=cc-x(n2-1)*x(n2)**2
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc-x(n1)*x(n1+1)**2
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc-x(n2)*x(n2+1)**2
         call exch(n1,n2)
         if(n1-1.ge.1) cc=cc+x(n1-1)*x(n1)**2
         if(n2-1.ge.1) cc=cc+x(n2-1)*x(n2)**2
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc+x(n1)*x(n1+1)**2
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc+x(n2)*x(n2+1)**2
      else if(n.eq.5) then ! lag 2,0  3-point-autocorrelation
         if(n1-2.ge.1) cc=cc-x(n1-2)**2*x(n1)
         if(n2-2.ge.1) cc=cc-x(n2-2)**2*x(n2)
         if(n1+2.le.nmax.and.n1+2.ne.n2) cc=cc-x(n1)**2*x(n1+2)
         if(n2+2.le.nmax.and.n2+2.ne.n1) cc=cc-x(n2)**2*x(n2+2)
         call exch(n1,n2)
         if(n1-2.ge.1) cc=cc+x(n1-2)**2*x(n1)
         if(n2-2.ge.1) cc=cc+x(n2-2)**2*x(n2)
         if(n1+2.le.nmax.and.n1+2.ne.n2) cc=cc+x(n1)**2*x(n1+2)
         if(n2+2.le.nmax.and.n2+2.ne.n1) cc=cc+x(n2)**2*x(n2+2)
      else if(n.eq.6) then ! lag 2,1  3-point-autocorrelation
         if(n1-2.ge.1) cc=cc-x(n1-2)*x(n1-1)*x(n1)
         if(n2-2.ge.1) cc=cc-x(n2-2)*x(n2-1)*x(n2)
         if(n1-1.ge.1.and.n1+1.le.nmax.and.n1+1.ne.n2) 
     .      cc=cc-x(n1-1)*x(n1)*x(n1+1)
         if(n2-1.ge.1.and.n2+1.le.nmax.and.n2+1.ne.n1) 
     .      cc=cc-x(n2-1)*x(n2)*x(n2+1)
         if(n1+2.le.nmax.and.n1+2.ne.n2.and.n1+1.ne.n2) 
     .      cc=cc-x(n1)*x(n1+1)*x(n1+2)
         if(n2+2.le.nmax.and.n2+2.ne.n1.and.n2+1.ne.n1) 
     .      cc=cc-x(n2)*x(n2+1)*x(n2+2)
         call exch(n1,n2)
         if(n1-2.ge.1) cc=cc+x(n1-2)*x(n1-1)*x(n1)
         if(n2-2.ge.1) cc=cc+x(n2-2)*x(n2-1)*x(n2)
         if(n1-1.ge.1.and.n1+1.le.nmax.and.n1+1.ne.n2) 
     .      cc=cc+x(n1-1)*x(n1)*x(n1+1)
         if(n2-1.ge.1.and.n2+1.le.nmax.and.n2+1.ne.n1) 
     .      cc=cc+x(n2-1)*x(n2)*x(n2+1)
         if(n1+2.le.nmax.and.n1+2.ne.n2.and.n1+1.ne.n2) 
     .      cc=cc+x(n1)*x(n1+1)*x(n1+2)
         if(n2+2.le.nmax.and.n2+2.ne.n1.and.n2+1.ne.n1) 
     .      cc=cc+x(n2)*x(n2+1)*x(n2+2)
      else if(n.eq.7) then ! lag 1,1,0  4-point-autocorrelation
         if(n1-1.ge.1) cc=cc-x(n1-1)**2*x(n1)**2
         if(n2-1.ge.1) cc=cc-x(n2-1)**2*x(n2)**2
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc-x(n1)**2*x(n1+1)**2
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc-x(n2)**2*x(n2+1)**2
         call exch(n1,n2)
         if(n1-1.ge.1) cc=cc+x(n1-1)**2*x(n1)**2
         if(n2-1.ge.1) cc=cc+x(n2-1)**2*x(n2)**2
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc+x(n1)**2*x(n1+1)**2
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc+x(n2)**2*x(n2+1)**2
      else if(n.eq.8) then ! lag 1,0,0  4-point-autocorrelation
         if(n1-1.ge.1) cc=cc-x(n1-1)*x(n1)**3
         if(n2-1.ge.1) cc=cc-x(n2-1)*x(n2)**3
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc-x(n1)*x(n1+1)**3
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc-x(n2)*x(n2+1)**3
         call exch(n1,n2)
         if(n1-1.ge.1) cc=cc+x(n1-1)*x(n1)**3
         if(n2-1.ge.1) cc=cc+x(n2-1)*x(n2)**3
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc+x(n1)*x(n1+1)**3
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc+x(n2)*x(n2+1)**3
      else if(n.eq.9) then ! lag 1,1,1  4-point-autocorrelation
         if(n1-1.ge.1) cc=cc-x(n1-1)**3*x(n1)
         if(n2-1.ge.1) cc=cc-x(n2-1)**3*x(n2)
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc-x(n1)**3*x(n1+1)
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc-x(n2)**3*x(n2+1)
         call exch(n1,n2)
         if(n1-1.ge.1) cc=cc+x(n1-1)**3*x(n1)
         if(n2-1.ge.1) cc=cc+x(n2-1)**3*x(n2)
         if(n1+1.le.nmax.and.n1+1.ne.n2) cc=cc+x(n1)**3*x(n1+1)
         if(n2+1.le.nmax.and.n2+1.ne.n1) cc=cc+x(n2)**3*x(n2+1)
      else                 ! no more items
         iugeneric=0
         return
      endif
      ccop(n)=cc
      call exch(n1,n2)
      end

      subroutine aver(cav,dc)
c      cav=max(cav,abs(dc))      ! max (L-infinity) norm
      cav=cav+abs(dc)           ! L-1 norm (sum of moduli)
c      cav=cav+dc**2             ! L-2 norm (sum of squares)
      end

