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
c   help.f
c   Utilities for usage message
c   author T. Schreiber (1998)
c===========================================================================
      subroutine whatido(text,iverb)
      character*72 progname
      character*(*) text

      call getarg(0,progname)
      call argdel(0)
      iverb=igetv(iverb)
      if(iv_io(iverb).eq.1) then
         write(istderr(),'()')  
         write(istderr(),'(a)') 
     .              "TISEAN 3.0.1 (C) R. Hegger, H. Kantz, T. Schreiber
     .(1998-2007)"
         write(istderr(),'()')  
         write(istderr(),'(a,a,a)') 
     .      progname(1:index(progname," ")-1), ": ", text
      endif
      if(lopt("h",1).eq.1) call usage()
      end

      subroutine whatineed(text)
      character*72 progname
      character*(*) text

      call getarg(0,progname)
      write(istderr(),'()') 
      write(istderr(),'(a,a,x,a)') 
     .   "Usage: ", progname(1:index(progname," ")-1),  text
      end

      subroutine popt(c,text)
      character*(*) c,text

      write(istderr(),'(5h    -,a,x,1h<,a,1h>)') c, text
      end

      subroutine ptext(text)
      character*(*) text

      write(istderr(),'(3x,a)') text
      end

      subroutine pout(text)
      character*(*) text

      write(istderr(),'(8h    -o <,a,a,1h>)') 
     .   "output file name, just -o means ", text
      end

      subroutine pall()

      call popt("V","verbosity level (0 = only fatal errors)")
      call popt("h","show this message")
      write(istderr(),'()')
      end

