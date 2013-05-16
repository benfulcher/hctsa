/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/* Author: Rainer Hegger Last modified: Jul 15, 1999 */
#include <stdlib.h>
#include <stdio.h>
#include "tisean_cec.h"

void check_alloc(void *pnt)
{
  if (pnt == NULL) {
    fprintf(stderr,"check_alloc: Couldn't allocate enough memory. Exiting\n");
    exit(CHECK_ALLOC_NOT_ENOUGH_MEMORY);
  }
}
