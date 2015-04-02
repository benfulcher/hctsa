##
##  Java Information Dynamics Toolkit (JIDT)
##  Copyright (C) 2012, Joseph T. Lizier
##  
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

def readIntsFile(filename):
	"Read a 2D array of int from a given file"
	with open(filename) as f:
 	    # Space separate numbers, one time step per line, each column is a variable
 	    array = []
	    for line in f: # read all lines
	    	if (line.startswith("%") or line.startswith("#")):
	    		# Assume this is a comment line
	    		continue;
	    	if (len(line.split()) == 0):
	    		# Line is empty
	    		continue;
		array.append([int(x) for x in line.split()])
	    return array
    
