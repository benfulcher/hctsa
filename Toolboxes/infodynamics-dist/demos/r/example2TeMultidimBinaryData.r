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

# = Example 2 - Transfer entropy on multidimensional binary data =

# Simple transfer entropy (TE) calculation on multidimensional binary data using the discrete TE calculator.

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# Change location of jar to match yours:
#  IMPORTANT -- If using the default below, make sure you have set the working directory
#   in R (e.g. with setwd()) to the location of this file (i.e. demos/r) !!
.jaddClassPath("../../infodynamics.jar")

# Create many columns in a multidimensional array (2 rows by 100 columns),
#  where the next time step (row 2) copies the value of the column on the left
#  from the previous time step (row 1):
twoDTimeSeriesRtime1 <- sample(0:1, 100, replace="TRUE")
twoDTimeSeriesRtime2 <- c(twoDTimeSeriesRtime1[100], twoDTimeSeriesRtime1[1:99])
twoDTimeSeriesR <- rbind(twoDTimeSeriesRtime1, twoDTimeSeriesRtime2)

# Create a TE calculator and run it:
teCalc<-.jnew("infodynamics/measures/discrete/TransferEntropyCalculatorDiscrete", 2L, 1L)
.jcall(teCalc,"V","initialise") # V for void return value
# Add observations of transfer across one cell to the right per time step:
twoDTimeSeriesJava <- .jarray(twoDTimeSeriesR, "[I", dispatch=TRUE)
.jcall(teCalc,"V","addObservations", twoDTimeSeriesJava, 1L)
result2D <- .jcall(teCalc,"D","computeAverageLocalOfObservations")
cat("The result should be close to 1 bit here, since we are executing copy operations of what is effectively a random bit to each cell here: ", result2D, "\n")

