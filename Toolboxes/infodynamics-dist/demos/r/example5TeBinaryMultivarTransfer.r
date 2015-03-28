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

# = Example 5 - Multivariate transfer entropy on binary data =

# Multivariate transfer entropy (TE) calculation on binary data using the discrete TE calculator:

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# Change location of jar to match yours:
#  IMPORTANT -- If using the default below, make sure you have set the working directory
#   in R (e.g. with setwd()) to the location of this file (i.e. demos/r) !!
.jaddClassPath("../../infodynamics.jar")

# Generate some random binary data.
numObservations <- 100
sourceArray<-matrix(sample(0:1,numObservations*2, replace=TRUE),numObservations,2)
sourceArray2<-matrix(sample(0:1,numObservations*2, replace=TRUE),numObservations,2)
# Destination variable takes a copy of the first bit of the source in bit 1,
#  and an XOR of the two bits of the source in bit 2:
destArray <- cbind( c(0L, sourceArray[1:numObservations-1,1]), # column 1
		    c(0L, 1L*xor(sourceArray[1:numObservations-1,1],
		                 sourceArray[1:numObservations-1,2]))) # column 2

# Convert the 2D arrays to Java format:
sourceArrayJava <- .jarray(sourceArray, "[I", dispatch=TRUE)
sourceArray2Java <- .jarray(sourceArray2, "[I", dispatch=TRUE)
destArrayJava <- .jarray(destArray, "[I", dispatch=TRUE)

# Create a TE calculator and run it:
teCalc<-.jnew("infodynamics/measures/discrete/TransferEntropyCalculatorDiscrete", 4L, 1L)
.jcall(teCalc,"V","initialise") # V for void return value
# We need to construct the joint values for the dest and source before we pass them in,
#  and need to use the matrix conversion routine when calling from Matlab/Octave:
mUtils<-.jnew("infodynamics/utils/MatrixUtils")
.jcall(teCalc,"V","addObservations",
		.jcall(mUtils,"[I","computeCombinedValues", sourceArrayJava, 2L),
		.jcall(mUtils,"[I","computeCombinedValues", destArrayJava, 2L))
result<-.jcall(teCalc,"D","computeAverageLocalOfObservations")
cat("For source which the 2 bits are determined from, result should be close to 2 bits : ", result, "\n")

.jcall(teCalc,"V","initialise")
.jcall(teCalc,"V","addObservations",
		.jcall(mUtils,"[I","computeCombinedValues", sourceArray2Java, 2L),
		.jcall(mUtils,"[I","computeCombinedValues", destArrayJava, 2L))
result2<-.jcall(teCalc,"D","computeAverageLocalOfObservations")
cat("For random source, result should be close to 0 bits in theory: ", result2, "\n");
cat("Result for random source is inflated towards 0.3 due to finite observation length ",
    .jcall(teCalc,"I","getNumObservations"), "\n",
    "One can verify that the answer is consistent with that from a\n",
    "random source by checking: teCalc.computeSignificance(1000); ans.pValue\n");

