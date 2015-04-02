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

# = Example 4 - Transfer entropy on continuous data using Kraskov estimators =

# Simple transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator.

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# Change location of jar to match yours:
#  IMPORTANT -- If using the default below, make sure you have set the working directory
#   in R (e.g. with setwd()) to the location of this file (i.e. demos/r) !!
.jaddClassPath("../../infodynamics.jar")

# Generate some random normalised data.
numObservations<-1000
covariance<-0.4
sourceArray<-rnorm(numObservations)
destArray = c(0, covariance*sourceArray[1:numObservations-1] + (1-covariance)*rnorm(numObservations-1, 0, 1))
sourceArray2<-rnorm(numObservations) # Uncorrelated source

# Create a TE calculator:
teCalc<-.jnew("infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov")
.jcall(teCalc,"V","setProperty", "k", "4") # Use Kraskov parameter K=4 for 4 nearest points

# Perform calculation with correlated source:
.jcall(teCalc,"V","initialise", 1L) # Use history length 1 (Schreiber k=1)
.jcall(teCalc,"V","setObservations", sourceArray, destArray)
result <- .jcall(teCalc,"D","computeAverageLocalOfObservations")
# Note that the calculation is a random variable (because the generated
#  data is a set of random variables) - the result will be of the order
#  of what we expect, but not exactly equal to it; in fact, there will
#  be a large variance around it.
cat("TE result ",  result, "nats; expected to be close to ", log(1/(1-covariance^2)), " nats for these correlated Gaussians\n")

# Perform calculation with uncorrelated source:
.jcall(teCalc,"V","initialise") # Initialise leaving the parameters the same
.jcall(teCalc,"V","setObservations", sourceArray2, destArray)
result2 <- .jcall(teCalc,"D","computeAverageLocalOfObservations")
cat("TE result ",  result2, "nats; expected to be close to 0 nats for uncorrelated Gaussians\n")

# We can also compute the local TE values for the time-series samples here:
#  (See more about utility of local TE in the CA demos)
localTE <- .jcall(teCalc,"[D","computeLocalOfPreviousObservations")
cat("Notice that the mean of locals", sum(localTE)/(numObservations-1),
	"nats equals the above result\n")

