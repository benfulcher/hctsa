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

from jpype import *
import random
import math

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Generate some random normalised data.
numObservations = 1000;
covariance=0.4;
# Source array of random normals:
sourceArray = [random.normalvariate(0,1) for r in xrange(numObservations)];
# Destination array of random normals with partial correlation to previous value of sourceArray
destArray = [0] + [sum(pair) for pair in zip([covariance*y for y in sourceArray[0:numObservations-1]], \
                                             [(1-covariance)*y for y in [random.normalvariate(0,1) for r in xrange(numObservations-1)]] ) ];
# Uncorrelated source array:
sourceArray2 = [random.normalvariate(0,1) for r in xrange(numObservations)];
# Create a TE calculator and run it:
teCalcClass = JPackage("infodynamics.measures.continuous.kraskov").TransferEntropyCalculatorKraskov
teCalc = teCalcClass();
teCalc.setProperty("NORMALISE", "true"); # Normalise the individual variables
teCalc.initialise(1); # Use history length 1 (Schreiber k=1)
teCalc.setProperty("k", "4"); # Use Kraskov parameter K=4 for 4 nearest points
# Perform calculation with correlated source:
teCalc.setObservations(JArray(JDouble, 1)(sourceArray), JArray(JDouble, 1)(destArray));
result = teCalc.computeAverageLocalOfObservations();
# Note that the calculation is a random variable (because the generated
#  data is a set of random variables) - the result will be of the order
#  of what we expect, but not exactly equal to it; in fact, there will
#  be a large variance around it.
print("TE result %.4f nats; expected to be close to %.4f nats for these correlated Gaussians" % \
    (result, math.log(1/(1-math.pow(covariance,2)))));
# Perform calculation with uncorrelated source:
teCalc.initialise(); # Initialise leaving the parameters the same
teCalc.setObservations(JArray(JDouble, 1)(sourceArray2), JArray(JDouble, 1)(destArray));
result2 = teCalc.computeAverageLocalOfObservations();
print("TE result %.4f nats; expected to be close to 0 nats for these uncorrelated Gaussians" % result2);


