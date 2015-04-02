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

# Import the JavaCall package:
using JavaCall;

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar";
# Start the JVM supplying classpath and heap size
#  (increase memory here if you get crashes due to not enough space)
JavaCall.init(["-Djava.class.path=$(jarLocation)", "-Xmx128M"]);

# Generate some random normalised data.
numObservations = 1000;
covariance=0.4;
sourceArray=randn(numObservations);
destArray = [0, covariance*sourceArray[1:numObservations-1] + (1-covariance)*randn(numObservations - 1)];
sourceArray2=randn(numObservations); # Uncorrelated source

# Create a TE calculator and run it:
teClass = @jimport infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;
teCalc=teClass(());
jcall(teCalc, "setProperty", Void, (JString, JString), "k", "4"); # Use Kraskov parameter K=4 for 4 nearest points
jcall(teCalc, "initialise", Void, (jint, ), 1); # Use history length 1 (Schreiber k=1)
# Perform calculation with correlated source:
jcall(teCalc, "setObservations", Void, (Array{jdouble,1}, Array{jdouble,1}),
				     sourceArray, destArray);
result = jcall(teCalc, "computeAverageLocalOfObservations", jdouble, ());
# Note that the calculation is a random variable (because the generated
#  data is a set of random variables) - the result will be of the order
#  of what we expect, but not exactly equal to it; in fact, there will
#  be a large variance around it.
@printf("TE result %.4f nats; expected to be close to %.4f nats for these correlated Gaussians\n",
    result, log(1/(1-covariance^2)));

# Perform calculation with uncorrelated source:
jcall(teCalc, "initialise", Void, ()); # Initialise leaving the parameters the same
jcall(teCalc, "setObservations", Void, (Array{jdouble,1}, Array{jdouble,1}),
				     sourceArray2, destArray);
result2 = jcall(teCalc, "computeAverageLocalOfObservations", jdouble, ());
@printf("TE result %.4f nats; expected to be close to 0 nats for these uncorrelated Gaussians\n", result2);

# We can also compute the local TE values for the time-series samples here:
#  (See more about utility of local TE in the CA demos)
localTE = jcall(teCalc, "computeLocalOfPreviousObservations", Array{jdouble,1}, ());
@printf("Notice that the mean of locals, %.4f nats, equals the previous result\n",
	sum(localTE)/(numObservations-1));

