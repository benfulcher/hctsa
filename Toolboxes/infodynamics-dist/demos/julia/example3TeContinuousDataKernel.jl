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

# = Example 3 - Transfer entropy on continuous data using kernel estimators =

# Simple transfer entropy (TE) calculation on continuous-valued data using the (box) kernel-estimator TE calculator.

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
teClass = @jimport infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel;
teCalc=teClass(());
jcall(teCalc, "setProperty", Void, (JString, JString), "NORMALISE", "true"); # Normalise the individual variables
jcall(teCalc, "initialise", Void, (jint, jdouble), 1, 0.5); # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
jcall(teCalc, "setObservations", Void, (Array{jdouble,1}, Array{jdouble,1}),
				     sourceArray, destArray);
# For copied source, should give something close to expected value for correlated Gaussians:
result = jcall(teCalc, "computeAverageLocalOfObservations", jdouble, ());
@printf("TE result %.4f bits; expected to be close to %.4f bits for these correlated Gaussians but biased upwards\n", result, log(1/(1-covariance^2))/log(2));

jcall(teCalc, "initialise", Void, ()); # Initialise leaving the parameters the same
jcall(teCalc, "setObservations", Void, (Array{jdouble,1}, Array{jdouble,1}),
				     sourceArray2, destArray);
# For random source, it should give something close to 0 bits
result2 = jcall(teCalc, "computeAverageLocalOfObservations", jdouble, ());
@printf("TE result %.4f bits; expected to be close to 0 bits for uncorrelated Gaussians but will be biased upwards\n", result2);

# We can get insight into the bias by examining the null distribution:
empiricalDistClass = @jimport infodynamics.utils.EmpiricalMeasurementDistribution;
nullDist = jcall(teCalc, "computeSignificance", empiricalDistClass, (jint,), 100);
@printf("Null distribution for unrelated source and destination (i.e. the bias) has mean %.4f and standard deviation %.4f\n",
	jcall(nullDist, "getMeanOfDistribution", jdouble, ()),
	jcall(nullDist, "getStdOfDistribution", jdouble, ()));

