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

# Example 6 - Mutual information calculation with dynamic specification of calculator

# This example shows how to write Julia code to take advantage of the
#  common interfaces defined for various information-theoretic calculators.
# Here, we use the common form of the infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  interface (which is never named here) to write common code into which we can plug
#  one of three concrete implementations (kernel estimator, Kraskov estimator or
#  linear-Gaussian estimator) by dynamically supplying the class name of
#  the concrete implementation.

# Import the JavaCall package:
using JavaCall;

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar";
# Start the JVM supplying classpath and heap size
#  (increase memory here if you get crashes due to not enough space)
JavaCall.init(["-Djava.class.path=$(jarLocation)", "-Xmx128M"]);

#---------------------
# 1. Properties for the calculation (these are dynamically changeable, you could
#    load them in from another properties file):
# The name of the data file (relative to this directory)
datafile = "../data/4ColsPairedNoisyDependence-1.txt";
# List of column numbers for variables 1 and 2:
#  (you can select any columns you wish to be contained in each variable)
variable1Columns = [1,2]; # array indices start from 1 in Julia
variable2Columns = [3,4];
# The name of the concrete implementation of the interface 
#  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  which we wish to use for the calculation.
# Note that one could use any of the following calculators (try them all!):
#  implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1"; % MI([1,2], [3,4]) = 0.36353
#  implementingClass = "infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel";
#  implementingClass = "infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian";
implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1";

#---------------------
# 2. Load in the data
data = readdlm(datafile, ' ', '\n')
# Pull out the columns from the data set which correspond to each of variable 1 and 2:
variable1 = data[:, variable1Columns];
variable2 = data[:, variable2Columns];

#---------------------
# 3. Dynamically instantiate an object of the given class:
# Since @jimport uses the provided text directly (rather than evaluating a variable)
#  we need to do some Julia interpolation of $implementingClass to get our
#  dynamic reference to the class.
miCalcClass = eval(:(@jimport $implementingClass));
miCalc = miCalcClass(());

#---------------------
# 4. Start using the MI calculator, paying attention to only
#  call common methods defined in the interface type
#  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  not methods only defined in a given implementation class.
# a. Initialise the calculator to use the required number of
#   dimensions for each variable:
jcall(miCalc, "initialise", Void, (jint,jint), length(variable1Columns), length(variable2Columns));
# b. Supply the observations to compute the PDFs from:

# For the moment we have to exit this example here::
@printf("We're stuck at this point until support for multidimensional arrays is included in JavaCall");
exit();

jcall(miCalc, "setObservations", Void, (Array{jdouble,2}, Array{jdouble,2}), variable1, variable2);
# c. Make the MI calculation:
miValue = jcall(miCalc, "computeAverageLocalOfObservations", jdouble, ());
@printf("MI calculator %s\ncomputed the joint MI as %.5f\n",
		implementingClass, miValue);

