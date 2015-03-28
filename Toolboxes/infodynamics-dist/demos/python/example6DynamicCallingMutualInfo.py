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

# = Example 6 - Mutual information calculation with dynamic specification of calculator =

# This example shows how to write Python code to take advantage of the
# common interfaces defined for various information-theoretic calculators.
# Here, we use the common form of the infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  interface (which is never named here) to write common code into which we can plug
#  one of three concrete implementations (kernel estimator, Kraskov estimator or
#  linear-Gaussian estimator) by dynamically supplying the class name of
#  the concrete implementation.

from jpype import *
import random
import string
import numpy
import readFloatsFile

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

#---------------------
# 1. Properties for the calculation (these are dynamically changeable):
# The name of the data file (relative to this directory)
datafile = '../data/4ColsPairedNoisyDependence-1.txt'
# List of column numbers for univariate time seres 1 and 2:
#  (you can select any columns you wish to be contained in each variable)
univariateSeries1Column = 0; # array indices start from 0 in python
univariateSeries2Column = 2;
# List of column numbers for joint variables 1 and 2:
#  (you can select any columns you wish to be contained in each variable)
jointVariable1Columns = [0,1] # array indices start from 0 in python
jointVariable2Columns = [2,3]
# The name of the concrete implementation of the interface 
#  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  which we wish to use for the calculation.
# Note that one could use any of the following calculators (try them all!):
#  implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1" # MI(1;3) as 0.10044, MI([1,2], [3,4]) = 0.36353 (NATS not bits)
#  implementingClass = "infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel"
#  implementingClass = "infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian"
implementingClass = "infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1"

#---------------------
# 2. Load in the data
data = readFloatsFile.readFloatsFile(datafile)
# As numpy array:
A = numpy.array(data)
# Pull out the columns from the data set for a univariate MI calculation:
univariateSeries1 = A[:,univariateSeries1Column];
univariateSeries2 = A[:,univariateSeries2Column];
# Pull out the columns from the data set for a multivariate MI calculation:
jointVariable1 = A[:,jointVariable1Columns]
jointVariable2 = A[:,jointVariable2Columns]

#--------------------
# 3. Dynamically instantiate an object of the given class:
# (in fact, all java object creation in python is dynamic - it has to be,
#  since the languages are interpreted. This makes our life slightly easier at this
#  point than it is in demos/java/example6 where we have to handle this manually)
indexOfLastDot = string.rfind(implementingClass, ".")
implementingPackage = implementingClass[:indexOfLastDot]
implementingBaseName = implementingClass[indexOfLastDot+1:]
miCalcClass = eval('JPackage(\'%s\').%s' % (implementingPackage, implementingBaseName))
miCalc = miCalcClass()

#--------------------
# 4. Start using the MI calculator, paying attention to only
#  call common methods defined in the interface type
#  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
#  not methods only defined in a given implementation class.
# a. Initialise the calculator for a univariate calculation:
miCalc.initialise(1, 1)
# b. Supply the observations to compute the PDFs from:
miCalc.setObservations(univariateSeries1, univariateSeries2)
# c. Make the MI calculation:
miUnivariateValue = miCalc.computeAverageLocalOfObservations()

#---------------------
# 5. Continue onto a multivariate calculation, still
#    only calling common methods defined in the interface type.
# a. Initialise the calculator for a multivariate calculation
#   to use the required number of dimensions for each variable:
miCalc.initialise(len(jointVariable1Columns), len(jointVariable2Columns))
# b. Supply the observations to compute the PDFs from:
miCalc.setObservations(jointVariable1, jointVariable2)
# c. Make the MI calculation:
miJointValue = miCalc.computeAverageLocalOfObservations()

print("MI calculator %s computed the univariate MI(%d;%d) as %.5f and joint MI([%s];[%s]) as %.5f\n" %
	(implementingClass, univariateSeries1Column, univariateSeries2Column, miUnivariateValue,
	str(jointVariable1Columns).strip('[]'), str(jointVariable2Columns).strip('[]'), miJointValue))

