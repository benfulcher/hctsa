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

# This example is important for Python users using JPype, because it shows how to handle multidimensional arrays from Python to Java.

from jpype import *
import random

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar"
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)

# Create many columns in a multidimensional array, e.g. for fully random values:
# twoDTimeSeriesOctave = [[random.randint(0,1) for y in xrange(2)] for x in xrange(10)] # for 10 rows (time-steps) for 2 variables

# However here we want 2 rows by 100 columns where the next time step (row 2) is to copy the
# value of the column on the left from the previous time step (row 1):
numObservations = 100
row1 = [random.randint(0,1) for r in xrange(numObservations)]
row2 = [row1[numObservations-1]] + row1[0:numObservations-1] # Copy the previous row, offset one column to the right
twoDTimeSeriesPython = []
twoDTimeSeriesPython.append(row1)
twoDTimeSeriesPython.append(row2)
twoDTimeSeriesJavaInt = JArray(JInt, 2)(twoDTimeSeriesPython); # 2 indicating 2D array

# Create a TE calculator and run it:
teCalcClass = JPackage("infodynamics.measures.discrete").TransferEntropyCalculatorDiscrete
teCalc = teCalcClass(2,1)
teCalc.initialise()
# Add observations of transfer across one cell to the right per time step:
teCalc.addObservations(twoDTimeSeriesJavaInt, 1)
result2D = teCalc.computeAverageLocalOfObservations()
print(('The result should be close to 1 bit here, since we are executing copy ' + \
      'operations of what is effectively a random bit to each cell here: %.3f ' + \
      'bits from %d observations') % (result2D, teCalc.getNumObservations()))

