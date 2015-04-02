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

# = Example 1 - Transfer entropy on binary data =

# Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator:

# Import the JavaCall package:
using JavaCall;

# Change location of jar to match yours:
jarLocation = "../../infodynamics.jar";
# Start the JVM supplying classpath and heap size
#  (increase memory here if you get crashes due to not enough space)
JavaCall.init(["-Djava.class.path=$(jarLocation)", "-Xmx128M"]);

# Generate some random binary data.
sourceArray=rand(0:1, 100);
destArray = [0; sourceArray[1:99]];
sourceArray2=rand(0:1, 100);

# Create a TE calculator and run it:
teClass = @jimport infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;
teCalc=teClass((jint, jint), 2, 1);
jcall(teCalc, "initialise", Void, ()); # This is how to indicate a void return type and arguments

# We can pass simple arrays of ints directly in:
jcall(teCalc, "addObservations", Void, (Array{jint,1}, Array{jint,1}),
				     sourceArray, destArray);
result = jcall(teCalc, "computeAverageLocalOfObservations", jdouble, ());
@printf("For copied source, result should be close to 1 bit : %.4f\n", result);

jcall(teCalc, "initialise", Void, ());
jcall(teCalc, "addObservations", Void, (Array{jint,1}, Array{jint,1}),
				     sourceArray2, destArray);
result2 = jcall(teCalc, "computeAverageLocalOfObservations", jdouble, ());
@printf("For random source, result should be close to 0 bits: %.4f\n", result2);

