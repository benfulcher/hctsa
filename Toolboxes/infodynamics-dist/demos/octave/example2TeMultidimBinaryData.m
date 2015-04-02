%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2012, Joseph T. Lizier
%%  
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%  
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%  
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

% = Example 2 - Transfer entropy on multidimensional binary data =

% Simple transfer entropy (TE) calculation on multidimensional binary data using the discrete TE calculator.

% This example is important for Octave users, because it shows how to handle multidimensional arrays from Octave to Java (this is not as simple as single dimensional arrays in example 1 - it requires using supplied scripts to convert the array).

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Create many columns in a multidimensional array (2 rows by 100 columns),
%  where the next time step (row 2) copies the value of the column on the left
%  from the previous time step (row 1):
twoDTimeSeriesOctave       = (rand(1, 100)>0.5)*1;
twoDTimeSeriesOctave(2, :) = [twoDTimeSeriesOctave(1,100), twoDTimeSeriesOctave(1, 1:99)];

% Things get a little tricky if we want to pass 2D arrays into Java.
% Unlike native Octave 1D arrays in Example 1, 
%  native Octave 2D+ arrays do not seem to get directly converted to java arrays,
%  so we use the supplied scripts to make the conversion (via org.octave.Matrix class in octave)
% Matlab handles the conversion automatically, so in Matlab this script just returns
%  the array that was passed in.
twoDTimeSeriesJavaInt = octaveToJavaIntMatrix(twoDTimeSeriesOctave);

% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 1);
teCalc.initialise();
% Add observations of transfer across one cell to the right per time step:
teCalc.addObservations(twoDTimeSeriesJavaInt, 1);
fprintf('The result should be close to 1 bit here, since we are executing copy operations of what is effectively a random bit to each cell here: ');
result2D = teCalc.computeAverageLocalOfObservations()

