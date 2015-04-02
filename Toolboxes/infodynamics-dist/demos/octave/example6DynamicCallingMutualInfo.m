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

% Example 6 - Mutual information calculation with dynamic specification of calculator

% This example shows how to write Matlab/Octave code to take advantage of the
%  common interfaces defined for various information-theoretic calculators.
% Here, we use the common form of the infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
%  interface (which is never named here) to write common code into which we can plug
%  one of three concrete implementations (kernel estimator, Kraskov estimator or
%  linear-Gaussian estimator) by dynamically supplying the class name of
%  the concrete implementation.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

%---------------------
% 1. Properties for the calculation (these are dynamically changeable, you could
%    load them in from another properties file):
% The name of the data file (relative to this directory)
datafile = '../data/4ColsPairedNoisyDependence-1.txt';
% List of column numbers for univariate time seres 1 and 2:
%  (you can select any columns you wish to be contained in each variable)
univariateSeries1Column = 1; % array indices start from 1 in octave/matlab
univariateSeries2Column = 3;
% List of column numbers for joint variables 1 and 2:
%  (you can select any columns you wish to be contained in each variable)
jointVariable1Columns = [1,2]; % array indices start from 1 in octave/matlab
jointVariable2Columns = [3,4];
% The name of the concrete implementation of the interface 
%  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
%  which we wish to use for the calculation.
% Note that one could use any of the following calculators (try them all!):
%  implementingClass = 'infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1'; % MI(1;3) as 0.10044, MI([1,2], [3,4]) = 0.36353 (NATS not bits)
%  implementingClass = 'infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel';
%  implementingClass = 'infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian';
implementingClass = 'infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1';

%---------------------
% 2. Load in the data
data = load(datafile);
% Pull out the columns from the data set for a univariate MI calculation:
univariateSeries1 = data(:, univariateSeries1Column);
univariateSeries2 = data(:, univariateSeries2Column);
% Pull out the columns from the data set for a multivariate MI calculation:
jointVariable1 = data(:, jointVariable1Columns);
jointVariable2 = data(:, jointVariable2Columns);

%---------------------
% 3. Dynamically instantiate an object of the given class:
% (in fact, all java object creation in octave/matlab is dynamic - it has to be,
%  since the languages are interpreted. This makes our life slightly easier at this
%  point than it is in demos/java/example6 where we have to handle this manually)
miCalc = javaObject(implementingClass);

%---------------------
% 4. Start using the MI calculator, paying attention to only
%  call common methods defined in the interface type
%  infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate
%  not methods only defined in a given implementation class.
% a. Initialise the calculator for a univariate calculation:
miCalc.initialise(1, 1);
% b. Supply the observations to compute the PDFs from:
miCalc.setObservations(octaveToJavaDoubleArray(univariateSeries1), octaveToJavaDoubleArray(univariateSeries2));
% c. Make the MI calculation:
miUnivariateValue = miCalc.computeAverageLocalOfObservations();

%---------------------
% 5. Continue onto a multivariate calculation, still
%    only calling common methods defined in the interface type.
% a. Initialise the calculator for a multivariate calculation
%   to use the required number of dimensions for each variable:
miCalc.initialise(length(jointVariable1Columns), length(jointVariable2Columns));
% b. Supply the observations to compute the PDFs from:
%    (Note the different method call octaveToJavaDoubleMatrix to convert a *multivariate*
%     time-series to Java format here)
miCalc.setObservations(octaveToJavaDoubleMatrix(jointVariable1), octaveToJavaDoubleMatrix(jointVariable2));
% c. Make the MI calculation:
miJointValue = miCalc.computeAverageLocalOfObservations();

fprintf('MI calculator %s\ncomputed the univariate MI(%d;%d) as %.5f and joint MI([%s];[%s]) as %.5f\n', ...
		implementingClass, univariateSeries1Column, univariateSeries2Column, miUnivariateValue, ...
		sprintf('%d,', jointVariable1Columns), sprintf('%d,', jointVariable2Columns), miJointValue);

