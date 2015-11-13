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

% = Example 3 - Transfer entropy on continuous data using kernel estimators =

% Simple transfer entropy (TE) calculation on continuous-valued data using the (box) kernel-estimator TE calculator.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Generate some random normalised data.
numObservations = 1000;
covariance=0.4;
sourceArray=randn(numObservations, 1);
destArray = [0; covariance*sourceArray(1:numObservations-1) + (1-covariance)*randn(numObservations - 1, 1)];
sourceArray2=randn(numObservations, 1); % Uncorrelated source
% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
teCalc.setProperty('NORMALISE', 'true'); % Normalise the individual variables
teCalc.initialise(1, 0.5); % Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
teCalc.setObservations(sourceArray, destArray);
% For copied source, should give something close to expected value for correlated Gaussians:
result = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f bits; expected to be close to %.4f bits for these correlated Gaussians but biased upwards\n', ...
    result, log(1/(1-covariance^2))/log(2));
teCalc.initialise(); % Initialise leaving the parameters the same
teCalc.setObservations(sourceArray2, destArray);
% For random source, it should give something close to 0 bits
result2 = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f bits; expected to be close to 0 bits for uncorrelated Gaussians but will be biased upwards\n', ...
    result2);

% We can get insight into the bias by examining the null distribution:
nullDist = teCalc.computeSignificance(100);
fprintf(['Null distribution for unrelated source and destination ', ...
	'(i.e. the bias) has mean %.4f and standard deviation %.4f\n'], ...
	nullDist.getMeanOfDistribution(), nullDist.getStdOfDistribution());

