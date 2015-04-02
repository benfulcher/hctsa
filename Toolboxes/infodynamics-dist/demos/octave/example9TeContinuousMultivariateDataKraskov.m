%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2014, Viola Priesemann, Joseph T. Lizier
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

% = Example 9 - Transfer entropy on continuous multivariate data using Kraskov estimators =

% Transfer entropy (TE) calculation on multivariate continuous-valued data using the Kraskov-estimator TE calculator.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar')

% Generate some random normalised data.
numObservations = 10000;
covariance=0.4; 

% Define the dimension of the states of the RVs
sourceDim = 2;  
destDim = 3;

sourceMVArray = randn(numObservations, sourceDim);
% Set first two columns of dest to copy source values
destMVArray  = [zeros(1,sourceDim); covariance*(sourceMVArray(1:numObservations-1,:)) + (1-covariance)*randn(numObservations-1, sourceDim)];
% Set a third colum to be randomised
destMVArray(:,3) = randn(numObservations, 1);
sourceMVArray2= randn(numObservations, sourceDim); % Uncorrelated source

% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov');
teCalc.initialise(1,sourceDim,destDim); % Use history length 1 (Schreiber k=1)
teCalc.setProperty('k', '4'); % Use Kraskov parameter K=4 for 4 nearest points
teCalc.setObservations(octaveToJavaDoubleMatrix(sourceMVArray), octaveToJavaDoubleMatrix(destMVArray));
% Perform calculation with correlated source:
result = teCalc.computeAverageLocalOfObservations();
% Note that the calculation is a random variable (because the generated
%  data is a set of random variables) - the result will be of the order
%  of what we expect, but not exactly equal to it; in fact, there will
%  be some variance around it. It will probably be biased down here
%  due to small correlations between the supposedly uncorrelated variables.
fprintf('TE result %.4f nats; expected to be close to %.4f nats for the two correlated Gaussians\n', ...
result, 2*log(1/(1-covariance^2)));

% Perform calculation with uncorrelated source:
teCalc.initialise(1,sourceDim,destDim); % Initialise leaving the parameters the same
teCalc.setObservations(octaveToJavaDoubleMatrix(sourceMVArray2), octaveToJavaDoubleMatrix(destMVArray));
result2 = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f nats; expected to be close to 0 nats for these uncorrelated Gaussians\n', result2);
clear teCalc

