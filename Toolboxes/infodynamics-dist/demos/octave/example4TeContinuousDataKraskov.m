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

% = Example 4 - Transfer entropy on continuous data using Kraskov estimators =

% Simple transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator.

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Generate some random normalised data.
numObservations = 1000;
covariance=0.4;
sourceArray=randn(numObservations, 1);
destArray = [0; covariance*sourceArray(1:numObservations-1) + (1-covariance)*randn(numObservations - 1, 1)];
sourceArray2=randn(numObservations, 1); % Uncorrelated source
% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
teCalc.setProperty('k', '4'); % Use Kraskov parameter K=4 for 4 nearest points
teCalc.initialise(1); % Use history length 1 (Schreiber k=1)
% Perform calculation with correlated source:
teCalc.setObservations(sourceArray, destArray);
result = teCalc.computeAverageLocalOfObservations();
% Note that the calculation is a random variable (because the generated
%  data is a set of random variables) - the result will be of the order
%  of what we expect, but not exactly equal to it; in fact, there will
%  be a large variance around it.
fprintf('TE result %.4f nats; expected to be close to %.4f nats for these correlated Gaussians\n', ...
    result, log(1/(1-covariance^2)));
% Perform calculation with uncorrelated source:
teCalc.initialise(); % Initialise leaving the parameters the same
teCalc.setObservations(sourceArray2, destArray);
result2 = teCalc.computeAverageLocalOfObservations();
fprintf('TE result %.4f nats; expected to be close to 0 nats for these uncorrelated Gaussians\n', result2);

% We can also compute the local TE values for the time-series samples here:
%  (See more about utility of local TE in the CA demos)
localTE = teCalc.computeLocalOfPreviousObservations();
fprintf('Notice that the mean of locals, %.4f nats, equals the previous result\n', ...
	sum(javaMatrixToOctave(localTE))/(numObservations-1));

