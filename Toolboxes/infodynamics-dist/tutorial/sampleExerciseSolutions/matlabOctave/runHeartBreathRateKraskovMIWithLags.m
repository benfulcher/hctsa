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

% function [miHeartToBreath] = runHeartBreathRateKraskovMIWithLags(lags, numSurrogates)
%
% runHeartBreathRateKraskovMI
% Version 1.0
% Joseph Lizier
% 3/2/2015
%
% Used to explore mutual information in the heart rate / breath rate example of Schreiber --
%  but estimates MI using Kraskov-Stoegbauer-Grassberger estimation.
%
% This script is used to address challenge 1, incorporating lags between heart and breath data.
%
% Note that the paths (to libraries, data, etc) are set assuming this code is run in
% the folder tutorial/sampleExercises/matlabOctave, relative to the main folder of the JIDT distribution.
%
%
% Inputs
% - lags - a scalar specifying a single, or vector specifying multiple, lag from heart to breath to evaluate MI with.
%    These lags may be negative, implying a (positive) lag from breath to heart.
% - surrogates - number of surrogates to compute MI for, showing the null distribution if there were no
%    relationship at the given lag.
% Outputs
% - miHeartToBreath - MI (heart ; breath) for each value of lags


function [miHeartToBreath] = runHeartBreathRateKraskovMIWithLags(lags, numSurrogates)

	tic;
	
	% Add Octave utilities to the path
	addpath('../../../demos/octave/');

	% Assumes the jar is three levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	if (nargin < 1)
		lags = [0];
	end

	if (nargin < 2)
		numSurrogates = 0;
	end

	data = load('../../../demos/data/SFI-heartRate_breathVol_bloodOx.txt');
	
	% Restrict to the samples that Schreiber mentions:
	data = data(2350:3550,:);
	
	% Separate the data from each column:
	heart = data(:,1);
	chestVol = data(:,2);
	bloodOx = data(:,3);
	timeSteps = length(heart);
	
	fprintf('MI for heart rate <-> breath rate for Kraskov estimation with %d samples:\n', timeSteps);

	% Using a KSG estimator for MI is the least biased way to run this:
	miCalc=javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');

	for lagIndex = 1:length(lags)
		lag = lags(lagIndex);
		% Compute an MI value for this lag
		
		miCalc.initialise(1,1); % univariate calculation
		miCalc.setProperty('k', '4');
		if (lag < 0)
			source = chestVol;
			target = heart;
		else
			source = heart;
			target = chestVol;
		end
		miCalc.setProperty('TIME_DIFF', sprintf("%d", abs(lag)));
		miCalc.setObservations(octaveToJavaDoubleArray(source), ...
					octaveToJavaDoubleArray(target));
		miHeartToBreath(lagIndex) = miCalc.computeAverageLocalOfObservations();
		fprintf('MI(lag=%d): = %.3f nats', lag, miHeartToBreath(lagIndex));
		if (numSurrogates > 0)
			% My observations suggest that most results here for lags -15:15 are statistically significant
			% (threshold for multiple comparisons corrected 0.05 level would be about 0.08 nats.)
			miHeartToBreathNullDist = miCalc.computeSignificance(numSurrogates);
			miHeartToBreathNullMean = miHeartToBreathNullDist.getMeanOfDistribution();
			miHeartToBreathNullStd = miHeartToBreathNullDist.getStdOfDistribution();
			fprintf(" (null = %.3f +/- %.3f)", miHeartToBreathNullMean, miHeartToBreathNullStd),
		end
		fprintf("\n");
	end
		
	tElapsed = toc;
	fprintf('Total runtime was %.1f sec\n', tElapsed);
	
	plot(lags, miHeartToBreath, 'rx', 'markersize', 10);
	set (gca,'fontsize',26);
	xlabel('Lag', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('MI (nats)', 'FontSize', 36, 'FontWeight', 'bold');
	print('heartBreathResults-kraskovMI.eps', '-depsc');
end

