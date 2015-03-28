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

% function [aisHeart, aisBreath] = activeInfoStorageHeartBreathRatesKraskov(kHistories, knn, numSurrogates)
%
% activeInfoStorageHeartBreathRatesKraskov
% Version 1.0
% Joseph Lizier
% 04/04/2014
%
% Used to explore active information storage in the heart rate / breath rate example of Schreiber --
%  estimated using Kraskov-Grassberger estimation.
%
% Usually plateaus of AIS indicate that the correct embedding is found; for Kraskov estimation,
%  a peak will indicate this (given that bias correction will pull down values for larger k.
%
%
% Inputs
% - kHistories - a vector of which embedded history lengths to evaluate for both variables
% - knn - a scalar specifying a single value of K nearest neighbours to evaluate AIS (Kraskov) with
% - numSurrogates - a scalar specifying the number of surrogates to evaluate AIS from null distribution
% Outputs
% - aisHeart - active information storage TE (heart -> breath) for each value of k nearest neighbours
% - aisBreath - TE (breath -> heart) for each value of k nearest neighbours


function [aisHeart, aisBreath] = activeInfoStorageHeartBreathRatesKraskov(kHistories, knn, numSurrogates)

	tic;
	
	% Add utilities to the path
	addpath('..');

	% Assumes the jar is two levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	if (nargin < 3)
		numSurrogates = 0;
	end

	data = load('../../data/SFI-heartRate_breathVol_bloodOx.txt');
	
	% Restrict to the samples that Schreiber mentions:
	data = data(2350:3550,:);
	
	% Separate the data from each column:
	heart = data(:,1);
	chestVol = data(:,2);
	bloodOx = data(:,3);
	timeSteps = length(heart);
	
	fprintf('AIS for heart rate and breath rate for Kraskov estimation with %d samples:\n', timeSteps);

	aisCalc=javaObject('infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov');

	for kIndex = 1:length(kHistories)
		kHistory = kHistories(kIndex);
		% Compute an AIS value for this embedding length for each variable:
				
		% Perform calculation for heart
		aisCalc.setProperty('k', sprintf('%d',knn));
		aisCalc.setProperty('NORMALISE', 'true');
		aisCalc.initialise(kHistory); % Use history length kHistory (Schreiber k)
		aisCalc.setObservations(octaveToJavaDoubleArray(heart(1:timeSteps)));
		aisHeart(kIndex) = aisCalc.computeAverageLocalOfObservations();
		if (numSurrogates > 0)
			aisHeartNullDist = aisCalc.computeSignificance(numSurrogates);
			aisHeartNullMean = aisHeartNullDist.getMeanOfDistribution();
			aisHeartNullStd = aisHeartNullDist.getStdOfDistribution();
		end
		
		% Perform calculation for breath
		aisCalc.setProperty('k', sprintf('%d',knn));
		aisCalc.setProperty('NORMALISE', 'true');
		aisCalc.initialise(kHistory); % Use history length kHistory (Schreiber k)
		aisCalc.setObservations(octaveToJavaDoubleArray(chestVol(1:timeSteps)));
		aisBreath(kIndex) = aisCalc.computeAverageLocalOfObservations();
		if (numSurrogates > 0)
			aisBreathNullDist = aisCalc.computeSignificance(numSurrogates);
			aisBreathNullMean = aisBreathNullDist.getMeanOfDistribution();
			aisBreathNullStd = aisBreathNullDist.getStdOfDistribution();
		end
		
		fprintf('AIS(k=%d,knns=%d): h = %.3f',  kHistory, knn, aisHeart(kIndex));
		if (numSurrogates > 0)
			fprintf(' (null = %.3f +/- %.3f)', aisHeartNullMean, aisHeartNullStd);
		end
		fprintf(', b = %.3f', aisBreath(kIndex));
		if (numSurrogates > 0)
			fprintf('(null = %.3f +/- %.3f)\n', aisBreathNullMean, aisBreathNullStd);
		else
			fprintf('\n');
		end
	end
	
	totaltime = toc;
	fprintf('Total runtime was %.1f sec\n', totaltime);
	
	hold off;
	plot(kHistories, aisHeart, 'rx', 'markersize', 10);
	hold on;
	plot(kHistories, aisBreath, 'bo', 'markersize', 10);
	hold off;
	legend(['AIS(Heart) '; 'AIS(Breath)'])
	set (gca,'fontsize',26);
	xlabel('embedding history k', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('AIS(k)', 'FontSize', 36, 'FontWeight', 'bold');
	print('heartBreathResults-kraskovAIS.eps', '-depsc');
end

