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

% function [teHeartToBreath, teBreathToHeart] = runHeartBreathRateKraskov(kHistory, lHistory, knns, numSurrogates)
%
% runHeartBreathRateKraskov
% Version 1.0
% Joseph Lizier
% 22/1/2014
%
% Used to explore information transfer in the heart rate / breath rate example of Schreiber --
%  but estimates TE using Kraskov-Stoegbauer-Grassberger estimation.
%
%
% Inputs
% - kHistory - destination embedding length
% - lHistory - source embedding length
% - knns - a scalar specifying a single, or vector specifying multiple, value of K nearest neighbours to evaluate TE (Kraskov) with.
% - numSurrogates - a scalar specifying the number of surrogates to evaluate TE from null distribution
% Outputs
% - teHeartToBreath - TE (heart -> breath) for each value of k nearest neighbours
% - teBreathToHeart - TE (breath -> heart) for each value of k nearest neighbours


function [teHeartToBreath, teBreathToHeart] = runHeartBreathRateKraskov(kHistory, lHistory, knns, numSurrogates)

	tic;
	
	% Add utilities to the path
	addpath('..');

	% Assumes the jar is three levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	if (nargin < 4)
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
	
	fprintf('TE for heart rate <-> breath rate for Kraskov estimation with %d samples:\n', timeSteps);

	% Using a KSG estimator for TE is the least biased way to run this:
	teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');

	for knnIndex = 1:length(knns)
		knn = knns(knnIndex);
		% Compute a TE value for knn nearest neighbours
				
		% Perform calculation for heart -> breath (lag 1)
		teCalc.initialise(kHistory,1,lHistory,1,1);
		teCalc.setProperty('k', sprintf('%d',knn));
		teCalc.setObservations(octaveToJavaDoubleArray(heart), ...
					octaveToJavaDoubleArray(chestVol));
		teHeartToBreath(knnIndex) = teCalc.computeAverageLocalOfObservations();
		if (numSurrogates > 0)
			teHeartToBreathNullDist = teCalc.computeSignificance(numSurrogates);
			teHeartToBreathNullMean = teHeartToBreathNullDist.getMeanOfDistribution();
			teHeartToBreathNullStd = teHeartToBreathNullDist.getStdOfDistribution();
		end
		
		% Perform calculation for breath -> heart (lag 1)
		teCalc.initialise(kHistory,1,lHistory,1,1);
		teCalc.setProperty('k', sprintf('%d',knn));
		teCalc.setObservations(octaveToJavaDoubleArray(chestVol), ...
					octaveToJavaDoubleArray(heart));
		teBreathToHeart(knnIndex) = teCalc.computeAverageLocalOfObservations();
		if (numSurrogates > 0)
			teBreathToHeartNullDist = teCalc.computeSignificance(numSurrogates);
			teBreathToHeartNullMean = teBreathToHeartNullDist.getMeanOfDistribution();
			teBreathToHeartNullStd = teBreathToHeartNullDist.getStdOfDistribution();
		end
		
		fprintf('TE(k=%d,l=%d,knn=%d): h->b = %.3f', kHistory, lHistory, knn, teHeartToBreath(knnIndex));
		if (numSurrogates > 0)
			fprintf(' (null = %.3f +/- %.3f)', teHeartToBreathNullMean, teHeartToBreathNullStd);
		end
		fprintf(', b->h = %.3f nats', teBreathToHeart(knnIndex));
		if (numSurrogates > 0)
			fprintf('(null = %.3f +/- %.3f)\n', teBreathToHeartNullMean, teBreathToHeartNullStd);
		else
			fprintf('\n');
		end
	end
		
	tElapsed = toc;
	fprintf('Total runtime was %.1f sec\n', tElapsed);
	
	hold off;
	plot(knns, teHeartToBreath, 'rx', 'markersize', 10);
	hold on;
	plot(knns, teBreathToHeart, 'mo', 'markersize', 10);
	hold off;
	legend(['TE(heart->breath)'; 'TE(breath->heart)']);
	set (gca,'fontsize',26);
	xlabel('K nearest neighbours', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('TE', 'FontSize', 36, 'FontWeight', 'bold');
	print('heartBreathResults-kraskovTE.eps', '-depsc');
end

