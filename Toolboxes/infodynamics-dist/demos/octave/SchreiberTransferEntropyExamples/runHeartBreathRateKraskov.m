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
% - kHistory - destination embedding length, or "auto" for an auto-embedding (Ragwitz criteria, which takes a minute or two to run)
% - lHistory - source embedding length, or "auto" for an auto-embedding (Ragwitz criteria, which takes a minute or two to run)
% - knns - a scalar specifying a single, or vector specifying multiple, value of K nearest neighbours to evaluate TE (Kraskov) with.
% - numSurrogates - a scalar specifying the number of surrogates to evaluate TE from null distribution (which further multiplies the runtime)
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
	% Set up for any potential auto-embedding:
	if (ischar(kHistory)) % Assume == 'auto'
		% we're auto-embedding at least the destination:
		if (ischar(lHistory)) % Assume == 'auto'
			% we're auto embedding both source and destination
			teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, ...
		                teCalc.AUTO_EMBED_METHOD_RAGWITZ);
		else
			% we're auto embedding destination only
			teCalc.setProperty(teCalc.PROP_AUTO_EMBED_METHOD, ...
		                teCalc.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY);
		end
		teCalc.setProperty(teCalc.PROP_K_SEARCH_MAX, '10');
		teCalc.setProperty(teCalc.PROP_TAU_SEARCH_MAX, '5');
	end

	for knnIndex = 1:length(knns)
		knn = knns(knnIndex);
		% Compute a TE value for knn nearest neighbours
				
		% Perform calculation for heart -> breath (lag 1)
		if (ischar(kHistory)) % Assume == 'auto'
			% we're auto-embedding at least the destination:
			teCalc.initialise();
		else
			% We're not auto embedding
			teCalc.initialise(kHistory,1,lHistory,1,1);
		end
		teCalc.setProperty('k', sprintf('%d',knn));
		teCalc.setObservations(octaveToJavaDoubleArray(heart), ...
					octaveToJavaDoubleArray(chestVol));
		teHeartToBreath(knnIndex) = teCalc.computeAverageLocalOfObservations();
		% Grab the embedding parameters (in case of auto-embedding), converting from Java to native strings
		kUsedHB = char(teCalc.getProperty(teCalc.K_PROP_NAME));
		kTauUsedHB = char(teCalc.getProperty(teCalc.K_TAU_PROP_NAME));
		lUsedHB = char(teCalc.getProperty(teCalc.L_PROP_NAME));
		lTauUsedHB = char(teCalc.getProperty(teCalc.L_TAU_PROP_NAME));
		% And compare to surrogates if required
		if (numSurrogates > 0)
			teHeartToBreathNullDist = teCalc.computeSignificance(numSurrogates);
			teHeartToBreathNullMean = teHeartToBreathNullDist.getMeanOfDistribution();
			teHeartToBreathNullStd = teHeartToBreathNullDist.getStdOfDistribution();
		end
		
		% Perform calculation for breath -> heart (lag 1)
		if (ischar(kHistory)) % Assume == 'auto'
			% we're auto-embedding at least the destination:
			teCalc.initialise();
		else
			% We're not auto embedding
			teCalc.initialise(kHistory,1,lHistory,1,1);
		end
		teCalc.setProperty('k', sprintf('%d',knn));
		teCalc.setObservations(octaveToJavaDoubleArray(chestVol), ...
					octaveToJavaDoubleArray(heart));
		teBreathToHeart(knnIndex) = teCalc.computeAverageLocalOfObservations();
		% Grab the embedding parameters (in case of auto-embedding)
		kUsedBH = char(teCalc.getProperty(teCalc.K_PROP_NAME));
		kTauUsedBH = char(teCalc.getProperty(teCalc.K_TAU_PROP_NAME));
		lUsedBH = char(teCalc.getProperty(teCalc.L_PROP_NAME));
		lTauUsedBH = char(teCalc.getProperty(teCalc.L_TAU_PROP_NAME));
		% And compare to surrogates if required
		if (numSurrogates > 0)
			teBreathToHeartNullDist = teCalc.computeSignificance(numSurrogates);
			teBreathToHeartNullMean = teBreathToHeartNullDist.getMeanOfDistribution();
			teBreathToHeartNullStd = teBreathToHeartNullDist.getStdOfDistribution();
		end
		
		fprintf('TE(k=%s,kTau=%s,l=%s,lTau=%s,knn=%d): h->b = %.3f', kUsedHB, kTauUsedHB, lUsedHB, lTauUsedHB, knn, teHeartToBreath(knnIndex));
		if (numSurrogates > 0)
			fprintf(' (null = %.3f +/- %.3f)', teHeartToBreathNullMean, teHeartToBreathNullStd);
		end
		fprintf('; TE(k=%s,kTau=%s,l=%s,lTau=%s,knn=%d): b->h = %.3f nats', kUsedBH, kTauUsedBH, lUsedBH, lTauUsedBH, knn, teBreathToHeart(knnIndex));
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

