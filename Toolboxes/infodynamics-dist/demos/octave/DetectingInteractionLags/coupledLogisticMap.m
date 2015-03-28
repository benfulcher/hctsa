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

% Examine transfer entropy from variable X to Y for coupled logistic maps X and Y
%  where the process is as defined in Example V.A of:
%  Pompe and Runge, Momentary information transfer as a coupling measure of time series,
%   PHYSICAL REVIEW E 83, 051122 (2011) :
%       fmod1(x) := 4 .* (x mod 1) .* (1 - (x mod 1));
% 	gYToX = cYToX .* Y(t - delayYtoX) + (1-cYToX) .* X(t - 1);
%	X(t) = fmod1(gYToX);
%	gXToY = cXToY .* X(t - delayXtoY) + (1-cXToY) .* Y(t - 1);
%	Y(t) = fmod1(gXToY);
%
% Specifically, we use delayYtoX = 5; and delayXtoY = 2; and check
% TE(X->Y, delay 1) and TE(X ->Y, delay 2), because Pompe and Runge reported
% that TE(X->Y, delay 1) was larger than TE(X ->Y, delay 2) despite the X -> Y
% coupling delay being 2.
% Their result appears to be a function of their measurement technique however
% (examing symbolic entropies, or the ordinal relationships between the variables).
% The main clue to this is that for their measure "momentary information transfer"
% (which conditions on the past of the source, in addition to the past of the
% destination), there is significantly non-zero information added by the source
% with lag 1, conditioned on the past of the source and the destination.
% Intuitively, their measure should be zero for lag 1, since that source variable
% should not be able to add any information about the destination (which is
% completely determined by its past and the source at lag 2).
% The source at lag 1 is only able to add information under their measure because
% it is estimated using ordinal values, and focussing on ordinals in the coupled
% logistic map process leaves uncertainty regarding the ordinal position of the
% next state of the destination given the ordinal relationships amongst its past
% and the source at lag 2. Some information regarding this uncertainty seems
% to be provided by the far past of the destination which is directly coupled to
% the source at lag 1 via a lag 5 relationship. This information is useful because
% of the memory in the destination variable.
% In a similar fashion, the transfer entropy at lag 1 contains much information 
% about what the source causally added at lag 2 (because of strong memory in the
% source). But it also contains more information about the far past of Y (beyond 
% the history or embedding length 1 used for the measure) which can be helpful
% to decode the process Y. If one is using ordinal relationships, it appears 
% that this information is more strongly contained in the source at last 1.
%
% However, since we know that Y(t) is causally determined from X(t-2) and Y(t-1)
% then we know that a proper information theoretic characterisation should
% find more information in X(t-2) about Y(t) given Y(t-1) than from X(t-1).
% So here we make such a proper characterisation, using Kraskov estimation
% rather than ordinals, which lose much relevant information about the processes.
% 

function coupledLogisticMap()

	% Constants as set by Pompe and Runge:
	delayYtoX = 5;
	delayXtoY = 2;
	cYToX = 0.2;
	cXToY = 0.5;
	T = 512;
	fprintf('For 1000 repeats, expect the calculations to take ~30 seconds ...\n');
	repeats = 1000; % General results visible for 100 repeats if you want to see them faster (~20 sec)
	k = 1; % history length
	
	% How many steps to randomly seed:
	seedSteps = max([k,delayYtoX,delayXtoY]);
	timeToKeep = T + seedSteps;
	
	% Transients before capturing observations -
	% Pompe and Runge don't state what they use;
	%  by playing with it, changing from 100 to 1000 has no discernable effect so we think this is long enough
	transientRunsOfT = 100;

	% The actual values of the results slightly change with the Kraskov K parameter,
	%  but we note that delay 2 is higher for K's we've tested
	KraskovK ='4'; % Use Kraskov parameter K=4 for 4 nearest points
	% For K=2 we get:
	%TE(X->Y,delay=1) = 0.7773 nats (+/- std 0.0959, stderr 0.0096)
	%TE(X->Y,delay=2) = 1.8759 nats (+/- std 0.0332, stderr 0.0033)

	
	tic;

	% Add utilities to the path
	addpath('..');

	% Assumes the jar is two levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	X = zeros(timeToKeep, repeats);
	X(1:seedSteps,:) = rand(seedSteps,repeats);
	Y = zeros(timeToKeep, repeats);
	Y(1:seedSteps,:) = rand(seedSteps,repeats);
	
	% Should run for some set up states to make sure we've reached stationarity:
	for ts = 1 : transientRunsOfT
		% Run the process for T time steps
		if (ts > 1)
			% We've just been running transients in X and Y, so copy the
			%  last seedSteps up into the first few rows:
			X(1:seedSteps,:) = X(size(X,1)-seedSteps+1:size(X,1),:);
			Y(1:seedSteps,:) = Y(size(Y,1)-seedSteps+1:size(Y,1),:);
		end
		% Now run the process from these initial states
		for n = 1 : T
			t = seedSteps + n;
			gYToX = cYToX .* Y(t - delayYtoX,:) + (1-cYToX) .* X(t - 1,:);
			X(t,:) = fmod1(gYToX);
			gXToY = cXToY .* X(t - delayXtoY,:) + (1-cXToY) .* Y(t - 1,:);
			Y(t,:) = fmod1(gXToY);
		end
	end
	
	% And compute the results on the last T steps:
	
	resultsLag1 = zeros(1,repeats);
	resultsLag2 = zeros(1,repeats);
	resultsLag3 = zeros(1,repeats);
	for r = 1 : repeats
		% Create a TE calculator and run it:
		%  (Our TE calculator is now using a single conditional MI calculator, so
		%   we replace the conditional MI calculator that was here before)
		teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
		% Perform calculation for X -> Y (lag 1)
		teCalc.initialise(k,1,1,1,1); % Use history length k (Schreiber k)
		teCalc.setProperty('k', KraskovK);
		teCalc.setObservations(octaveToJavaDoubleArray(X(seedSteps:size(X,1),r)), ...
					octaveToJavaDoubleArray(Y(seedSteps:size(Y,1),r)));
		resultsLag1(r) = teCalc.computeAverageLocalOfObservations();
		% Perform calculation for X -> Y (lag 2)
		teCalc.initialise(k,1,1,1,2); % Use history length k (Schreiber k)
		teCalc.setProperty('k', KraskovK);
		teCalc.setObservations(octaveToJavaDoubleArray(X(seedSteps:size(X,1),r)), ...
					octaveToJavaDoubleArray(Y(seedSteps:size(Y,1),r)));
		resultsLag2(r) = teCalc.computeAverageLocalOfObservations();
		% Perform calculation for X -> Y (lag 3)
		teCalc.initialise(k,1,1,1,3); % Use history length k (Schreiber k)
		teCalc.setProperty('k', KraskovK);
		teCalc.setObservations(octaveToJavaDoubleArray(X(seedSteps:size(X,1),r)), ...
					octaveToJavaDoubleArray(Y(seedSteps:size(Y,1),r)));
		resultsLag3(r) = teCalc.computeAverageLocalOfObservations();
		
		% Kernel estimator returns the correct ordering of lag 1 and 2 for 
		% reasonably tight values of the kernel width (<~ 0.45 normalised units)
		% At larger kernel widths, the ordering becomes incorrect - it seems that
		%  because larger kernel widths mean we have imprecision in observations that
		%  we are grouping together, then we start to see the same effect as we did with
		%  the symbolic coding!!
		%teCalc=javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
		%kernelWidth = '0.25'; % normalised units
		%% Perform calculation for X -> Y (lag 1)
		%teCalc.initialise(k); % Use history length k (Schreiber k)
		%teCalc.setProperty('EPSILON', kernelWidth);
		%teCalc.setObservations(X(seedSteps:size(X,1),r), Y(seedSteps:size(Y,1),r));
		%resultsLag1(r) = teCalc.computeAverageLocalOfObservations();
		%% Perform calculation for X -> Y (lag 2)
		%teCalc.initialise(k); % Use history length k (Schreiber k)
		%teCalc.setProperty('EPSILON', kernelWidth);
		%teCalc.setObservations(X(seedSteps-1:size(X,1)-1,r), Y(seedSteps:size(Y,1),r));
		%resultsLag2(r) = teCalc.computeAverageLocalOfObservations();
	end
	fprintf('TE(X->Y,delay=1) = %.4f nats (+/- std %.4f, stderr %.4f) or %.4f bits\n', ...
		mean(resultsLag1), std(resultsLag1), std(resultsLag1)./sqrt(repeats), mean(resultsLag1)./log(2) );
	fprintf('TE(X->Y,delay=2) = %.4f nats (+/- std %.4f, stderr %.4f) or %.4f bits\n', ...
		mean(resultsLag2), std(resultsLag2), std(resultsLag2)./sqrt(repeats), mean(resultsLag2) ./log(2));
	fprintf('TE(X->Y,delay=3) = %.4f nats (+/- std %.4f, stderr %.4f) or %.4f bits\n', ...
		mean(resultsLag3), std(resultsLag3), std(resultsLag3)./sqrt(repeats), mean(resultsLag3) ./log(2));
	
	toc;
end


function r = fmod1(x)
	x = mod(x,1);
	r = 4 .* x .* (1-x);
end

