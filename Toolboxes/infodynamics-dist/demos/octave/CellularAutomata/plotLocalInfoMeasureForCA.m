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

% function plotLocalInfoMeasureForCA(neighbourhood, base, rule, cells, timeSteps, measureId, measureParams, options)
% 
% Plot one run of the given CA and compute and plot a local information dynamics profile for it
%
% Inputs:
% - neighbourhood - neighbourhood size for the rule (ECA has neighbourhood 3).
%     For an even size neighbourhood (meaning a different number of neighbours on each side of the cell),
%     we take an extra cell from the lower cell indices (i.e. from the left).
% - base - number of discrete states for each cell (for binary states this is 2)
% - rule - supplied as either:
%     a. an integer rule number if <= 2^31 - 1 (Wolfram style; e.g. 110, 54 are the complex ECA rules)
%     b. a HEX string, e.g. phi_par from Mitchell et al. is '0xfeedffdec1aaeec0eef000a0e1a020a0' (note: the leading 0x is not required)
% - cells - number of cells in the CA
% - timeSteps - number of rows to execute the CA for (including the random initial row)
% - measureId - which local info dynamics measure to plot - can be a string or an integer as follows:
%    - 'active', 0 - active information storage (requires measureParams.k)
%    - 'transfer', 1 - pairwise or apparent transfer entropy (requires measureParams.k and j)
%    - 'transfercomplete', 2 - complete transfer entropy (requires measureParams.k and j), conditioning on all other sources
%    - 'separable', 3 - separable information (requires measureParams.k)
%    - 'entropy', 4 - excess entropy (requires no params)
%    - 'entropyrate', 5 - excess entropy (requires measureParams.k)
%    - 'excess', 6 - excess entropy (requires measureParams.k)
%    - 'all', -1 - plot all measures
% - measureParams - a structure containing options as described for each measure above:
%    - measureParams.k - history length for information dynamics measures
%        (for excess entropy, predictive information formulation, this is the history length and future length)
%    - measureParams.j - we measure information transfer across j cells to the right per time step
% - options - a stucture containing a range of other options, i.e.:
%    - plotOptions - structure as defined for the plotLocalInfoValues function
%    - seed - state for the random number generator used to set the initial condition of the CA (use this
%       for reproducibility of plots, or to produce profiles for several different measures of the same CA raw states).
%       We set rand('state', options.seed) if options.seed is supplied, and restore the previous seed afterwards.
%    - initialState - supply the initial state of the CA. If this option exists, no other random initial state is set.
%    - plotRawCa - default true
%    - saveImages - whether to save the plots or not (default false)
%    - saveImagesFormat - 'eps' or 'pdf' (Default eps)
%    - movingFrameSpeed - moving frame of reference's cells/time step, as in Lizier & Mahoney paper (default 0)

function [caStates, localValues] = plotLocalInfoMeasureForCA(neighbourhood, base, rule, cells, timeSteps, measureId, measureParams, options)

	tic

	if (nargin < 8)
		options = {};
	end
	if not(isfield(options, 'plotOptions'))
		options.plotOptions = {}; % Create it ready for plotRawCa etc
	end
	if not(isfield(options, 'saveImages'))
		options.saveImages = false;
	end
	if not(isfield(options, 'saveImagesFormat'))
		options.saveImagesFormat = 'eps';
	end
	if not(isfield(options, 'plotRawCa'))
		options.plotRawCa = true;
	end

	% Add utilities to the path
	addpath('..');

	% Assumes the jar is two levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	% Simulate and plot the CA
	if (isfield(options, 'initialState'))
		% User has supplied an initial state for the CA
		caStates = runCA(neighbourhood, base, rule, cells, timeSteps, false, options.initialState);
	elseif (isfield(options, 'seed'))
		previousSeed = rand('state');
		caStates = runCA(neighbourhood, base, rule, cells, timeSteps, false, options.seed);
		rand('state', previousSeed);
	else
		caStates = runCA(neighbourhood, base, rule, cells, timeSteps, false);
	end
	if (options.plotRawCa)
		plotRawCa(caStates, rule, options.plotOptions, options.saveImages, options.saveImagesFormat);
	end
	if (ischar(rule))
		ruleString = rule;
	else
		ruleString = sprintf('%d', rule);
	end
	if (strcmp(options.saveImagesFormat, 'eps'))
		printDriver = 'epsc'; % to force colour
		fontSize = 32;
	else
		printDriver = options.saveImagesFormat;
		fontSize = 13;
	end
	figNum = 2;
	toc
	% The offsets of the parents (see runCA for how this is computed, especially for even neighbourhood):
	fullSetOfParents = ceil(-neighbourhood / 2) : ceil(-neighbourhood / 2) + (neighbourhood-1);

	if (isfield(options, 'movingFrameSpeed'))
		% User has requested us to evaluate information dynamics with a moving frame of reference
		% (see Lizier and Mahoney paper).
		% The shift for each row of the CA is the negative of the movingFrameSpeed, since
		%  frame moving at 1 cell/time step is same as next row being shifted backwards by 1 cell:
		caStates = accumulateShift(caStates, -options.movingFrameSpeed);
		% Also account for the moved frame of reference in the offsets of parents to the destination:
		fullSetOfParents = fullSetOfParents - options.movingFrameSpeed;
		if (isfield(measureParams, 'j'))
			% Also account for the moved frame of reference in the j parameter for transfer across j cells
			measureParams.j = measureParams.j - options.movingFrameSpeed;
		end
	end

	% convert the states to a format usable by java:
	caStatesJInts = octaveToJavaIntMatrix(caStates);
	toc
	
	plottedOne = false;

	% Make the local information dynamics measurement(s)
	
	%============================
	% Active information storage
	if ((ischar(measureId) && (strcmpi('active', measureId) || strcmpi('all', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 0) || (measureId == -1))))
		% Compute active information storage
		activeCalc = javaObject('infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete', base, measureParams.k);
		activeCalc.initialise();
		activeCalc.addObservations(caStatesJInts);
		avActive = activeCalc.computeAverageLocalOfObservations();
		fprintf('Average active information storage = %.4f\n', avActive);
		javaLocalValues = activeCalc.computeLocalFromPreviousObservations(caStatesJInts);
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states:
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-active-k%d.%s', ruleString, measureParams.k, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end
	
	%============================
	% Apparent transfer entropy
	if ((ischar(measureId) && (strcmpi('transfer', measureId) || strcmpi('all', measureId) || strcmpi('apparenttransfer', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 1) || (measureId == -1))))
		% Compute apparent transfer entropy
		if (measureParams.j == 0)
			error('Cannot compute transfer entropy from a cell to itself (setting measureParams.j == 0)');
		end
		transferCalc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', base, measureParams.k);
		transferCalc.initialise();
		transferCalc.addObservations(caStatesJInts, measureParams.j);
		avTransfer = transferCalc.computeAverageLocalOfObservations();
		fprintf('Average apparent transfer entropy (j=%d) = %.4f\n', measureParams.j, avTransfer);
		javaLocalValues = transferCalc.computeLocalFromPreviousObservations(caStatesJInts, measureParams.j);
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states:
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-transfer-k%d-j%d.%s', ruleString, measureParams.k, measureParams.j, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end
	
	%============================
	% Complete transfer entropy, conditioning on all other sources
	if ((ischar(measureId) && (strcmpi('transfercomplete', measureId) || strcmpi('completetransfer', measureId) || strcmpi('all', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 2) || (measureId == -1))))
		% Compute complete transfer entropy
		if (measureParams.j == 0)
			error('Cannot compute transfer entropy from a cell to itself (setting measureParams.j == 0)');
		end
		transferCalc = javaObject('infodynamics.measures.discrete.ConditionalTransferEntropyCalculatorDiscrete', ...
			base, measureParams.k, neighbourhood - 2);
		transferCalc.initialise();
		% Offsets of all parents can be included here - even 0 and j, these will be eliminated internally:
		transferCalc.addObservations(caStatesJInts, measureParams.j, octaveToJavaIntArray(fullSetOfParents));
		avTransfer = transferCalc.computeAverageLocalOfObservations();
		fprintf('Average complete transfer entropy (j=%d) = %.4f\n', measureParams.j, avTransfer);
		javaLocalValues = transferCalc.computeLocalFromPreviousObservations(caStatesJInts, ...
					measureParams.j, octaveToJavaIntArray(fullSetOfParents));
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states:
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-transferComp-k%d-j%d.%s', ruleString, measureParams.k, measureParams.j, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end
	
	%============================
	% Separable information
	if ((ischar(measureId) && (strcmpi('separable', measureId) || strcmpi('all', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 3) || (measureId == -1))))
		% Compute separable information
		separableCalc = javaObject('infodynamics.measures.discrete.SeparableInfoCalculatorDiscrete', ...
			base, measureParams.k, neighbourhood - 1);
		separableCalc.initialise();
		% Offsets of all parents can be included here - even 0 and j, these will be eliminated internally:
		separableCalc.addObservations(caStatesJInts, octaveToJavaIntArray(fullSetOfParents));
		avSeparable = separableCalc.computeAverageLocalOfObservations();
		fprintf('Average separable information = %.4f\n', avSeparable);
		javaLocalValues = separableCalc.computeLocalFromPreviousObservations(caStatesJInts, ...
					octaveToJavaIntArray(fullSetOfParents));
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states:
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-separable-k%d.%s', ruleString, measureParams.k, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end

	%============================
	% Entropy
	if ((ischar(measureId) && (strcmpi('entropy', measureId) || strcmpi('all', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 4) || (measureId == -1))))
		% Compute entropy
		entropyCalc = javaObject('infodynamics.measures.discrete.EntropyCalculatorDiscrete', ...
			base);
		entropyCalc.initialise();
		entropyCalc.addObservations(caStatesJInts);
		avEntropy = entropyCalc.computeAverageLocalOfObservations();
		fprintf('Average entropy = %.4f\n', avEntropy);
		javaLocalValues = entropyCalc.computeLocalFromPreviousObservations(caStatesJInts);
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states):
			% (Note for entropy, the shifts to and back don't make any difference)
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-entropy.%s', ruleString, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end

	%============================
	% Entropy rate
	if ((ischar(measureId) && (strcmpi('entropyrate', measureId) || strcmpi('all', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 5) || (measureId == -1))))
		% Compute entropy rate
		entRateCalc = javaObject('infodynamics.measures.discrete.EntropyRateCalculatorDiscrete', base, measureParams.k);
		entRateCalc.initialise();
		entRateCalc.addObservations(caStatesJInts);
		avEntRate = entRateCalc.computeAverageLocalOfObservations();
		fprintf('Average entropy rate = %.4f\n', avEntRate);
		javaLocalValues = entRateCalc.computeLocalFromPreviousObservations(caStatesJInts);
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states:
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-entrate-k%d.%s', ruleString, measureParams.k, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end
	
	%============================
	% Excess entropy
	if ((ischar(measureId) && (strcmpi('excess', measureId) || strcmpi('all', measureId))) || ...
	    (not(ischar(measureId)) && ((measureId == 6) || (measureId == -1))))
		% Compute excess entropy
		excessEntropyCalc = javaObject('infodynamics.measures.discrete.PredictiveInformationCalculatorDiscrete', base, measureParams.k);
		excessEntropyCalc.initialise();
		excessEntropyCalc.addObservations(caStatesJInts);
		avExcessEnt = excessEntropyCalc.computeAverageLocalOfObservations();
		fprintf('Average excess entropy = %.4f\n', avExcessEnt);
		javaLocalValues = excessEntropyCalc.computeLocalFromPreviousObservations(caStatesJInts);
		localValues = javaMatrixToOctave(javaLocalValues);
		if (isfield(options, 'movingFrameSpeed'))
			% User has requested us to evaluate information dynamics with a moving frame of reference
			% (see Lizier and Mahoney paper).
			% Need to shift the computed info dynamics back (to compensate for earlier shift to CA states:
			localValues = accumulateShift(localValues, options.movingFrameSpeed);
		end
		toc
		figure(figNum)
		figNum = figNum + 1;
		plotLocalInfoValues(localValues, options.plotOptions);
		if (options.saveImages)
			set(gca, 'fontsize', fontSize);
			colorbar('fontsize', fontSize);
			print(sprintf('figures/%s-excessentropy-k%d.%s', ruleString, measureParams.k, options.saveImagesFormat), sprintf('-d%s', printDriver));
		end
		plottedOne = true;
	end

	if (not(plottedOne))
		error(sprintf('Supplied measureId %s did not match any measurement types', measureId));
	end

	toc
end

% Perform a circular shift of each row of the matrix, shifting the first row by zero,
% the second by shiftPerRow, the third by 2*shiftPerRow, and so on
function shiftedMatrix = accumulateShift(matrix, shiftPerRow)
	% Allocate required space up front (makes this much faster):
	shiftedMatrix = zeros(size(matrix,1), size(matrix, 2));
	shiftedMatrix(1,:) = matrix(1,:);
	for r = 2 : size(matrix,1)
		% circshift operates on shifting rows, so we transpose the input and output to it:
		shiftedMatrix(r,:) = circshift(matrix(r,:)', shiftPerRow*(r-1))'; 
	end
end

