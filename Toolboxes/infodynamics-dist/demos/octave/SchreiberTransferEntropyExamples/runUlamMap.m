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

% function [teValues1to2, teValues2to1] = runUlamMap()
%
% Version 1.0
% Joseph Lizier
% 1/8/14
%
% Used to explore information transfer in the Ulam Map example of Schreiber's paper,
%  recreating figure 2 in that paper.
% 
% The code should take 2-5 minutes run time for all repeats and couplings.

function [teValues1to2, teValues2to1] = runUlamMap()

	% Add utilities to the path (needed for octaveToJavaIntMatrix)
	addpath('..');

	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	% Number of cells (M), transients, iterates and number of runs
	%  used by Schreiber:
	M = 100;
	transientLength = 100000;
	iterates = 10000;
	couplings = 0:0.02:1; % It seems Schreiber used a step size of 0.02
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Schreiber's paper states a kernel width of 0.2 was used, but I'm fairly sure
	%  this should actually be 0.3
	kernelWidth = 0.3;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	tic;

	% Construct for kernel estimation, k=1 history length
	teCalc=javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
	teCalc.setProperty('NORMALISE', 'false'); % Normalise the individual variables. Not mentioned in the paper.
	teCalc.setProperty('DYN_CORR_EXCL', '100'); % Dynamic correlation exclusion time - Schreiber uses 100 steps
	
	teValues1to2 = zeros(1, length(couplings));
	teValues2to1 = zeros(1, length(couplings));

	for couplingIndex = 1:length(couplings)
		coupling = couplings(couplingIndex);

		% Initialise first row randomly (values should be in range [-2,2])
		transientMapValues = -2 + rand(1, M) .* 4;
		% Run transients - no need to keep the transient values
		for t = 2 : transientLength
			transientMapValues = ulamMap(coupling .* [ transientMapValues(M) , transientMapValues(1:M-1)] + ...
						     (1 - coupling) .* transientMapValues);
		end

		% Run iterates - now keep the iterated map values
		mapValues = zeros(iterates, M);
		mapValues(1,:) = transientMapValues;
		for t = 2 : iterates
			mapValues(t, :) = ulamMap(coupling .* [ mapValues(t-1, M) , mapValues(t-1,1:M-1)] + ...
						  (1 - coupling) .* mapValues(t-1,:));
		end
		
		% Take only the first two columns, as Schreiber does
		x1 = mapValues(:,1);
		x2 = mapValues(:,2);
		
		% Compute TE(1 -> 2)
		teCalc.initialise(1, kernelWidth); % Use history length k=1, kernel width supplied
		teCalc.setObservations(x1, x2);
		teValues1to2(couplingIndex) = teCalc.computeAverageLocalOfObservations();
		% One could attempt bias correction, though this makes no difference for 10k points here:
		% teValues1to2(couplingIndex) = teCalc.computeAverageLocalOfObservationsWithCorrection();

		% Compute TE(2 -> 1)
		teCalc.initialise(1, kernelWidth); % Use history length k=1, kernel width supplied
		teCalc.setObservations(x2, x1);
		teValues2to1(couplingIndex) = teCalc.computeAverageLocalOfObservations();
		% One could attempt bias correction, though this makes no difference for 10k points here:
		%teValues2to1(couplingIndex) = teCalc.computeAverageLocalOfObservationsWithCorrection();
		fprintf('coupling=%.3f: te1to2=%.4f, te2to1=%.4f\n', coupling, teValues1to2(couplingIndex), teValues2to1(couplingIndex));
	end
	
	% Plot results 
	hold off;
	plot(couplings, teValues1to2, '-r');
	hold on;
	plot(couplings, teValues2to1, '-g');
	% And add points we extracted from Schreiber's plot for TE(1->2):
	schreiberResults1to2 =  load('SchreiberExample2.txt');
	plot(schreiberResults1to2(:,1), schreiberResults1to2(:,2), '-b');
	hold off;
	set (gca,'fontsize',26);
	xlabel('coupling', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('TE (bits)', 'FontSize', 36, 'FontWeight', 'bold');
	axis([0 1 -0.2 2.2]);
	h = legend('TE(1->2)', 'TE(2->1)', 'Schreiber(1->2)', 'Location', 'NorthWest');
	set (h,'fontsize',12);
	print('ulamMapResults.eps', '-depsc');

	totaltime = toc;
	fprintf('Total runtime was %.1f sec\n', totaltime);
end

function xOut = ulamMap(xIn)
	xOut = 2 - xIn .^ 2;
end

