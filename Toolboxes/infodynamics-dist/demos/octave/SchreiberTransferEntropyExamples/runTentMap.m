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

% function Values = runTentMap()
%
% Version 1.0
% Joseph Lizier
% 1/8/14
%
% Used to explore information transfer in the testMap example of Schreiber's paper,
%  recreating figure 1 in that paper.
% 
% The code should take 15-30 minutes run time for all repeats and couplings.

function teValues = runTentMap()

	% Add utilities to the path (needed for octaveToJavaIntMatrix)
	addpath('..');

	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	% Number of cells (M), transients, iterates and number of runs
	%  used by Schreiber:
	M = 100;
	transientLength = 100000;
	iterates = 100000;
	numberOfRuns = 10;
	couplings = 0:0.002:0.05;

	tic;
	
	% Construct for binary values, k=1 history length
	teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 1);
	
	teValues = zeros(numberOfRuns, length(couplings));
	for couplingIndex = 1:length(couplings)
		coupling = couplings(couplingIndex);
		for r = 1:numberOfRuns
			% Construct the map matrix with maps in columns, time in rows

			% Initialise first row randomly
			transientMapValues = rand(1, M);
			% Run transients - no need to keep the transient values
			for t = 2 : transientLength
				transientMapValues = tentMap(coupling .* [ transientMapValues(M) , transientMapValues(1:M-1)] + ...
							     (1 - coupling) .* transientMapValues);
			end
	
			% Run iterates - now keep the iterated map values
			mapValues = zeros(iterates, M);
			mapValues(1,:) = transientMapValues;
			for t = 2 : iterates
				mapValues(t, :) = tentMap(coupling .* [ mapValues(t-1, M) , mapValues(t-1,1:M-1)] + ...
							  (1 - coupling) .* mapValues(t-1,:));
			end
		
			% Take binary values (need the *1 for now for Java conversion from boolean)
			binaryValues = (mapValues >= 0.5)*1;
		
			% Compute TE
			teCalc.initialise();
			% Add observations for TE across 1 column per time step:
			teCalc.addObservations(octaveToJavaIntMatrix(binaryValues), 1);
			teValues(r, couplingIndex) = teCalc.computeAverageLocalOfObservations();
			fprintf('teValues(r=%d, coupling=%.3f)=%.4f\n', r, coupling, teValues(r, couplingIndex));
		end
	end

	meanTes = mean(teValues);
	stdTes = std(teValues);

	% Plot results 
	hold off;
	h = errorbar(couplings, meanTes, stdTes./sqrt(numberOfRuns));
	set(h, 'LineStyle', 'none');
	set(h, 'MarkerEdgeColor', 'r');
	set(h, 'Marker', 'x');
	set(h, 'markersize', 14);
	hold on;
	% Add the curve that Schreiber fitted:
	plot(couplings, 0.77^2*couplings.^2 ./ log(2), '-r');
	hold off;
	set (gca,'fontsize',26);
	xlabel('coupling', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('TE (bits)', 'FontSize', 36, 'FontWeight', 'bold');
	axis([0 0.05 0 0.0025]);
	print('tentMapResults.eps', '-depsc');

	totaltime = toc;
	fprintf('Total runtime was %.1f sec\n', totaltime);
end

function xOut = tentMap(xIn)
	xOut = xIn .* 2 .* (xIn < 0.5) + (2 - 2 .* xIn) .* (xIn >= 0.5);
end
