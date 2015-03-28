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

% Compute the transfer entropy (te) and the Pompe-Runge Momentary information transfer (mit)
%  in an example where we have short-term source memory (with decay).
%
% The model is:
%  X_{n-1} + noise -> X_{n}   (noise represents new information entering the source)
%  X_{n} -> Y_{n+1}           (Y_{n+1} copies directly from X_{n} though is correlated to X_{n-1})
% We expect to find the connection X_{n} -> Y_{n+1} dominating, since it is
%  causal, if our measure of transfer is correct. Indeed, the proof of Wibral
%  et al. states that transfer entropy will be maximised here. 
%  We show that the Pompe-Runge Momentary information transfer is actually
%  maximised for lag 2: X_{n-1} -> Y_{n+1}, for many noise levels instead
%  (because for X_{n} -> Y_{n+1} it conditions out the memory from X_{n-1}),
%  showing that it is not maximised
%  for the correct delay even in such simple unidirectional coupling.
%
% This script takes 1-2 minutes to run and generate the plots.
% 
% NOTE: You may need to increase the Java heap space in Matlab for this
%  to work (you will get a "java.lang.OutOfMemoryError: Java heap space"
%  error if this is a problem).
%
% Inputs:
% - savePlot - true if you want eps files of the plots saved
%
function transferWithSourceMemory(savePlot)

	tic;

	% Add utilities to the path
	addpath('..');

	% Assumes the jar is two levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	if (nargin < 1)
		savePlot = false;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This section is here for description only:
	% Possible states of source X
	xValues =        [0, 1, 2, 3];
	% X is self-mapped to its own next state stochastically;
	% - with probability 1-delta, one of the following occurs
	%   with a 50-50 probability (ordinary mapping):
	xMapping1 =      [0, 0, 1, 1];
	xMapping2 =      [2, 2, 3, 3];
	% - or with probability delta, one of the following occurs
	%   with a 50-50 probability (noisy mapping):
	%  (note: noisy mapping removes the memory component)
	xNoisyMapping1 = [1, 1, 0, 0];
	xNoisyMapping2 = [3, 3, 2, 2];
	% Y is mapped from X as:
	yMapping =       [0, 1, 0, 1];
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Empirical results
	N = 1000000; % number of observations
	% noise level (i.e. percent times that X_{n} is altered away from X_{n-1}
	deltas = 0:0.01:0.5;
	teXnToYnplus1 = zeros(1,length(deltas));
	teXnminus1ToYnplus1 = zeros(1,length(deltas));
	mitXnToYnplus1 = zeros(1,length(deltas));
	mitXnminus1ToYnplus1 = zeros(1,length(deltas));
	for deltaIndex = 1:length(deltas)
		delta = deltas(deltaIndex);
		
		% We can implement the (Y,X) process much faster by considering
		%  X as a joint variable (X2,X1) where X2 is updated at random at
		%  each time step, and X1 copies the previous value of X2 with
		%  probability delta, or inverts it with probability 1-delta:
		x2 = rand(N+1,1)<0.5;
		noisyMapping = rand(N,1)<delta;
		x1 = [rand()<0.5; x2(1:N) .* (1 - noisyMapping) + not(x2(1:N)) .* noisyMapping];
		x = x2*2+x1;
		y = [rand() < 0.5; x1(1:N)];
		
		% Compute TEs
		teCalc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 4, 1);
		teCalc.initialise();
		teCalc.addObservations(x, y);
		teXnToYnplus1(deltaIndex) = teCalc.computeAverageLocalOfObservations();
		teCalc.initialise();
		teCalc.addObservations(x(1:length(x)-1), y(2:length(y)));
		teXnminus1ToYnplus1(deltaIndex) = teCalc.computeAverageLocalOfObservations();
		
		% Compute MITs using a conditional TE calculator, adding the past of the source to the conditionals
		compTeCalc = javaObject('infodynamics.measures.discrete.ConditionalTransferEntropyCalculatorDiscrete', 4, 1, 1);
		compTeCalc.initialise();
		% We need to additionally condition on the past of X:
		compTeCalc.addObservations(octaveToJavaIntArray(x(2:length(x))), ...
					octaveToJavaIntArray(y(2:length(y))), ...
					octaveToJavaIntArray(x(1:length(x)-1)));
		mitXnToYnplus1(deltaIndex) = compTeCalc.computeAverageLocalOfObservations();
		compTeCalc.initialise();
		% We need to shift x forward here to investigate the lag, and condition on the past of Y
		compTeCalc.addObservations(octaveToJavaIntArray(x(2:length(x)-1)), ...
					octaveToJavaIntArray(y(3:length(y))), ...
					octaveToJavaIntArray(x(1:length(x)-2)));
		mitXnminus1ToYnplus1(deltaIndex) = compTeCalc.computeAverageLocalOfObservations();
		
		fprintf('delta=%.3f, TE(X_{n} -> Y_{n+1})=%.3f, TE(X_{n-1} -> Y_{n+1})=%.3f, MIT(X_{n} -> Y_{n+1})=%.3f, MIT(X_{n-1} -> Y_{n+1})=%.3f\n', ...
			delta, teXnToYnplus1(deltaIndex), teXnminus1ToYnplus1(deltaIndex), ...
			mitXnToYnplus1(deltaIndex), mitXnminus1ToYnplus1(deltaIndex));
	end	

	% Analytic results:
	teXnToYnplus1_an = ones(1, length(deltas));
	teXnminus1ToYnplus1_an = 1 + (1 - deltas).*log2(1 - deltas) + ...
				deltas .* log2(deltas);
	mitXnToYnplus1_an = 1 - teXnminus1ToYnplus1_an;
	mitXnminus1ToYnplus1_an = teXnminus1ToYnplus1_an;

	% Plot the TEs
	% plot types: h - open diamonds, d - closed diamonds, s - closed squares, p - open squares
	markersize = 15;
	figure;
	if (exist ('OCTAVE_VERSION', 'builtin'))
		% Make the full plots only on octave
		if (savePlot)
			set(gca, 'fontsize', 32); % do this first to get fontsize right for the key
		end
		plot(deltas, teXnToYnplus1, 'x1;Empirical: TE_{SPO}(X \rightarrow Y, 1);', 'markersize', markersize);
		hold on;
		plot(deltas, teXnToYnplus1_an, 'p1;Analytic: TE_{SPO}(X \rightarrow Y, 1);', 'markersize', markersize); % h plots open diamonds
		plot(deltas, teXnminus1ToYnplus1, '+2;Empirical: TE_{SPO}(X \rightarrow y, 2);', 'markersize', markersize);
		plot(deltas, teXnminus1ToYnplus1_an, 'o2;Analytic: TE_{SPO}(X \rightarrow Y, 2);', 'markersize', markersize);
		hold off;
		% Make figures ok for plotting:
		legend('location', 'east');
	else
		% We're on Matlab - just make a quick plot
		plot(deltas, teXnToYnplus1, 'rx', 'markersize', markersize); % Empirical: TE_{SPO}(X \rightarrow Y, 1);
		hold on;
		plot(deltas, teXnToYnplus1_an, 'rp', 'markersize', markersize); % Analytic: TE_{SPO}(X \rightarrow Y, 1)
		plot(deltas, teXnminus1ToYnplus1, 'g+', 'markersize', markersize); % Empirical: TE_{SPO}(X \rightarrow y, 2)
		plot(deltas, teXnminus1ToYnplus1_an, 'go', 'markersize', markersize); % Analytic: TE_{SPO}(X \rightarrow Y, 2)
		hold off;
		% Make figures ok for plotting:
		legend('Empirical: TE_{SPO}(X \rightarrow Y, 1)', 'Analytic: TE_{SPO}(X \rightarrow Y, 1)', ...
			'Empirical: TE_{SPO}(X \rightarrow y, 2)', 'Analytic: TE_{SPO}(X \rightarrow Y, 2)', ...
			'Location', 'East');
	end
	if (savePlot)
		xlabel('\eta', 'fontsize', 32);
		ylabel('Information (bits)', 'fontsize', 32);
		print('te.eps', '-depsc');
	else
		% do these without fontsize to have more stable displays in octave
		xlabel('\eta');
		ylabel('Information (bits)');
	end

	% Plot the MIT's
	figure;
	if (exist ('OCTAVE_VERSION', 'builtin'))
		% Make the full plots only on octave
		if (savePlot)
			set(gca, 'fontsize', 32); % do this first to get fontsize right for the key
		end
		plot(deltas, mitXnToYnplus1, 'x1;Empirical: MIT(X \rightarrow Y, 1);', 'markersize', markersize);
		hold on;
		plot(deltas, mitXnToYnplus1_an, 'p1;Analytic: MIT(X \rightarrow Y, 1);', 'markersize', markersize);
		plot(deltas, mitXnminus1ToYnplus1, '+2;Empirical: MIT(X \rightarrow Y, 2);', 'markersize', markersize);
		plot(deltas, mitXnminus1ToYnplus1_an, 'o2;Analytic: MIT(X \rightarrow Y, 2);', 'markersize', markersize);
		hold off;
		% Make figures ok for plotting:
		legend('location', 'east')
	else
		% We're on Matlab - just make a quick plot
		plot(deltas, mitXnToYnplus1, 'rx', 'markersize', markersize); % Empirical: MIT(X \rightarrow Y, 1)
		hold on;
		plot(deltas, mitXnToYnplus1_an, 'rp', 'markersize', markersize); % Analytic: MIT(X \rightarrow Y, 1)
		plot(deltas, mitXnminus1ToYnplus1, 'g+', 'markersize', markersize); % Empirical: MIT(X \rightarrow Y, 2)
		plot(deltas, mitXnminus1ToYnplus1_an, 'go', 'markersize', markersize); % Analytic: MIT(X \rightarrow Y, 2)
		hold off;
		% Make figures ok for plotting:
		legend('Empirical: MIT(X \rightarrow Y, 1)', 'Analytic: MIT(X \rightarrow Y, 1)', ...
			'Empirical: MIT(X \rightarrow Y, 2)', 'Analytic: MIT(X \rightarrow Y, 2)', ...
			'Location', 'East');
	end
	if (savePlot)
		xlabel('\eta', 'fontsize', 32);
		ylabel('Information (bits)', 'fontsize', 32);
		print('mit.eps', '-depsc');
	else
		% do these without fontsize to have more stable displays in octave
		xlabel('\eta');
		ylabel('Information (bits)');
	end

	toc;
end

