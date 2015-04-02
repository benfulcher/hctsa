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

% function [mis] checkMiNullDistribution(repeats, observations, bias1, bias2)
%
% Return and plot a sample of MI values for a given number of observations and the given biases of the individual
%  variables, with no relation between the first and second variable.
%
% Inputs:
% - repeats - number of surrogate MI calculations to perform. Best to have at least 1000 here.
% - observations - length of time series for each calculation
% - bias1 - proportion of values for variable 1 to be assigned 0 instead of 1. 0.5 means 50-50 probability
% - bias2 - proportion of values for variable 2 to be assigned 0 instead of 1. 0.5 means 50-50 probability
% - useToolkit - whether or not to use the toolkit to generate the bootstrapped distribution. Default is false, and it is better not to do this, since it fixes the bias precisely for each time series sample, restricting the number of possible MI values that could be measured. While the results from setting this still follow the analytic distribution, they do so in a stepped manner (since there is a restricted set of allowable MI values)

function [mis] = checkMiDiscreteNullDistribution(repeats, observations, bias1, bias2, useToolkit)

	addpath('..');
	javaaddpath('../../../infodynamics.jar');

	if (nargin < 5)
		useToolkit = false;
	end

	if (useToolkit)
		% 1. We can let the toolkit compute the distribution of MIs for us,
		%  by first supplying the variables with the correct bias.
		% It is NOT recommended to use this approach however, as described in the header comments.
		miCalc=javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2);
		x = [zeros(1,observations*bias1) ones(1,observations*(1-bias1))]; % Create data with exact bias
		y = [zeros(1,observations*bias2) ones(1,observations*(1-bias2))]; % Create data with exact bias
		fprintf('Creating surrogates from x and y of lengths %d and %d with biases %.3f and %.3f \n', length(x), length(y), sum(x == 0)/length(x), sum(y == 0)/length(y));
		miCalc.initialise();
		miCalc.addObservations(octaveToJavaIntArray(x), octaveToJavaIntArray(y));
		nullDist = miCalc.computeSignificance(repeats); % Compute null distribution
		mis = javaMatrixToOctave(nullDist.distribution);	
	else
		% OR
		% 2. We could compute the bootstrapped distribution of MIs by bootstrapping ourselves:
		mis = zeros(1, repeats);
		miCalc=javaObject('infodynamics.measures.discrete.MutualInformationCalculatorDiscrete', 2);
		for s = 1 : repeats
			x = (rand(1,observations) < bias1)*1; % Create data with sampled bias
			y = (rand(1,observations) < bias2)*1; % Create data with sampled bias
			miCalc.initialise();
			miCalc.addObservations(octaveToJavaIntArray(x), octaveToJavaIntArray(y));
			mis(s) = miCalc.computeAverageLocalOfObservations();
		end
	end

	bins = 100;

	% Got all the samples now
	[pdfY,pdfX] = hist(mis, bins);
	pdfY = pdfY / repeats; % Normalise the sum of the pdf to 1.
	% Compute CDF
	cdfY = cumsum(pdfY);
	figure(1);
	plot(pdfX, cdfY, 'rx', 'MarkerSize', 10);
	hold on;
	% Compare to the theoretical \chi^2 distribution
	plot(pdfX, chi2cdf(2*observations*pdfX*log(2),1), 'go', 'MarkerSize', 10);
	% title('CDFs of MI');
	hold off;
	legend(['Bootstrapped'; 'Chi^2 analytic'], 'Location', 'East');
	axis([0 max(pdfX) 0 1]);
	% Draw a line at the alpha=0.05 cutoff
	line([0 max(pdfX)], [0.95 0.95], 'linewidth', 2, 'color', 'blue');
	set (gca,'fontsize',26);
	xlabel('MI', 'FontSize', 36, 'FontWeight', 'bold');
	ylabel('CDF(MI)', 'FontSize', 36, 'FontWeight', 'bold');
	print(sprintf('cdfMiDiscrete-N%d-p-%.2f-%.2f.eps', observations, bias1, bias2), '-depsc');

end

