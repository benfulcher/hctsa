function out = DN_SimpleFit(x, dmodel, numBins)
% DN_SimpleFit  Fit a simple model to the distribution of time-series values.
%
% Uses the 'fit' function from Matlab's Curve Fitting Toolbox to fit a simple
% parametric model to an estimate of the distribution of values in the time
% series, ignoring their temporal ordering.
%
% The distribution of time-series values is estimated using either a
% kernel-smoothed density via the Matlab function ksdensity with the default
% width parameter, or by a histogram with a specified number of bins, numBins.
%
% Deprecated: fits of time-series models (sinusoids or Fourier series) against
% time have moved to SP_SinusoidFit, since they depend on the temporal ordering
% of the data. For backward compatibility, the time-series models ('sin1',
% 'sin2', 'sin3', 'fourier1', 'fourier2', 'fourier3') are still accepted here:
% the call is passed on to SP_SinusoidFit (with identical outputs) after issuing
% a one-time 'hctsa:deprecated' warning.
%
% ---INPUTS:
% x, the input time series
%
% dmodel, the distribution model to fit:
%           (i) 'gauss1'
%           (ii) 'gauss2'
%           (iii) 'exp1'
%           (iv) 'power1'
%       (The deprecated time-series models 'sin1', 'sin2', 'sin3', 'fourier1',
%       'fourier2', 'fourier3' are dispatched to SP_SinusoidFit.)
%
% numBins, the number of bins for a histogram-estimate of the distribution of
%       time-series values. If numBins = 0, uses ksdensity instead of histogram.
%
% ---OUTPUTS: the goodness of fit, R^2, root mean square error, the
% autocorrelation of the residuals, and a runs test on the residuals.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------

% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

% ------------------------------------------------------------------------------
%% Deprecated: time-series models now live in SP_SinusoidFit
% ------------------------------------------------------------------------------
persistent warnedTSmodel % warn only once per session
TSmodels = {'sin1', 'sin2', 'sin3', 'fourier1', 'fourier2', 'fourier3'}; % time-series models
if ischar(dmodel) && any(strcmp(TSmodels, dmodel))
	if isempty(warnedTSmodel)
		warning('hctsa:deprecated', ['DN_SimpleFit no longer fits time-series models; ' ...
				'''%s'' is dispatched to SP_SinusoidFit (this warning is shown once per session).'], dmodel);
		warnedTSmodel = true;
	end
	out = SP_SinusoidFit(x, dmodel);
	return
end

% ------------------------------------------------------------------------------
%% Fit the model
% ------------------------------------------------------------------------------
distModels = {'gauss1', 'gauss2', 'exp1', 'power1'}; % valid distribution models

if any(strcmp(distModels, dmodel)) % valid DISTRIBUTION model name
	if nargin < 3 || isempty(numBins) % haven't specified numBins
		numBins = 'sqrt'; % use sqrt of number of data points
	end

	% Compute the histogram counts:
	if ischar(numBins) % specify a binning method
		[dny, binEdges] = histcounts(x, 'BinMethod', numBins);
		dnx = mean([binEdges(1:end - 1); binEdges(2:end)]);
	elseif numBins == 0 % use ksdensity instead of a histogram
		[dny, dnx] = ksdensity(x);
	else
		[dny, binEdges] = histcounts(x, numBins);
		dnx = mean([binEdges(1:end - 1); binEdges(2:end)]);
	end

	% Both must be column vectors:
	if size(dnx, 2) > size(dnx, 1)
		dnx = dnx';
		dny = dny';
	end

	% Fit the distribution model:
	try
		[cfun, gof, output] = fit(dnx, dny, dmodel);
	catch emsg % this model can't even be fitted OR license problem...
		if strcmp(emsg.identifier, 'curvefit:fit:nanComputed') ...
				|| strcmp(emsg.identifier, 'curvefit:fit:infComputed')
			fprintf(1, 'Error fitting the model ''%s'' to this data:\n%s', dmodel, emsg.message);
			out = NaN; return
		elseif strcmp(emsg.message, 'Power functions cannot be fit to non-positive xdata.') ...
				|| strcmp(emsg.identifier, 'curvefit:fit:powerFcnsRequirePositiveData')
			fprintf(1, 'The model ''%s'' can not be applied to non-positive data\n', dmodel);
			out = NaN; return
		else
			error('Error fitting %s to the data distribution\n%s', dmodel, emsg.message);
		end
	end

else
	error('Invalid distribution model ''%s'' specified', dmodel);
end

% ------------------------------------------------------------------------------
%% Compute the outputs into a structure
% ------------------------------------------------------------------------------
out.r2 = gof.rsquare; % rsquared
out.adjr2 = gof.adjrsquare; % degrees of freedom-adjusted rsquared (not currently registered
                             % by any mop -- redundant with r2 for these fixed-order fits)

out.rmse = gof.rmse; % root mean square error
out.resAC1 = CO_AutoCorr(output.residuals, 1, 'Fourier'); % autocorrelation of residuals at lag 1
out.resAC2 = CO_AutoCorr(output.residuals, 2, 'Fourier'); % autocorrelation of residuals at lag 2
out.resruns = HT_HypothesisTest(output.residuals, 'runstest'); % runs test on residuals -- outputs p-value

end
