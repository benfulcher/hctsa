function out = DN_SimpleFit(x, dmodel, numBins)
% DN_SimpleFit   Fits a simple curve to the distribution of the values.
%
% Fits a simple parametric curve to an estimate of the distribution of values in
% the time series, ignoring their temporal ordering. The distribution is estimated
% either as a histogram, with a specified number of bins, or as a
% kernel-smoothed density (ksdensity, with its default width). The outputs
% measure the goodness of fit, and test the residuals (in order of increasing
% value) for remaining structure. The curve is fitted by least squares to the
% density, in a deterministic way that does not depend on random starts or on an
% optimizer's default settings (BF_FitDensityCurve).
%
% ---INPUTS:
% x, the input time series
% dmodel, the distribution model to fit:
%           (i) 'gauss1': a single Gaussian
%           (ii) 'gauss2': a sum of two Gaussians
%           (iii) 'exp1': an exponential, a*exp(b*x)
%           (iv) 'power1': a power law, a*x^b (cannot be fit if any bin center
%                   is not positive; NaN is returned)
% NaN is also returned if there are no more bins than parameters of the model.
% numBins, how to estimate the distribution (default: 'sqrt'):
%           a text option: the name of a binning rule for histcounts,
%                   e.g., 'sqrt' uses the square root of the number of data
%                   points as the number of bins
%           a positive integer: the number of bins in the histogram
%           0: use ksdensity instead of a histogram
%
% ---OUTPUTS:
% r2, the R^2 goodness of fit
% adjr2, R^2 adjusted for the number of fitted parameters
% rmse, the root-mean-square error of the fit, in units of probability density of
%       the standardized series (the fit is to the density: the histogram counts
%       divided by the number of points and the bin width, or the ksdensity
%       estimate; the error is multiplied by the standard deviation of x), so it
%       does not depend on the length or the scale of the series
% resAC1, resAC2, the autocorrelation of the residuals, in order of
%       increasing value, at lags 1 and 2
% resrunsz, the signed z-statistic of a runs test on the residuals, in order of
%       increasing value (BF_RunsZ): negative when the residuals have fewer runs
%       about their median than expected for a random order
%
% ---NOTES:
% Fits of time-series models (sinusoids or Fourier series) against time have
% moved to SP_SinusoidFit, since they depend on the temporal ordering of the
% data. For backward compatibility, the time-series models ('sin1', 'sin2',
% 'sin3', 'fourier1', 'fourier2', 'fourier3') are still accepted here: the
% call is passed on to SP_SinusoidFit (with identical outputs) after issuing a
% one-time 'hctsa:deprecated' warning.

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

	% Compute the distribution (histogram, normalized to a probability density):
	if ischar(numBins) % specify a binning method
		[dny, binEdges] = histcounts(x, 'BinMethod', numBins);
		dnx = mean([binEdges(1:end - 1); binEdges(2:end)]);
		dny = dny / (sum(dny) * mean(diff(binEdges))); % counts -> probability density
	elseif numBins == 0 % use ksdensity instead of a histogram
		[dny, dnx] = ksdensity(x);
	else
		[dny, binEdges] = histcounts(x, numBins);
		dnx = mean([binEdges(1:end - 1); binEdges(2:end)]);
		dny = dny / (sum(dny) * mean(diff(binEdges))); % counts -> probability density
	end

	% Both must be column vectors:
	if size(dnx, 2) > size(dnx, 1)
		dnx = dnx';
		dny = dny';
	end

	% Fit the distribution model, by least squares for the density (BF_FitDensityCurve):
	if strcmp(dmodel, 'power1') && any(dnx <= 0)
		fprintf(1, 'The model ''%s'' can not be applied to non-positive data\n', dmodel);
		out = NaN; return
	end
	fitModels = {'gauss1', 'gauss', 3; 'gauss2', 'gauss2', 6; 'exp1', 'exp', 2; 'power1', 'power', 2}; % model, curve, number of parameters
	iModel = find(strcmp(fitModels(:, 1), dmodel));
	dnyFit = BF_FitDensityCurve(dnx, dny, fitModels{iModel, 2});

	% Residuals (in order of increasing value) and goodness of fit as in the Curve Fitting
	% Toolbox: R^2, R^2 adjusted for the number of fitted parameters, and the root-mean-square
	% error from the residual sum of squares divided by the degrees of freedom of the error
	res = dny - dnyFit;
	sse = sum(res.^2);
	sstot = sum((dny - mean(dny)).^2);
	dfe = length(dny) - fitModels{iModel, 3}; % degrees of freedom of the error
	if dfe < 1 % no more bins than parameters: the fit is not meaningful
		out = NaN; return
	end
	r2 = 1 - sse/sstot;
	adjr2 = 1 - (1 - r2)*(length(dny) - 1)/dfe;
	rmse = sqrt(sse/dfe);

else
	error('Invalid distribution model ''%s'' specified', dmodel);
end

% ------------------------------------------------------------------------------
%% Compute the outputs into a structure
% ------------------------------------------------------------------------------
out.r2 = r2; % rsquared
out.adjr2 = adjr2; % degrees of freedom-adjusted rsquared (not currently registered
                   % by any mop -- redundant with r2 for these fixed-order fits)

% Root mean square error. The fit was done directly to the probability density, so
% multiplying by std(x) expresses it in density units of the standardized series,
% which does not grow with the length of the series (histogram counts do) or
% depend on the scale of x:
out.rmse = rmse * std(x);

% Remaining structure in the residuals, in order of increasing value:
[out.resAC1, out.resAC2, out.resrunsz] = BF_ResidualStats(res, sstot);

end
