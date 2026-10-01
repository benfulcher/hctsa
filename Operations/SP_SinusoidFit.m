function out = SP_SinusoidFit(y, model)
% SP_SinusoidFit   Fit sinusoids or a Fourier series to the time series.
%
% Fits a sum of 1-3 sinusoids, or a Fourier series with 1-3 terms, to the time
% series y as a function of its time index t = 1:N (the values are fitted in the
% order in which they occur, so the result depends on the temporal ordering of
% the data). Uses the 'fit' function from Matlab's Curve Fitting Toolbox.
%
% The fitted models are:
%   'sinK':     a sum of K sinusoids, sum_{i=1..K} a_i*sin(b_i*t + c_i), with
%               free amplitudes a_i, frequencies b_i, and phases c_i,
%   'fourierK': a K-term Fourier series,
%               a0 + sum_{i=1..K} (a_i*cos(i*w*t) + b_i*sin(i*w*t)),
%               with a single fitted fundamental frequency w (so that the terms
%               are harmonically related).
%
% Goodness of fit is summarized using the root mean square error (and R^2), and
% the residuals of the fit are then characterized by their autocorrelation at
% lags 1 and 2 and by a runs test, which together reveal remaining temporal
% structure that the model has not captured (e.g., whether a periodic component
% has been fully explained).
%
% ---INPUTS:
% y, the input time series (a vector; row vectors are converted to columns)
% model, the model to fit:
%       (i) 'sin1': a single sinusoid
%       (ii) 'sin2': a sum of two sinusoids
%       (iii) 'sin3': a sum of three sinusoids
%       (iv) 'fourier1': a Fourier series with one term
%       (v) 'fourier2': a Fourier series with two terms
%       (vi) 'fourier3': a Fourier series with three terms
%
% ---OUTPUTS: a structure containing
% r2, the R^2 of the fit
% adjr2, the degrees-of-freedom-adjusted R^2
% rmse, the root mean square error of the fit
% resAC1, the autocorrelation of the residuals at lag 1 (using the 'Fourier'
%         method of CO_AutoCorr)
% resAC2, the autocorrelation of the residuals at lag 2
% resruns, the p-value of a runs test on the residuals (HT_IndependenceTests,
%         'runstest')
% If the model cannot be fitted (NaN or Inf computed by the model function), NaN
% is returned instead of a structure.
%
% ---NOTES:
% This function holds the time-series-model branch of the former DN_SimpleFit,
% from which it was split because the distribution of values is unaffected by
% temporal ordering, whereas these fits are not. r2 and adjr2 are not registered
% for the sin1/sin2/sin3 mops (rmse is).

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
%% Fit the model
% ------------------------------------------------------------------------------
TSmodels = {'sin1', 'sin2', 'sin3', 'fourier1', 'fourier2', 'fourier3'}; % valid time-series models

if ~ismember(model, TSmodels)
	error('Invalid time-series model ''%s'' specified', model);
end

if size(y, 2) > size(y, 1)
	y = y';
end % y must be a column vector
t = (1:length(y))'; % Time variable for equal sampling of the univariate time series
try
	[cfun, gof, output] = fit(t, y, model); % fit the model
catch emsg % this model can't even be fitted OR license problem
	if strcmp(emsg.message, 'NaN computed by model function.') || strcmp(emsg.message, 'Inf computed by model function.')
		fprintf(1, 'The model %s failed for this data -- returning NaNs for all fitting outputs\n', model);
		out = NaN; return
	else
		error('Unexpected error fitting ''%s'' to the time series', model)
	end
end

% ------------------------------------------------------------------------------
%% Compute the outputs into a structure
% ------------------------------------------------------------------------------
out.r2 = gof.rsquare; % rsquared (not currently registered by the sin1/sin2/sin3 mops,
                       % which register rmse instead)
out.adjr2 = gof.adjrsquare; % degrees of freedom-adjusted rsquared (not currently registered
                             % by any mop -- redundant with r2 for these fixed-order fits)

out.rmse = gof.rmse; % root mean square error
out.resAC1 = CO_AutoCorr(output.residuals, 1, 'Fourier'); % autocorrelation of residuals at lag 1
out.resAC2 = CO_AutoCorr(output.residuals, 2, 'Fourier'); % autocorrelation of residuals at lag 2
out.resruns = HT_IndependenceTests(output.residuals, 'runstest'); % runs test on residuals -- outputs p-value

end
