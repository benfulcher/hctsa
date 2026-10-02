function out = CO_AddNoise(y, tau, amiMethod, extraParam, randomSeed)
% CO_AddNoise   How the automutual information of the series falls as noise is added.
%
% Adds independent Gaussian noise of standard deviation eta to the (z-scored)
% series, for 50 noise levels evenly spaced from eta = 0 to eta = 3, drawing
% fresh noise at each level, and measures the automutual information (AMI) at
% lag tau at each level. The AMI falls as the noise swamps the signal; the outputs
% describe the resulting curve of AMI against eta, including a fit to an
% exponential decay.
%
% The AMI can be estimated using histograms with extraParam bins (implemented in
% CO_HistogramAMI) or using the Information Dynamics Toolkit (IN_AutoMutualInfo).
% The AMI is in nats for all methods (JIDT's Kraskov estimator, like its Gaussian
% estimator, uses natural logarithms).
%
% This algorithm is quite different from, but was based on the idea of, noise
% titration, presented in Poon and Barahona (2001).
%
% ---INPUTS:
% y, the input time series (should be z-scored)
% tau, the time delay for computing the AMI (a number of samples, or 'ac' for the
%       first zero-crossing of the autocorrelation function of y)
% amiMethod, the method for computing the AMI:
%       * 'std1', 'std2', 'quantiles', 'even': histogram-based estimation
%         (see CO_HistogramAMI)
%       * 'gaussian', 'kernel', 'kraskov1', 'kraskov2': estimation using JIDT
%         (see IN_AutoMutualInfo)
%       Default: 'even'.
% extraParam, a parameter of the estimator: the number of bins for the histogram
%       methods (CO_HistogramAMI), or the number of nearest neighbors for the
%       Kraskov methods (IN_AutoMutualInfo)
% randomSeed, how to reset the random seed, using BF_ResetSeed, for reproducible
%       results
%
% ---OUTPUTS: statistics of the AMI as a function of noise level (50 levels):
% pdec, the proportion of steps on which the AMI decreases
% meanch, the mean change in AMI per step
% ac1, ac2, the autocorrelation of the sequence of AMI values, at lags 1 and 2
% firstUnder75, firstUnder50, firstUnder25: the noise level at which the AMI
%       first falls below 75%, 50% and 25% of its value without noise (the
%       largest noise level, 3, if it never does)
% ami_at_5, ami_at_10, ami_at_15, ami_at_20: the AMI at the first noise level at or
%       above eta = 0.5, 1, 1.5 and 2
% pcrossmean, the proportion of steps on which the AMI curve crosses its mean
% fitexpa, fitexpb, fitexpr2, fitexpadjr2, fitexprmse: the amplitude a, rate b,
%       R^2, adjusted R^2 and root-mean-square error of a fit of a * exp(b * eta)
%       (requires the Curve Fitting Toolbox)
% fitlina, fitlinb, linfit_mse: the slope, intercept and mean squared error of a
%       straight-line fit
%
% ---REFERENCES:
% C.-S. Poon and M. Barahona, "Titration of chaos with added noise", Proc. Natl.
% Acad. Sci. USA 98(13), 7107 (2001).

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

% -------------------------------------------------------------------------------
% Preliminary checks

% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = false; % plot outputs to figure

% -------------------------------------------------------------------------------
%% Check inputs:

% Expecting a z-scored input time series:
BF_iszscored(y);

if nargin < 2
	tau = []; % set default in CO_HistogramAMI
end
% Set tau to minimum of autocorrelation function
if ~isempty(tau) && ischar(tau) && ismember(tau, {'ac', 'tau'})
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
end
if nargin < 3
	amiMethod = 'even'; % using evenly spaced bins in CO_HistogramAMI
end
if nargin < 4
	extraParam = []; % number of bins for CO_HistogramAMI
end
if nargin < 5
	randomSeed = [];
end

% -------------------------------------------------------------------------------
%% Preliminaries
% -------------------------------------------------------------------------------

% Set up noise range:
BF_ResetSeed(randomSeed); % reset the random seed if specified
noiseRange = linspace(0, 3, 50); % compare properties across this noise range
numRepeats = length(noiseRange);

% ------------------------------------------------------------------------------
%% Compute the automutual information across a range of noise levels
% ------------------------------------------------------------------------------
% An independent, uncorrelated Gaussian noise vector is drawn at each noise
% level, so that each point on the curve is an independent sample of the
% AMI-vs-noise relationship rather than all points sharing one noise draw's
% idiosyncrasies (rescaled by increasing standard deviation).
amis = zeros(numRepeats, 1); % preassign
switch amiMethod
	case {'std1', 'std2', 'quantiles', 'even'}
		% histogram-based methods using my naive implementation in CO_Histogram
		for i = 1:numRepeats
			noise = randn(size(y)); % fresh uncorrelated additive noise at this level
			amis(i) = CO_HistogramAMI(y + noiseRange(i) * noise, tau, amiMethod, extraParam);
			if isnan(amis(i))
				error('Error computing AMI: Time series too short (?)');
			end
		end
	case {'gaussian', 'kernel', 'kraskov1', 'kraskov2'}
		for i = 1:numRepeats
			noise = randn(size(y)); % fresh uncorrelated additive noise at this level
			amis(i) = IN_AutoMutualInfo(y + noiseRange(i) * noise, tau, amiMethod, extraParam);
			if isnan(amis(i))
				error('Error computing AMI: Time series too short (?)');
			end
		end
end

% -------------------------------------------------------------------------------
%% Output statistics
% -------------------------------------------------------------------------------

% Proportion decreases:
out.pdec = sum(diff(amis) < 0) / (numRepeats - 1);

% Mean change in AMI:
out.meanch = mean(diff(amis));

% Autocorrelation of AMIs:
out.ac1 = CO_AutoCorr(amis, 1, 'Fourier');
out.ac2 = CO_AutoCorr(amis, 2, 'Fourier');

% Noise level required to reduce ami to proportion x of its initial value:
firstUnderVals = [0.75, 0.5, 0.25];
for i = 1:length(firstUnderVals)
	out.(sprintf('firstUnder%u', firstUnderVals(i) * 100)) = ...
					firstUnder_fn(firstUnderVals(i) * amis(1), noiseRange, amis);
end

% AMI at actual noise levels: 0.5, 1, 1.5 and 2
noiseLevels = [0.5, 1, 1.5, 2];
for i = 1:length(noiseLevels)
	out.(sprintf('ami_at_%u', noiseLevels(i) * 10)) = ...
			amis(find(noiseRange >= noiseLevels(i), 1, 'first'));
end

% Count number of times the AMI function crosses its mean
out.pcrossmean = sum(BF_SignChange(amis - mean(amis))) / (numRepeats - 1);

% -------------------------------------------------------------------------------
% Fit exponential decay (using Curve Fitting Toolbox)
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [amis(1) -1]);
f = fittype('a*exp(b*x)', 'options', s);
[c, gof] = fit(noiseRange', amis, f);

% Output statistics on fit to an exponential decay
out.fitexpa = c.a;
out.fitexpb = c.b;
out.fitexpr2 = gof.rsquare;
out.fitexpadjr2 = gof.adjrsquare;
out.fitexprmse = gof.rmse;

% ------------------------------------------------------------------------------
% Fit linear function:
p = polyfit(noiseRange', amis, 1);
out.fitlina = p(1); % gradient
out.fitlinb = p(2); % intercept
linfit = polyval(p, noiseRange);
out.linfit_mse = mean((linfit' - amis).^2);

% -------------------------------------------------------------------------------
% Plot output:
if doPlot
	figure('color', 'w'); box('on');
	cc = BF_GetColorMap('set1', 2, 1);
	% figure('color','w');
	hold on; box('on')
	plot(noiseRange, c.a * exp(c.b * noiseRange), 'color', cc{2}, 'linewidth', 2)
	plot(noiseRange, amis, '.-', 'color', cc{1})
	xlabel('\eta'); ylabel('AMI_1')
end

% -------------------------------------------------------------------------------
function firsti = firstUnder_fn(x, m, p)
	% Find the value of m for the first time p goes under the threshold, x
	% p and m vectors of the same length

	firsti = m(find(p < x, 1, 'first'));

	% If it never goes under -- saturate as m at the maximum
	% (could be NaN, but this is more interpretable/comparable)
	if isempty(firsti)
		firsti = m(end);
	end

end

end
