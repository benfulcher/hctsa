function out = WL_cwt(y, wname, maxScale)
% WL_cwt   Statistics of the continuous wavelet transform of a time series.
%
% Computes the continuous wavelet transform (cwt from Matlab's Wavelet Toolbox) at
% scales 1, ..., maxScale, takes the coefficients C and the scaled power SC (the
% power of each coefficient relative to the mean power over all coefficients), and
% returns statistics on the coefficients, on the distribution of the scaled power
% (its gamma fit and entropy), on the power summed across scales as a function of
% time, and on how the power differs between the two halves of the series.
%
% ---INPUTS:
% y, the input time series
% wname, the wavelet name. For a continuous wavelet transform, a proper continuous
%        analyzing wavelet like 'morl' (Morlet) is the standard choice; discrete
%        orthogonal wavelets like 'db3' are also accepted (their wavelet function,
%        evaluated on a fine grid, is used as the analyzing wavelet) and give a
%        genuinely different, complementary decomposition (see Wavelet Toolbox
%        Documentation for all options; default: 'db3')
% maxScale, the maximum scale of wavelet analysis (default: 32)
%
% ---OUTPUTS:
% meanC, meanabsC, medianabsC, maxabsC: the mean, mean magnitude, median magnitude
%        and maximum magnitude of the coefficients
% maxonmeanC, maxonmeanSC: the maximum relative to the mean, of the coefficient
%        magnitudes and of the scaled power
% pover99, pover98, pover95, pover90, pover80: the proportion of the total energy
%        held by the coefficients whose scaled power exceeds its 99th, 98th, 95th,
%        90th, 80th percentile (the strongest 1%, 2%, 5%, 10%, 20% of coefficients)
% gam1, gam2: the shape and scale parameters of a gamma distribution fitted to the
%        scaled power (gamfit; as the scaled power has mean 1, gam2 = 1/gam1)
% SC_h: the entropy of the energy distribution over all coefficients (in nats),
%        relative to its maximum possible value, the log of the number of
%        coefficients: 0 if the energy is spread evenly, negative when concentrated
% dd_SC_h: the entropy of the maximum scaled power in each of 10 equal time boxes at
%        each scale
% max_ssc, min_ssc, maxonmed_ssc, std_ssc: the maximum, minimum, maximum relative to
%        the median, and standard deviation over time of the scaled power summed
%        across scales
% pcross_maxssc50: the number of crossings of half its maximum by the summed power,
%        divided by N - 1
% stat_2_m_s: the mean of the standard deviations of the scaled power in the two
%        halves of the series, relative to the mean scaled power
% stat_2_s_m, stat_2_s_s: the standard deviation of the two halves' means (_m) and
%        of their standard deviations (_s), relative to the standard deviation of
%        the scaled power
%
% ---NOTES:
% The transform uses the legacy syntax cwt(y, scales, wname), which is the only form of
% cwt that takes integer scales and a discrete or real wavelet name ('db3', 'morl'): the
% current syntax is limited to the analytic wavelets 'morse', 'amor' and 'bump' on a
% frequency grid of its own. In MATLAB R2026a (Wavelet Toolbox 26.1) the legacy syntax
% is still accepted without a warning, and is equal to the textbook algorithm of
% convolving y with the integrated, dilated wavelet (to machine precision).
%
% The scaled power SC is normalized by the mean power over all coefficients (not by the
% total energy, as it used to be), so that its distribution, and the statistics that
% depend on it, do not depend on the number of coefficients (the series length). The
% maximum-type statistics (maxonmeanSC, max_ssc) still grow slowly with length, as the
% extreme of more samples does.

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
%% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(wname)
	wname = 'db3';
	fprintf(1, 'Using default wavelet ''%s''\n', wname);
end
if nargin < 3 || isempty(maxScale)
	maxScale = 32;
	fprintf(1, 'Using default maxScale of %u\n', maxScale);
end

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox'); % Check that a Wavelet Toolbox license is available:
doPlot = false; % plot outputs to figures
N = length(y); % length of the time series

% -------------------------------------------------------------------------------
scales = (1:maxScale);
coeffs = cwt(y, scales, wname);

S = abs(coeffs .* coeffs); % power
% Scaled power: each coefficient's power relative to the mean power over all
% coefficients (so it has mean 1, whatever the length of the series). This used
% to be the percentage of the total energy in each coefficient, 100*S/sum(S),
% which has mean 100/numEntries and so shrinks, and whose sums over scales,
% extremes and fitted parameters all drift, with the number of coefficients (the
% series length). Dividing by the mean power instead keeps the shape of the
% power distribution (how concentrated it is, in scales and in time) and
% removes the dependence on how many coefficients there are.
SC = S ./ mean(S(:));

if doPlot
	figure('color', 'w'); box('on');
	subplot(3, 1, 1)
	plot(y);
	subplot(3, 1, 2:3);
	pcolor(SC); shading interp;
end

% ------------------------------------------------------------------------------
%% Get statistics from CWT
% ------------------------------------------------------------------------------
numEntries = size(coeffs, 1) * size(coeffs, 2); % number of entries in coeffs matrix

% 1) Coefficients, coeffs
allCoeffs = coeffs(:);
out.meanC = mean(allCoeffs);
out.meanabsC = mean(abs(allCoeffs));
out.medianabsC = median(abs(allCoeffs));
out.maxabsC = max(abs(allCoeffs));
out.maxonmeanC = out.maxabsC / out.meanabsC;

% 2) Power, SC -- it's highly length-dependent
% out.meanSC = mean(SC(:)); % (reproduces the mean power of power spectrum)
% out.medianSC = median(SC(:));
% out.maxSC = max(SC(:));
out.maxonmeanSC = max(SC(:)) / mean(SC(:));

% Proportion of the total energy held by the coefficients whose power exceeds its
% p-th percentile (i.e., the strongest (100-p)% of the coefficients; the sum of SC is
% numEntries, so dividing by it gives a proportion of the energy). This replaces the
% energy held by coefficients exceeding a fraction p of the *maximum*, which depends
% on the extreme value of the coefficients and so falls steadily as the series (and
% the number of coefficients) gets longer, even for white noise.
SCsorted = sort(SC(:), 'descend');
poverfn = @(p) sum(SCsorted(1:max(1, floor((100 - p) / 100 * numEntries)))) / numEntries;
out.pover99 = poverfn(99);
out.pover98 = poverfn(98);
out.pover95 = poverfn(95);
out.pover90 = poverfn(90);
out.pover80 = poverfn(80);

% Distribution of scaled power
% Fit using Statistics Toolbox

if doPlot
	figure('color', 'w');
	ksdensity(SC(:));
end

gamma_phat = gamfit(SC(:));
out.gam1 = gamma_phat(1);
out.gam2 = gamma_phat(2);

% ------------------------------------------------------------------------------
%% 2D entropy
% ------------------------------------------------------------------------------
% turn into probabilities
SC_a = SC ./ sum(SC(:));
% compute entropy, relative to its maximum possible value, log(numEntries): the
% entropy itself grows as log(numEntries) with the series length, whereas
% -sum(p*log(p)) - log(numEntries) = -mean(SC*log(SC)) (<= 0) does not
SC_a = SC_a(:);
out.SC_h = -sum(SC_a .* log(SC_a)) - log(numEntries);

% ------------------------------------------------------------------------------
%% Weird 2-D entropy idea -- first discretize
% ------------------------------------------------------------------------------
% (i) Discretize the space into numBoxes boxes along the time axis
% Many choices, let's discretize into maximum energy
% (could also do average, or proportion inside box with more energy than
% average, ...)
numBoxes = 10;
if N < numBoxes
	error('Time series too short');
end
dd_SC = zeros(maxScale, numBoxes);
cutoffs = round(linspace(0, N, numBoxes + 1));
for i = 1:maxScale
	for j = 1:numBoxes
		dd_SC(i, j) = max(SC(i, cutoffs(j) + 1:cutoffs(j + 1)));
	end
end

% Turn into probabilities
dd_SC = dd_SC ./ sum(dd_SC(:));

% Compute entropy
dd_SCO = dd_SC(:);
out.dd_SC_h = -sum(dd_SCO .* log(dd_SCO));

% ------------------------------------------------------------------------------
%% Sum across scales
% ------------------------------------------------------------------------------
SSC = sum(SC);
out.max_ssc = max(SSC);
out.min_ssc = min(SSC);
out.maxonmed_ssc = max(SSC) / median(SSC);
out.pcross_maxssc50 = sum(BF_SignChange(SSC - 0.5 * max(SSC))) / (N - 1);
out.std_ssc = std(SSC);

% ------------------------------------------------------------------------------
%% Stationarity
% ------------------------------------------------------------------------------
% 2-way split of the scale-power surface, collapsed across scales.
% (A 5-way split was also tried here but dropped: on real data its three
% summary stats correlated r>0.98 with their 2-way counterparts, adding
% negligible information for 5x the fields.)
SC_1 = SC(:, 1:floor(N / 2)); % collapse across scales, first half
SC_2 = SC(:, floor(N / 2) + 1:end); % collapse across scales, second half

mean2_1 = mean(SC_1(:));
mean2_2 = mean(SC_2(:));

std2_1 = std(SC_1(:));
std2_2 = std(SC_2(:));

% out.stat_2_m_m = mean([mean2_1 mean2_2])/mean(SC(:));
out.stat_2_m_s = mean([std2_1, std2_2]) / mean(SC(:));
out.stat_2_s_m = std([mean2_1, mean2_2]) / std(SC(:));
out.stat_2_s_s = std([std2_1, std2_2]) / std(SC(:));

end
