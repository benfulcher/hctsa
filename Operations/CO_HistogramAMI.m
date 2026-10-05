function out = CO_HistogramAMI(y, tau, meth, numBins)
% CO_HistogramAMI   Automutual information of a time series, from a 2-d histogram.
%
% Estimates the mutual information (in nats) between y(t) and y(t+tau) from the
% joint histogram of the two, with a Miller-Madow (Panzeri-Treves) correction for
% the upward bias of the plug-in estimate (about (Mxy - Mx - My + 1)/(2n) nats, for
% Mxy, Mx, My occupied joint and marginal bins and n samples). The corrected
% estimate is not clamped at zero, so it can be slightly negative for independent
% data.
%
% ---INPUTS:
% y, the input time series
% tau, the time lag, or a vector of lags (1 by default). Can also be 'ac' or 'tau',
%      to use the first zero-crossing of the autocorrelation function.
% meth, how the bin edges are chosen (default 'even'):
%       'even': evenly spaced bins through the range of the time series,
%       'std1', 'std2': numBins bins across -1 to 1 or -2 to 2 (the series is
%                       assumed to be z-scored), plus an extra bin beyond each
%                       edge that the series exceeds, to include outliers,
%       'quantiles': equiprobable bins chosen using quantiles (fewer bins if
%                    values are tied at a quantile).
% numBins, the number of bins (default 10)
%
% ---OUTPUTS:
% For a single lag, a scalar: the automutual information at that lag.
% For several lags, a structure with fields ami1, ami2, ami3, ami4, ami5, ...: the automutual information at
%       the first, second, ... lag in tau.

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
%% Check Inputs:
% ------------------------------------------------------------------------------
% Time-lag, tau
if nargin < 2 || isempty(tau)
	tau = 1;  % time-lag of 1
end
if ischar(tau) && ismember(tau, {'ac', 'tau'})
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
	if isnan(tau)
		out = NaN; return
	end
end

if nargin < 3 || isempty(meth)
	meth = 'even'; % default
end

if nargin < 4 || isempty(numBins)
	numBins = 10; % default number of bins: 10
end

% Number of options:
% remove outliers first?, number of bins, range of bins, bin sizes

% ------------------------------------------------------------------------------
%% Bins for the data:
% ------------------------------------------------------------------------------

% same for both -- assume same distribution (true for stationary processes,
% or small lags)
% (the edges are explicit; see BF_HistEdges and BF_QuantileEdges, which place
% edges so that a value exactly on one always falls in the same bin)
switch meth
	case 'even'
		b = BF_HistEdges(y, numBins); % through the range of the time series

	case 'std1' % bins out to +/- 1 std
		b = BF_HistEdges(y, numBins, [-1, 1]);
		if min(y) < b(1)
			b = [min(y) - 0.1, b];
		end
		if max(y) > b(end)
			b = [b, max(y) + 0.1];
		end

	case 'std2' % bins out to +/- 2 std
		b = BF_HistEdges(y, numBins, [-2, 2]);
		if min(y) < b(1)
			b = [min(y) - 0.1, b];
		end
		if max(y) > b(end)
			b = [b, max(y) + 0.1];
		end

	case 'quantiles' % use quantiles with ~equal number in each bin
		b = BF_QuantileEdges(y, numBins);

	otherwise
		error('Unknown method ''%s''', meth)
end

% Sometimes bins can be added (e.g., with std1 and std2), so need to redefine numBins
numBins = length(b) - 1; % number of bins (-1 since b defines edges)

% ------------------------------------------------------------------------------
% Form the time-delay vectors y1 and y2
% ------------------------------------------------------------------------------
amis = zeros(length(tau), 1);
for i = 1:length(tau)
	y1 = y(1:end - tau(i));
	y2 = y(1 + tau(i):end);

	% (1) Joint distribution of y1 and y2
	% (rows index y1 bins, columns index y2 bins)
	pij = histcounts2(y1, y2, b, b);
	pij = pij / sum(pij(:)); % joint
	pi = sum(pij, 2); % marginal over y1 (rows)
	pj = sum(pij, 1); % marginal over y2 (columns)

	pii = pi * ones(1, numBins);
	pjj = ones(numBins, 1) * pj;

	r = (pij > 0); % Defining the range in this way, we set log(0) = 0
	amis(i) = sum(pij(r) .* log(pij(r) ./ pii(r) ./ pjj(r)));

	% Miller-Madow (Panzeri-Treves) bias correction.
	%
	% The plug-in estimator above is biased *upwards* by approximately
	% (Mxy - Mx - My + 1)/(2n) nats, where the M are the numbers of *occupied*
	% joint and marginal bins and n is the sample size. For 10 bins that is
	% about 81/(2n): +0.41 nats at n=100 against +0.013 at n=3200. Since the
	% true AMI of an independent series is zero, the whole measured value was
	% this bias, decaying as 1/n -- i.e. the field tracked time-series length
	% rather than any property of the data. Subtracting the correction removes
	% the leading 1/n term.
	%
	% The result is deliberately not clamped at zero: for genuinely independent
	% data the corrected estimate should scatter symmetrically about 0, whereas
	% clamping would pile such series up on a single value.
	nSamples = numel(y1);
	Mxy = sum(r(:)); % occupied joint bins
	Mx = sum(pi > 0); % occupied y1 marginal bins
	My = sum(pj > 0); % occupied y2 marginal bins
	amis(i) = amis(i) - (Mxy - Mx - My + 1) / (2 * nSamples);
end

if length(tau) == 1
	out = amis;
else
	for i = 1:length(tau)
		out.(sprintf('ami%u', i)) = amis(i);
	end
end

end
