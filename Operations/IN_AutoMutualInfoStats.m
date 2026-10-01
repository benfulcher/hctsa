function out = IN_AutoMutualInfoStats(y, maxTau, estMethod, extraParam)
% IN_AutoMutualInfoStats   Statistics on the automutual information function of a time series.
%
% Computes the automutual information (AMI) at each time lag from 1 to maxTau
% (using IN_AutoMutualInfo), and returns the AMI values and statistics on their
% pattern across lags: mean and spread, the first local minimum, the number and
% spacing of local maxima and minima, the proportion of crossings of the mean,
% median and 10th and 90th percentiles, and the lag-1 autocorrelation of the AMI
% function. With a Gaussian estimator the AMI is -0.5*log(1 - r^2) for the Pearson
% correlation r between the series and its lagged copy; the other estimators are
% implemented in the Java Information Dynamics Toolkit (JIDT).
%
% ---INPUTS:
% y, column vector of time series data
% maxTau, the maximal time delay (default: ceil(N/4), for series length N; it is
%    reduced to ceil(N/2) if larger than that)
% estMethod, the estimation method for the AMI (default: 'kernel'), one of
%    'gaussian', 'kernel', 'kraskov1', 'kraskov2'; cf. IN_AutoMutualInfo
% extraParam, an extra parameter of the estimator (default: none); for
%    'kraskov1' and 'kraskov2', the number of nearest neighbors, as a string
%    (default: '4'); cf. IN_AutoMutualInfo
%
% ---OUTPUTS:
% A structure with the following fields (the AMI at each lag is returned as ami1,
% ami2, ..., up to maxTau; NaN where the series is too short for that lag):
% ami1, ami2, ami3, ami4, ami5, ami6, ami7, ami8, ami9, ami10, ami11, ami12,
% ami13, ami14, ami15, ami16, ami17, ami18, ami19, ami20, ami21, ami22,
% ami23, ami24, ami25, ami26, ami27, ami28, ami29, ami30, ami31, ami32,
% ami33, ami34, ami35, ami36, ami37, ami38, ami39, ami40
% mami, the mean of the AMI across lags
% stdami, the standard deviation of the AMI across lags
% pextrema, the number of local extrema (peaks and troughs) of the AMI function,
%    as a proportion of the number of lags
% fmmi, the lag of the first local minimum of the AMI function (the number of lags
%    if there is none)
% sumami_fmmi, the sum of the AMI from lag 1 to fmmi
% pmaxima, the number of intervals between successive local maxima, divided by
%    floor(number of lags/2)
% modeperiodmax, the most common spacing between successive local maxima (NaN if
%    fewer than two maxima)
% pmodeperiodmax, the proportion of spacings between successive local maxima that
%    equal modeperiodmax
% pminima, the number of intervals between successive local minima, divided by
%    floor(number of lags/2)
% modeperiodmin, the most common spacing between successive local minima (NaN if
%    fewer than two minima)
% pmodeperiodmin, the proportion of spacings between successive local minima that
%    equal modeperiodmin
% pcrossmean, the proportion of successive lags at which the AMI function crosses
%    its mean
% pcrossmedian, ... crosses its median
% pcrossq10, ... crosses its 10th percentile
% pcrossq90, ... crosses its 90th percentile
% amiac1, the lag-1 autocorrelation of the AMI function (CO_AutoCorr)

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
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------

% maxTau: the maximum time delay to investigate
if nargin < 2 || isempty(maxTau)
	maxTau = ceil(N / 4);
end
maxTau0 = maxTau;

% Don't go above N/2
maxTau = min(maxTau, ceil(N / 2));

% Estimation method:
if nargin < 3
	estMethod = '';
end

% extraParam
if nargin < 4
	extraParam = [];
end

% ------------------------------------------------------------------------------
%% Get the AMI data:
% ------------------------------------------------------------------------------
ami = IN_AutoMutualInfo(y, 1:maxTau, estMethod, extraParam);

% Convert structure to a vector
ami = struct2cell(ami);
ami = [ami{:}];

% ------------------------------------------------------------------------------
% Output the raw values:
% ------------------------------------------------------------------------------
for i = 1:maxTau0
	if i <= maxTau
		out.(sprintf('ami%u', i)) = ami(i);
	else % we've trimmed the maximum back because time series is too short
		out.(sprintf('ami%u', i)) = NaN;
	end
end

% ------------------------------------------------------------------------------
% Output statistics:
% ------------------------------------------------------------------------------
lami = length(ami);

% Mean and std of automutual information over this range of time delays:
out.mami = mean(ami);
out.stdami = std(ami);

% First minimum of mutual information across range
dami = diff(ami);
extremai = find(dami(1:end - 1) .* dami(2:end) < 0);
out.pextrema = length(extremai) / (lami - 1);

% First local MINIMUM specifically (extremai above contains both minima and
% maxima; a minimum is where ami is still decreasing into the sign change,
% i.e., dami(i) < 0, and it occurs at ami-index i+1 since dami(i) = ami(i+1)-ami(i)):
minimai = extremai(dami(extremai) < 0) + 1;
if isempty(minimai)
	out.fmmi = lami; % no local minimum found in this range
else
	out.fmmi = min(minimai);
end

% Integrated AMI up to the first minimum: how much automutual information
% accumulates before the series first decorrelates, as opposed to mami's
% flat average over the whole (often plateaued) lag range:
out.sumami_fmmi = sum(ami(1:out.fmmi));

% ----Look for periodicities in local maxima
maximai = find(dami(1:end - 1) > 0 & dami(2:end) < 0) + 1;
dmaximai = diff(maximai);
% Is there a big peak in dmaxima?
% (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pmaxima = length(dmaximai) / floor(lami / 2);
if isempty(dmaximai) % fewer than 2 local maxima
	out.modeperiodmax = NaN;
	out.pmodeperiodmax = NaN;
else
	out.modeperiodmax = mode(dmaximai);
	out.pmodeperiodmax = sum(dmaximai == mode(dmaximai)) / length(dmaximai);
end

% ----Look for periodicities in local minima
minimai = find(dami(1:end - 1) < 0 & dami(2:end) > 0) + 1;
dminimai = diff(minimai);
% Is there a big peak in dminima?
% (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pminima = length(dminimai) / floor(lami / 2);
if isempty(dminimai) % fewer than 2 local maxima
	out.modeperiodmin = NaN;
	out.pmodeperiodmin = NaN;
else
	out.modeperiodmin = mode(dminimai);
	out.pmodeperiodmin = sum(dminimai == mode(dminimai)) / length(dminimai);
end

% ----Number of crossings at mean/median level, percentiles
out.pcrossmean = mean(BF_SignChange(ami - mean(ami)));
out.pcrossmedian = mean(BF_SignChange(ami - median(ami)));
out.pcrossq10 = mean(BF_SignChange(ami - quantile(ami, 0.1)));
out.pcrossq90 = mean(BF_SignChange(ami - quantile(ami, 0.9)));

% ac1
out.amiac1 = CO_AutoCorr(ami, 1, 'Fourier');

end
