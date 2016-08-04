function out = SY_SpreadRandomLocal(y,l,numSegs,randomSeed)
% SY_SpreadRandomLocal  Bootstrap-based stationarity measure.
%
% numSegs time-series segments of length l are selected at random from the time
% series and in each segment some statistic is calculated: mean, standard
% deviation, skewness, kurtosis, ApEn(1,0.2), SampEn(1,0.2), AC(1), AC(2), and the
% first zero-crossing of the autocorrelation function.
% Outputs summarize how these quantities vary in different local segments of the
% time series.
%
%---INPUTS:
% y, the input time series
%
% l, the length of local time-series segments to analyze as a positive integer.
%    Can also be a specified character string:
%       (i) 'ac2': twice the first zero-crossing of the autocorrelation function
%       (ii) 'ac5': five times the first zero-crossing of the autocorrelation function
%
% numSegs, the number of randomly-selected local segments to analyze
%
% randomSeed, the input to BF_ResetSeed to control reproducibility
%
%---OUTPUTS: the mean and also the standard deviation of this set of 100 local
% estimates.

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

doPlot = 0; % set to 1 to plot outputs to figure

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(l)
    l = 100; % by default use 100 samples
end

if ischar(l)
    taug = CO_FirstZero(y,'ac'); % tau (global)
    switch l
    case 'ac2'
        l = 2*taug;
    case 'ac5'
        l = 5*taug;
    otherwise
        error('Unknown specifier ''%s''',l);
    end

    % Very short l for this sort of time series:
    if l < 5
        warning(['This time series has a very short correlation length;\nSetting ' ...
            'l=%u means that changes estimates will be difficult to compare...'],l);
    end
end

if nargin < 3 || isempty(numSegs)
    numSegs = 100; % 100 segments by default
end

% Check the parameters are appropriate for the length of the input time series:
N = length(y); % the length of the time series
if l > 0.9*N % operation is not suitable -- time series is too short
	warning(['SY_SpreadRandomLocal: This time series (N = %u) ' ...
                                'is too short to use l = %.1f\n'],N,l)
    out = NaN; return % NaN means not suitable
end

if nargin < 4
    randomSeed = []; % use default random seed
end

% ------------------------------------------------------------------------------
% numSegs segments, each of length segl data points

numFeat = 8; % number of features
qs = zeros(numSegs,numFeat);

% Reset random seed, for reproducibility:
BF_ResetSeed(randomSeed);

for j = 1:numSegs
    % pick a range
    % in this implementation, ranges CAN overlap

    ist = randi(N-1-l,1); % random start point (not exceeding the endpoint)
    ifh = ist+l-1; % finish index
    rs = ist:ifh; % sample range (from starting to finishing index)
    ysub = y(rs); % subsection of the time series

    taul = CO_FirstZero(ysub,'ac');

    qs(j,1) = mean(ysub); % mean
    qs(j,2) = std(ysub); % standard deviation
    qs(j,3) = skewness(ysub); % skewness
    qs(j,4) = kurtosis(ysub); % kurtosis
    entropyStruct = EN_SampEn(ysub,1,0.15);
    qs(j,5) = entropyStruct.quadSampEn1; % SampEn_1_01
    qs(j,6) = CO_AutoCorr(ysub,1,'Fourier'); % AC1
    qs(j,7) = CO_AutoCorr(ysub,2,'Fourier'); % AC2
    qs(j,8) = taul;
end

% ------------------------------------------------------------------------------
% Plot some output?:
% ------------------------------------------------------------------------------
if doPlot
    figure('color','w');
    subplot(2,1,1); hold on;
    plot(y,'k');
    plot(ists,y(ists),'.r');
    title('time series')
    subplot(2,1,2); plot(qs(:,1),'b'); title('local means')
end

% ------------------------------------------------------------------------------
% Take mean and standard deviation of this set of local time-series statistics:
% ------------------------------------------------------------------------------
% Can think of this as a big bootstrapped distribution of the timeseries at
% a scale given by the length l

fs = zeros(numFeat,2);
fs(:,1) = nanmean(qs); % the mean value of the feature across subsegments of the time series
fs(:,2) = nanstd(qs); % the spread of the feature across subsegments of the time series

out.meanmean = fs(1,1);
out.meanstd = fs(2,1);
out.meanskew = fs(3,1);
out.meankurt = fs(4,1);
out.meansampen1_015 = fs(5,1);
out.meanac1 = fs(6,1);
out.meanac2 = fs(7,1);
out.meantaul = fs(8,1);

out.stdmean = fs(1,2);
out.stdstd = fs(2,2);
out.stdskew = fs(3,2);
out.stdkurt = fs(4,2);
out.stdsampen1_015 = fs(5,2);
out.stdac1 = fs(6,2);
out.stdac2 = fs(7,2);
out.stdtaul = fs(8,2);

end
