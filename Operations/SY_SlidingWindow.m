function out = SY_SlidingWindow(y,windowStat,acrossWinStat,numSeg,incMove)
% SY_SlidingWindow  Sliding window measures of stationarity.
%
% This function is based on sliding a window along the time series, measuring
% some quantity in each window, and outputting some summary of this set of local
% estimates of that quantity.
%
% Another way of saying it: calculate 'windowStat' in each window, and computes
% 'acrossWinStat' for the set of statistics calculated in each window.
%
%---INPUTS:
%
% y, the input time series
%
% windowStat, the measure to calculate in each window:
%               (i) 'mean', mean
%               (ii) 'std', standard deviation
%               (iii) 'ent', distribution entropy
%               (iv) 'mom3', skewness
%               (v) 'mom4', kurtosis
%               (vi) 'mom5', the fifth moment of the distribution
%               (vii) 'lillie', the p-value for a Lilliefors Gaussianity test
%               (viii) 'AC1', the lag-1 autocorrelation
%               (ix) 'apen', Approximate Entropy
%
% acrossWinStat, controls how the obtained sequence of local estimates is
%                   compared (as a ratio to the full time series):
%                       (i) 'std': standard deviation
%                       (ii) 'ent' histogram entropy
%                       (iii) 'apen': Approximate Entropy, ApEn(1,0.2)
%                               cf. "Approximate entropy as a measure of system
%                               complexity", S. M. Pincus, P. Natl. Acad. Sci.
%                               USA 88(6) 2297 (1991)
%
% numSeg, the number of segments to divide the time series up into, thus
%       controlling the window length
%
% incMove, the increment to move the window at each iteration, as 1/fraction of the
%       window length (e.g., incMove = 2, means the window moves half the length of the
%       window at each increment)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

doPlot = 0; % plot outputs

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(windowStat)
    windowStat = 'mean'; % measure within each window
end
if nargin < 3 || isempty(acrossWinStat)
    acrossWinStat = 'std'; % measure across all windows
end
if nargin < 4 || isempty(numSeg)
    numSeg = 5;
end
if nargin < 5 || isempty(incMove)
    incMove = 2;
end

% ------------------------------------------------------------------------------

winLength = floor(length(y)/numSeg); % size of window
inc = floor(winLength/incMove); % increment to move at each step
if inc == 0; inc = 1; end % increment rounded down to zero, prop it up

numSteps = (floor((length(y)-winLength)/inc)+1);
qs = zeros(numSteps,1);

switch windowStat
    case 'mean' % Sliding window mean
        for i = 1:numSteps
            qs(i) = mean(y((i-1)*inc + 1:(i-1)*inc + winLength));
        end
    case 'std' % Sliding window std
        for i = 1:numSteps
            qs(i) = std(y((i-1)*inc + 1:(i-1)*inc + winLength));
        end
    case 'ent' % Sliding window distributional entropy
        for i = 1:numSteps
            ksstats = DN_FitKernelSmooth(y((i-1)*inc + 1:(i-1)*inc + winLength),'entropy');
            qs(i) = ksstats.entropy;
        end
    case 'apen' % Sliding window ApEn
        for i = 1:numSteps
            qs(i) = EN_ApEn(y((i-1)*inc + 1:(i-1)*inc + winLength),1,0.2);
        end
    case 'mom3' % Third moment
        for i = 1:numSteps
            qs(i) = DN_Moments(y((i-1)*inc + 1:(i-1)*inc + winLength),3);
        end
    case 'mom4' % Fourth moment
        for i = 1:numSteps
            qs(i) = DN_Moments(y((i-1)*inc + 1:(i-1)*inc + winLength),4);
        end
    case 'mom5' % Fifth moment
        for i = 1:numSteps
            qs(i) = DN_Moments(y((i-1)*inc + 1:(i-1)*inc + winLength),5);
        end
    case 'lillie' % Lilliefors test
        for i = 1:numSteps
            qs(i) = HT_DistributionTest(y((i-1)*inc + 1:(i-1)*inc + winLength),'lillie','norm');
        end
    case 'AC1' % Lag-1 autocorrelation
        for i = 1:numSteps
            qs(i) = CO_AutoCorr(y((i-1)*inc + 1:(i-1)*inc + winLength),1,'Fourier');
        end
    otherwise
        error('Unknown statistic ''%s''',windowStat)
end

% ------------------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------------------
if doPlot
    figure('color','w'); box('on');
    plot(round(winLength/2):inc:(numSteps-1)*inc+round(winLength/2),qs,'r');
end

% ------------------------------------------------------------------------------
% Compute the output statistic
% ------------------------------------------------------------------------------
switch acrossWinStat
    case 'std'
        out = std(qs)/std(y);
    case 'apen'
        out = EN_ApEn(qs,1,0.2); % ApEn of the sliding window measures
    case 'ent'
        kssimpouts = DN_FitKernelSmooth(qs); % get a load of statistics from kernel-smoothed distribution
        out = kssimpouts.entropy; % distributional entropy
    otherwise
        error('Unknown statistic: ''%s''',acrossWinStat)
end

end
