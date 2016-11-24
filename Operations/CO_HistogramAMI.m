function out = CO_HistogramAMI(y,tau,meth,numBins)
% CO_HistogramAMI      The automutual information of the distribution using histograms.
%
% The approach used to bin the data is provided.
%
%---INPUTS:
%
% y, the input time series
%
% tau, the time-lag (1 by default)
%
% meth, the method of computing automutual information:
%           (i) 'even': evenly-spaced bins through the range of the time series,
%           (ii) 'std1', 'std2': bins that extend only up to a multiple of the
%                                standard deviation from the mean of the time
%                                series to exclude outliers,
%           (iii) 'quantiles': equiprobable bins chosen using quantiles.
%
% numBins, the number of bins, required by some methods, meth (see above)
%
%---OUTPUT: the automutual information calculated in this way.

% Uses the hist2 function (renamed NK_hist2.m here) by Nedialko Krouchev, obtained
% from Matlab Central,
% http://www.mathworks.com/matlabcentral/fileexchange/12346-hist2-for-the-people
% [[hist2 for the people by Nedialko Krouchev, 20 Sep 2006 (Updated 21 Sep 2006)]]

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
% Time-lag, tau
if nargin < 2 || isempty(tau)
    tau = 1;  % time-lag of 1
end
if ischar(tau) && ismember(tau,{'ac','tau'})
    tau = CO_FirstZero(y,'ac');
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
%% Bins by standard deviation (=1)
% ------------------------------------------------------------------------------

% same for both -- assume same distribution (true for stationary processes,
% or small lags)
switch meth
    case 'even'
        b = linspace(min(y)-0.1,max(y)+0.1,numBins+1); % +0.1 to make sure all points included

    case 'std1' % std bins up to 1
        b = linspace(-1,1,numBins+1);
        if min(y) < -1; b = [min(y)-0.1, b]; end
        if max(y) > 1; b = [b, max(y)+0.1]; end

    case 'std2'
        b = linspace(-2,2,numBins+1);
        if min(y) < -2; b = [min(y)-0.1, b]; end
        if max(y) > 2; b = [b, max(y)+0.1]; end

    case 'quantiles' % use quantiles with ~equal number in each bin
        b = quantile(y,linspace(0,1,numBins+1));
        b(1) = b(1) - 0.1; b(end) = b(end) + 0.1;

    otherwise
        error('Unknown method ''%s''',meth)
end

nb = length(b) - 1; % number of bins (-1 since b defines edges)

% ------------------------------------------------------------------------------
% Form the time-delay vectors y1 and y2
% ------------------------------------------------------------------------------
amis = zeros(length(tau),1);
for i = 1:length(tau)
    y1 = y(1:end-tau(i));
    y2 = y(1+tau(i):end);

    % (1) Joint distribution of y1 and y2
    pij = NK_hist2(y1,y2,b,b);
    pij = pij(1:nb,1:nb); % joint
    pij = pij/sum(sum(pij)); % joint
    pi = sum(pij,1); % marginal
    pj = sum(pij,2); % other marginal

    % Old-fashioned method (should give same result):
    % pi = histc(y1,b); pi = pi(1:nb); pi = pi/sum(pi); % marginal
    % pj = histc(y2,b); pj= pj(1:nb); pj = pj/sum(pj); % other marginal

    pii = ones(nb,1)*pi;
    pjj = pj*ones(1,nb);

    r = (pij > 0); % Defining the range in this way, we set log(0) = 0
    amis(i) = sum(pij(r).*log(pij(r)./pii(r)./pjj(r)));
end

if length(tau)==1
    out = amis;
else
    for i = 1:length(tau)
        out.(sprintf('ami%u',i)) = amis(i);
    end
end

end
