function out = ST_SimpleStats(x,whatStat)
% ST_SimpleStats   Basic statistics about an input time series
%
%---INPUTS:
% x, the input time series
%
% whatStat, the statistic to return:
%          (i) 'zcross': the proportionof zero-crossings of the time series
%                        (z-scored input thus returns mean-crossings)
%          (ii) 'maxima': the proportion of the time series that is a local maximum
%          (iii) 'minima': the proportion of the time series that is a local minimum
%          (iv) 'pmcross': the ratio of the number of times that the (ideally
%                          z-scored) time-series crosses +1 (i.e., 1 standard
%                          deviation above the mean) to the number of times
%                          that it crosses -1 (i.e., 1 standard deviation below
%                          the mean)
%          (v) 'zsczcross': the ratio of zero crossings of raw to detrended
%                           time series where the raw has zero mean

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

N = length(x); % length of the time series

switch whatStat
    case 'zcross'
        % Proportion of zero-crossings of the time series
        % (in the case of z-scored input, crosses its mean)
        xch = x(1:end-1).*x(2:end);
        out = sum(xch < 0)/N;

    case 'maxima'
        % Proportion of local maxima in the time series
        dx = diff(x);
        out = sum(dx(1:end-1) > 0 & dx(2:end) < 0)/(N-1);

    case 'minima'
        % Proportion of local minima in the time series
        dx = diff(x);
        out = sum(dx(1:end-1) < 0 & dx(2:end) > 0)/(N-1);

    case 'pmcross'
        % ratio of times cross 1 to -1
        c1sig = sum(BF_sgnchange(x-1)); % num times cross 1
        c2sig = sum(BF_sgnchange(x+1)); % num times cross -1
        if c2sig == 0
            out = NaN;
        else
            out = c1sig/c2sig;
        end

    case 'zsczcross'
        % ratio of zero crossings of raw to detrended time series
        % where the raw has zero mean
        x = zscore(x);
        xch = x(1:end-1).*x(2:end);
        h1 = sum(xch < 0); % num of zscross of raw series
        y = detrend(x);
        ych = y(1:end-1).*y(2:end);
        h2 = sum(ych < 0); % of detrended series
        if h1 == 0;
            out = NaN;
        else
            out = h2/h1;
        end

    otherwise
        error('Unknown statistic ''%s''',whatStat);
end

end
