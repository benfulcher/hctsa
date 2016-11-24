function out = SY_Trend(y)
% SY_Trend  Quantifies various measures of trend in a time series.
%
%---INPUT:
% y, the input time series.
%
%---OUTPUTS:
% Linearly detrends the time series using detrend, and returns the ratio of
% standard deviations before and after the linear detrending. If a strong linear
% trend is present in the time series, this operation should output a low value.
%
% Also fits a line and gives parameters from that fit, as well as statistics on
% a cumulative sum of the time series.

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

if ~BF_iszscored(y)
    warning('The input time series should be z-scored')
end

N = length(y);

% Ratio of std before and after linear detrending:
out.stdRatio = std(detrend(y)) / std(y);

% Linear fit:
[out.gradient,out.intercept] = LinearFit(1:N,y);

% Stats on the cumulative sum:
yC = cumsum(y);
out.meanYC = mean(yC);
out.stdYC = std(yC);
[out.gradientYC,out.interceptYC] = LinearFit(1:N,yC);

% Mean cumsum in first and second half of the time series:
out.meanYC12 = mean(yC(1:floor(N/2)));
out.meanYC22 = mean(yC(floor(N/2)+1:end));

% ------------------------------------------------------------------------------
function [m,b] = LinearFit(xData,yData)
    if size(xData,1) ~= size(yData,1);
        yData = yData';
    end
    coeff = polyfit(xData,yData,1);
    m = coeff(1); b = coeff(2);
end
% ------------------------------------------------------------------------------

end
