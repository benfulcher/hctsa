function out = SY_KPSStest(y,lags)
% SY_KPSStest   The KPSS stationarity test.
%
% The KPSS stationarity test, of Kwiatkowski, Phillips, Schmidt, and Shin:
% "Testing the null hypothesis of stationarity against the alternative of a
% unit root: How sure are we that economic time series have a unit root?"
% Kwiatkowski, Denis and Phillips, Peter C. B. and Schmidt, Peter and Shin, Yongcheol
% J. Econometrics, 54(1-3) 159 (2002)
%
% Uses the function kpsstest from Matlab's Econometrics Toolbox. The null
% hypothesis is that a univariate time series is trend stationary, the
% alternative hypothesis is that it is a non-stationary unit-root process.
%
% The code can implemented for a specific time lag, tau. Alternatively, measures
% of change in p-values and test statistics will be outputted if the input is a
% vector of time lags.
%
%---INPUTS:
% y, the input time series
% lags, can be either a scalar (returns basic test statistic and p-value), or
%                   vector (returns statistics on changes across these time lags)

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
%% Check that an Econometrics license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('econometrics_toolbox')

%-------------------------------------------------------------------------------
% Check inputs
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(lags)
    lags = 0;
end

% ------------------------------------------------------------------------------
%% (1) Perform the test(s)
% ------------------------------------------------------------------------------
% Temporarily turn off warnings for the test statistic being too big or small
warning('off','econ:kpsstest:StatTooSmall')
warning('off','econ:kpsstest:StatTooBig')
[~, pValue, stat] = kpsstest(y,'lags',lags);
warning('on','econ:kpsstest:StatTooSmall')
warning('on','econ:kpsstest:StatTooBig')

% ------------------------------------------------------------------------------
%% (2) Return statistics on outputs of test(s)
% ------------------------------------------------------------------------------
if length(lags) > 1
    % Return statistics on outputs
    out.maxpValue = max(pValue);
    out.minpValue = min(pValue);
    out.maxstat = max(stat);
    out.minstat = min(stat);
    out.lagmaxstat = lags(stat == max(stat)); % lag at max test statistic
    out.lagminstat = lags(stat == min(stat));
else
    % return the statistic and pvalue
    out.stat = stat;
    out.pValue = pValue;
end

end
