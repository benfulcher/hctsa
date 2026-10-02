function out = SY_KPSStest(y, lags)
% SY_KPSStest   The KPSS stationarity test.
%
% The KPSS stationarity test, of Kwiatkowski, Phillips, Schmidt, and Shin, using
% the function kpsstest from MATLAB's Econometrics Toolbox. The null hypothesis is
% that a univariate time series is trend stationary; the alternative hypothesis
% is that it is a non-stationary unit-root process. Large test statistics and
% small p-values reject stationarity.
%
% For a scalar lag, the test statistic and p-value are returned. If the input is
% a vector of lags, statistics on how the test statistics and p-values change
% across these lags are returned instead.
%
% ---INPUTS:
% y, the input time series
%
% lags, the number of autocovariance lags used in the long-run variance estimate.
%       Either a scalar (returns the test statistic and p-value) or a vector
%       (returns statistics on changes across these lags). Default: 0.
%
% ---OUTPUTS:
% For scalar lags:
% stat, the KPSS test statistic
% pValue, the p-value of the test
% For vector lags:
% maxpValue, minpValue: the maximum and minimum p-value across lags
% maxstat, minstat: the maximum and minimum test statistic across lags
% lagmaxstat, lagminstat: the lag at which the test statistic is largest and
%       smallest (the first such lag in the case of a tie)
%
% ---REFERENCES:
% Kwiatkowski, Denis and Phillips, Peter C. B. and Schmidt, Peter and Shin,
% Yongcheol, "Testing the null hypothesis of stationarity against the alternative
% of a unit root: How sure are we that economic time series have a unit root?",
% J. Econometrics 54(1-3), 159-178 (1992). DOI: 10.1016/0304-4076(92)90104-Y

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
%% Check that an Econometrics license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('econometrics_toolbox');

% -------------------------------------------------------------------------------
% Check inputs
% -------------------------------------------------------------------------------
if nargin < 2 || isempty(lags)
	lags = 0;
end

% ------------------------------------------------------------------------------
%% (1) Perform the test(s)
% ------------------------------------------------------------------------------
% Temporarily turn off warnings for the test statistic being too big or small
warning('off', 'econ:kpsstest:StatTooSmall')
warning('off', 'econ:kpsstest:StatTooBig')
[~, pValue, stat] = kpsstest(y, 'lags', lags);
warning('on', 'econ:kpsstest:StatTooSmall')
warning('on', 'econ:kpsstest:StatTooBig')

% ------------------------------------------------------------------------------
%% (2) Return statistics on outputs of test(s)
% ------------------------------------------------------------------------------
if length(lags) > 1
	% Return statistics on outputs
	out.maxpValue = max(pValue);
	out.minpValue = min(pValue);
	out.maxstat = max(stat);
	out.minstat = min(stat);
	% (index with the FIRST extremum: an exact tie between two lags would
	% otherwise make these fields vectors rather than scalars)
	out.lagmaxstat = lags(find(stat == max(stat), 1, 'first')); % lag at max test statistic
	out.lagminstat = lags(find(stat == min(stat), 1, 'first'));
else
	% return the statistic and pvalue
	out.stat = stat;
	out.pValue = pValue;
end

end
