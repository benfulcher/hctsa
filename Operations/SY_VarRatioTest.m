function out = SY_VarRatioTest(y,periods,IIDs)
% SY_VarRatioTest   Variance ratio test for random walk.
%
% Implemented using the vratiotest function from Matlab's Econometrics Toolbox.
%
% The test assesses the null hypothesis of a random walk in the time series,
% which is rejected for some critical p-value.
%
%---INPUTS:
% y, the input time series
%
% periods, a vector (or scalar) of period(s)
%
% IIDs, a vector (or scalar) representing boolean values indicating whether to
%       assume independent and identically distributed (IID) innovations for
%       each period.

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
%% Check that an Econometrics Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('econometrics_toolbox')

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
% Can set step sizes for random walk, and also change the null hypothesis
% to include non IID random walk increments

% periods, e.g., could be [2,4,6,8,2,4,6,8]
if nargin < 2 || isempty(periods)
    periods = 2;
end

% IIDs, e.g., could be [1,1,1,1,0,0,0,0]
if nargin < 3 || isempty(IIDs)
    IIDs = 0;
end
IIDs = logical(IIDs);

% ------------------------------------------------------------------------------
%% Perform the test:
% ------------------------------------------------------------------------------
[h, pValue, stat, ~, ratio] = vratiotest(y,'period',periods,'IID',IIDs);

if length(h) == 1
    % Summarize the single test performed
    out.pValue = pValue;
    out.stat = stat;
    out.ratio = ratio;

else
   % Return statistics on multiple outputs for multiple periods/IIDs
   out.maxpValue = max(pValue);
   out.minpValue = min(pValue);
   out.meanpValue = mean(pValue);

   imaxp = find(pValue == max(pValue),1,'first');
   iminp = find(pValue == min(pValue),1,'first');
   out.periodmaxpValue = periods(imaxp);
   out.periodminpValue = periods(iminp);
   out.IIDperiodmaxpValue = IIDs(imaxp);
   out.IIDperiodminpValue = IIDs(iminp);

   out.meanstat = mean(stat);
   out.maxstat = max(stat);
   out.minstat = min(stat);
end

end
