function out = DT_IsSeasonal(y)
% DT_IsSeasonal     A simple test of seasonality.
%
% Fits a 'sin1' model to the time series using fit function from the Curve Fitting
% Toolbox. The output is binary: 1 if the goodness of fit, R^2, exceeds 0.3 and
% the amplitude of the fitted periodic component exceeds 0.5, and 0 otherwise.
%
%---INPUTS:
% y, the input time series
%
%---OUTPUT: Binary: 1 (= seasonal), 0 (= non-seasonal)

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

%-------------------------------------------------------------------------------
%% Preliminaries
%-------------------------------------------------------------------------------
% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

% Make sure the input time series, y, is a column vector
if size(y,2) > size(y,1)
    y = y';
end
N = length(y); % length of input time series
r = (1:N)'; % range over which to fit

%-------------------------------------------------------------------------------
%% Fit a sinusoidal model using the Curve-Fitting Toolbox
%-------------------------------------------------------------------------------
[cfun, gof] = fit(r,y,'sin1'); % fits the following form: a1*sin(b1*x+c1)

%-------------------------------------------------------------------------------
%% Two conditions for determining whether time series contains periodicities:
%-------------------------------------------------------------------------------
% Condition 1: fit is ok
th_fit = 0.3; % r2>th_fit

% Condition 2: amplitude is not too small
th_ampl = 0.5; % a1>th_ampl

if gof.rsquare > th_fit && abs(cfun.a1 > th_ampl)
    out = 1; % test thinks the time series has strong periodicities
else
    out = 0; % test thinks the time series doesn't have any strong periodicities
end

end
