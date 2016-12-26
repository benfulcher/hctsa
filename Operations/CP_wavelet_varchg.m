function out = CP_wavelet_varchg(y, wname, level, maxnchpts, minDelay)
% CP_wavelet_varchg     Variance change points in a time series.
%
% Finds variance change points using functions from Matlab's Wavelet Toolbox,
% including the primary function wvarchg, which estimates the change points in
% the time series.
%
%---INPUTS:
%
% y, the input time series
%
% wname, the name of the mother wavelet to analyze the data with: e.g., 'db3',
%           'sym2', cf. Wavelet Toolbox Documentation for details
%
% level, the level of wavelet decomposition
%
% maxnchpts, the maximum number of change points
%
% minDelay, the minimum delay between consecutive change points (can be
%           specified as a proportion of the time-series length, e.g., 0.02
%           ensures that change points are separated by at least 2% of the
%           time-series length)
%
%
%---OUTPUT:
% The optimal number of change points.

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

% Check that a Wavelet Toolbox license is available:
BF_CheckToolbox('wavelet_toolbox');

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

if nargin < 2 || isempty(wname)
    wname = 'db3'; % default wavelet
end

if nargin < 3 || isempty(level)
   level = 3; % level of wavelet decomposition
end
if strcmp(level,'max')
    level = wmaxlev(N,wname);
end

if nargin < 4 || isempty(maxnchpts)
   maxnchpts = 5; % maximum number of change points
end

if nargin < 5 || isempty(minDelay)
     minDelay = 0.01; % 1% of the time series length
end
if (minDelay > 0) && (minDelay < 1)
   minDelay = ceil(minDelay*N);
end

if wmaxlev(N, wname) < level
    error('Chosen level, %u, is too large for this wavelet on this signal. Sorry.', level);
end

% The aim of this example is to recover the change points in signal y.
% In addition, this example illustrates how the GUI tools propose change point
% locations for interval dependent de-noising thresholds.

% ------------------------------------------------------------------------------
% 1. Recover a noisy signal by suppressing an approximation.
% ------------------------------------------------------------------------------

% Perform a single-level wavelet decomposition :
[c, l] = wavedec(y,level,wname);

% Reconstruct detail at the same level.
det = wrcoef('d',c,l,wname,level);

% ------------------------------------------------------------------------------
% 2. Replace 2% of the greatest (absolute) values by the mean
% ------------------------------------------------------------------------------
% % in order to remove almost all the signal.
x = sort(abs(det));
v2p100 = x(fix(length(x)*0.98));
det(abs(det) > v2p100) = mean(det);

% ------------------------------------------------------------------------------
% 3. Use wvarchg to estimate the change points
% ------------------------------------------------------------------------------
try
    [~, kopt, ~] = wvarchg(det, maxnchpts, minDelay);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:nomem')
       error('Not enough memory.');
    end
end

% Return the number of change points found
out = kopt;

end
