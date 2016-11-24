function out = WL_scal2frq(y, wname, amax, delta)
% WL_scal2frq   Frequency components in a periodic time series
%
% Estimates frequency components using functions from Matlab's Wavelet Toolbox,
% including the scal2frq function.
%
%---INPUTS:
% y, the input time series
%
% wname, the name of the mother wavelet to analyze the data with: e.g., 'db3',
%           'sym2', cf. Wavelet Toolbox Documentation for details
%
% amax, the maximum scale / level (can be 'max' to set according to wmaxlev)
%
% delta, the sampling period
%
%---OUTPUTS: the level with the highest energy coefficients, the dominant
% period, and the dominant pseudo-frequency.
%
% Adapted from example in Matlab Wavelet Toolbox documentation. It's kind of a
% weird idea to apply the method to generic time series.

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
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox')

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
doplot = 0; % plot outputs to figure
N = length(y); % length of the time series

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(wname)
    fprintf(1,'Wavelet not specified -- using the default db3 wavelet\n');
    wname = 'db3';
end

if nargin < 3 || isempty(amax)
    amax = 5; % maximum 'scale'
end
maxlevel = wmaxlev(N,wname); % maximum level for this time-series length
if strcmp(amax,'max') % set to maximum for this wavelet
    amax = wmaxlev(N,wname);
end

if nargin < 4 || isempty(delta)
    delta = 1; % the sampling period
end

if maxlevel < amax
    fprintf(1,'Chosen level (%u) is too large for this wavelet on this signal...',amax);
    amax = maxlevel;
    fprintf(1,' changed to maximum level computed with wmaxlev: %u\n',amax);
end

% ------------------------------------------------------------------------------
%% Do your thing:
% ------------------------------------------------------------------------------
% This example demonstrates that, starting from the periodic function
% x(t) = 5*sin(5t) + 3*sin(2t) + 2*sin(t), the scal2frq function translates
% the scales corresponding to the maximum values of the CWT coefficients
% to pseudo-frequencies ([0.796 0.318 0.159]), which are near to the true
% frequencies ([5, 2, 1] / (2*pi) =~ [0.796 0.318 0.159]).

% delta = 0.1;
% wname = 'coif3';

% Define scales.
scales = 1:amax;
a = 2.^scales;

% Compute associated pseudo-frequencies.
f = scal2frq(a, wname, delta);

% Compute associated pseudo-periods.
per = 1./f;

% Decompose the time series at level specified as maximum
[c, l] = wavedec(y, amax, wname);

% Estimate standard deviation of detail coefficients.
stdc = wnoisest(c,l,scales);

if doplot
    figure('color','w'); box('on');
    plot(stdc,'k') % plot them
end

% Compute identified period.
[~, jmax] = max(stdc); % level with highest energy coefficients

out.lmax = jmax; % level with highest energy coefficients
out.period = per(jmax); % output dominant period
out.pf = f(jmax); % output dominant pseudo-frequency

end
