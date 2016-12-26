function out = WL_coeffs(y, wname, level)
% WL_coeffs     Wavelet decomposition of the time series.
%
% Performs a wavelet decomposition of the time series using a given wavelet at a
% given level and returns a set of statistics on the coefficients obtained.
%
% Uses Matlab's Wavelet Toolbox.
%
%---INPUTS:
% y, the input time series
%
% wname, the wavelet name, e.g., 'db3' (see Wavelet Toolbox Documentation for
%                                       all options)
%
% level, the level of wavelet decomposition

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
BF_CheckToolbox('wavelet_toolbox')

% ------------------------------------------------------------------------------
%% Check Inputs
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

if wmaxlev(N,wname) < level
    fprintf(1,'Chosen level is too large for this wavelet on this signal\n');
    out = NaN;
    return
end

% ------------------------------------------------------------------------------
%% Perform a single-level wavelet decomposition
% (Recover a noisy signal by suppressing an approximation)

[c, l] = wavedec(y,level,wname);

% Reconstruct detail
det = wrcoef('d',c,l,wname,level); % detail this level

det_s = sort(abs(det),'descend'); % sorted detail coefficient magnitudes

% plot(det_s);

% ------------------------------------------------------------------------------
%% Return statistics
out.mean_coeff = mean(det_s);
out.max_coeff = max(det_s);
out.med_coeff = median(det_s);

% Decay rate stats ('where below _ maximum' = 'wb_m')
out.wb99m = findMyThreshold(0.99);
out.wb90m = findMyThreshold(0.90);
out.wb75m = findMyThreshold(0.75);
out.wb50m = findMyThreshold(0.50);
out.wb25m = findMyThreshold(0.25);
out.wb10m = findMyThreshold(0.10);
out.wb1m = findMyThreshold(0.01);

% ------------------------------------------------------------------------------
function propt = findMyThreshold(x)
    % where drops below proportion x of maximum
    propt = find(det_s < x*max(det_s),1,'first') / N;
            % (as a proportion of time-series length)
    if isempty(propt)
        propt = NaN;
    end
end
% ------------------------------------------------------------------------------

end
