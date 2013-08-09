% WL_coeffs
% 
% Performs a wavelet decomposition of the time series using a given wavelet at a
% given level and returns a set of statistics on the coefficients obtained.
% 
% Uses Matlab's Wavelet Toolbox.
% 
% INPUTS:
% y, the input time series
% 
% wname, the wavelet name, e.g., 'db3' (see Wavelet Toolbox Documentation for
%                                       all options)
% 
% level, the level of wavelet decomposition
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = WL_coeffs(y, wname, level)
% Ben Fulcher, 23/1/2010

%% Check that a Wavelet Toolbox license exists:
a = license('test','wavelet_toolbox');
if a==0
    error('This function requires Matlab''s Wavelet Toolbox');
end
% Try to check out a license:
[lic_free,~] = license('checkout','wavelet_toolbox');
if lic_free == 0
    error('Could not obtain a license for Matlab''s Wavelet Toolbox');
end

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
    out = NaN; return
end

% 1. Recover a noisy signal by suppressing an
% approximation.

%% Perform a single-level wavelet decomposition 
[c, l] = wavedec(y,level,wname);

% Reconstruct detail
det = wrcoef('d',c,l,wname,level); % detail this level

det_s = sort(abs(det),'descend'); % sorted detail coefficient magnitudes

% plot(det_s);

%% Return statistics
out.mean_coeff = mean(det_s);
out.max_coeff = max(det_s);
out.med_coeff = median(det_s);

% decay rate stats ('where below _ maximum' = 'wb_m')
out.wb99m = findmythreshold(0.99);
out.wb90m = findmythreshold(0.90);
out.wb75m = findmythreshold(0.75);
out.wb50m = findmythreshold(0.50);
out.wb25m = findmythreshold(0.25);
out.wb10m = findmythreshold(0.10);
out.wb1m = findmythreshold(0.01);

    function propt = findmythreshold(x)
        % where drops below proportion x of maximum
        propt = find(det_s < x*max(det_s),1,'first') / N;
                % (as a proportion of time-series length)
        if isempty(i), propt = NaN; end
    end

end