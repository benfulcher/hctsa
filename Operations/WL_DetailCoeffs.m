function out = WL_DetailCoeffs(y, wname, maxlevel)
% WL_DetailCoeffs   Detail coefficients of a wavelet decomposition.
%
% Compares the detail coefficients obtained at each level of the wavelet
% decomposition from 1 to the maximum possible level for the wavelet given the
% length of the input time series (computed using wmaxlev from
% Matlab's Wavelet Toolbox).
%
%---INPUTS:
% y, the input time series
%
% wname, the name of the mother wavelet to analyze the data with: e.g., 'db3',
%           'sym2', cf. Wavelet Toolbox Documentation for details
%
% maxlevel, the maximum wavelet decomposition level (can also set to 'max' to be
%               that determined by wmaxlev)
%
%---OUTPUTS:
% Statistics on the detail coefficients.

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
BF_CheckToolbox('wavelet_toolbox');
% ------------------------------------------------------------------------------

doPlot = 0; % plot outputs to figure

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

if nargin < 2 || isempty(wname)
    wname = 'db3'; % default wavelet
end
if nargin < 3 || isempty(maxlevel)
   maxlevel = 20; % maximum wavelet decomposition level
end
if strcmp(maxlevel,'max')
    maxlevel = wmaxlev(N,wname);
end

if wmaxlev(N, wname) < maxlevel
    fprintf(1,'Chosen wavelet level is too large for the %s wavelet for this signal of length N = %u\n',wname,N);
    maxlevel = wmaxlev(N,wname);
    fprintf(1,'Using a wavelet level of %u instead.\n',maxlevel);
end

% ------------------------------------------------------------------------------
%% Perform a single-level wavelet decomposition
% ------------------------------------------------------------------------------
means = zeros(maxlevel,1); % mean detail coefficient magnitude at each level
medians = zeros(maxlevel,1); % median detail coefficient magnitude at each level
maxs = zeros(maxlevel,1); % max detail coefficient magnitude at each level

for k = 1:maxlevel
    level = k;

    [c, l] = wavedec(y,level,wname);
    % Reconstruct detail at this level
    det = wrcoef('d',c,l,wname,level);

    means(k) = mean(abs(det));
    medians(k) = median(abs(det));
    maxs(k) = max(abs(det));
end

% ------------------------------------------------------------------------------
%% Plot the bad boy
% ------------------------------------------------------------------------------
if doPlot
    subplot(5,1,1:2); title('signal')
    plot(y);
    subplot(5,1,3); title('means');
    plot(means)
    subplot(5,1,4); title('medians');
    plot(medians)
    subplot(5,1,5); title('maxs');
    plot(maxs);
end

% ------------------------------------------------------------------------------
%% Return statistics on detail coefficients
% ------------------------------------------------------------------------------
% Sort
means_s = sort(means,'descend');
medians_s = sort(medians,'descend');
maxs_s = sort(maxs,'descend');

% What is the maximum across these levels
out.max_mean = means_s(1);
out.max_median = medians_s(1);
out.max_max = maxs_s(1);

% stds
out.std_mean = std(means);
out.std_median = std(medians);
out.std_max = std(maxs);

% At what level is the maximum
out.wheremax_mean = find(means == means_s(1),1,'first');
out.wheremax_median = find(medians == medians_s(1),1,'first');
out.wheremax_max = find(maxs == maxs_s(1),1,'first');

% Size of maximum (relative to next maximum)
out.max1on2_mean = means_s(1)/means_s(2);
out.max1on2_median = medians_s(1)/medians_s(2);
out.max1on2_max = maxs_s(1)/maxs_s(2);

% Where sum of values to left equals sum of values to right
% Measure of centrality
out.wslesr_mean = SUB_slosr(means);
out.wslesr_median = SUB_slosr(medians);
out.wslesr_max = SUB_slosr(maxs);

% What's the correlation between maximum and median
r = corrcoef(maxs,medians);
out.corrcoef_max_medians = r(1,2);

% ------------------------------------------------------------------------------
function meIsGorilla = SUB_slosr(xx)
    theMaxLevel = length(xx);
    slosr = zeros(theMaxLevel-2,1);
    for i = 2:theMaxLevel-1
        slosr(i-1) = sum(xx(1:i-1))/sum(xx(i+1:end));
    end
    absm1 = abs(slosr-1); % how close to 1 (the same sum on either side) each is
    meIsGorilla = find(absm1 == min(absm1),1,'first') + 1;
end
% ------------------------------------------------------------------------------

end
