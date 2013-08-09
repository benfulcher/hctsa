% CO_StickAngles
% 
% Analyzes line-of-sight angles between time-series points where each
% time-series value is treated as a stick protruding from an opaque baseline
% level. Statistics are returned on the raw time series, where sticks protrude
% from the zero-level, and the z-scored time series, where sticks
% protrude from the mean level of the time series.
% 
% INPUTS:
% y, the input time series
% 
% Outputs are returned on the obtained sequence of angles, theta, reflecting the
% maximum deviation a stick can rotate before hitting a stick representing
% another time point. Statistics include the mean and spread of theta,
% the different between positive and negative angles, measures of symmetry of
% the angles, stationarity, autocorrelation, and measures of the distribution of
% these stick angles.
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

function out = CO_StickAngles(y)
% Ben Fulcher, September 2009

doplot = 0; % can plot output
N = length(y);

ix = cell(2,1); %indicies for positive(1) and negative(2) entries of time series vector
ix{1} = find(y >= 0); % bias here -- 'look up' if on 'ground'
ix{2} = find(y < 0);

n = zeros(2,1);
n(1) = length(ix{1})-1; % minus one because the last point has no next one to compare to
n(2) = length(ix{2})-1; % minus one for same reason

angles = cell(2,1); % stores the angles
angles{1} = zeros(n(1),1); % positives (above axis)
angles{2} = zeros(n(2),1); % negatives (below axis)

% First positive time points: store in angles_p
for j = 1:2
    for i = 1:n(j)
        % find the next time series point with the same sign as the current one:
        angles{j}(i) = (y(ix{j}(i+1))-y(ix{j}(i)))/(ix{j}(i+1)-ix{j}(i));
    end
    angles{j} = atan(angles{j});
end
allangles = vertcat(angles{:});

if doplot
    figure('color','w')
    % A few options of what to plot:
    hold off; plot(angles{1},'.k'); hold on
    plot(angles{2},'.r'); hold off;
    [yp, xp] = ksdensity(angles{1});
    [yn, xn] = ksdensity(angles{2});
    plot(xp,yp,'r'); hold on; plot(xn,yn,'b');
    hist(angles{1},50);
end


% quite interesting
%% Basic stats?
% on raw values
out.std_p = std(angles{1});
out.mean_p = mean(angles{1});
out.median_p = median(angles{1});

out.std_n = std(angles{2});
out.mean_n = mean(angles{2});
out.median_n = median(angles{2});

out.std = std(allangles);
out.mean = mean(allangles);
out.median = median(allangles);

%% Difference between positive and negative angles
% Return difference in densities
ksx = linspace(min(allangles),max(allangles),200);
if ~isempty(angles{1}) && ~isempty(angles{2})
    ksy1 = ksdensity(angles{1},ksx); % spans the range of full extent (of both positive and negative angles)
    ksy2 = ksdensity(angles{2},ksx); % spans the range of full extent (of both positive and negative angles)
    out.pnsumabsdiff = sum(abs(ksy1-ksy2));
else
    out.pnsumabsdiff = NaN;
end


%% How symmetric is the distribution of angles?
% on raw outputs
% difference between ksdensities of positive and negative portions
if ~isempty(angles{1});
    maxdev=max(abs(angles{1}));
    ksy1 = ksdensity(angles{1},linspace(-maxdev,maxdev,201));
    out.symks_p = sum(abs(ksy1(1:100)-fliplr(ksy1(102:end))));
    out.ratmean_p = mean(angles{1}(angles{1}>0))/mean(angles{1}(angles{1}<0));
else
    out.symks_p = NaN; out.ratmean_p = NaN;
end

if ~isempty(angles{2})
    maxdev=max(abs(angles{2}));
    ksy2 = ksdensity(angles{2},linspace(-maxdev,maxdev,201));
    out.symks_n = sum(abs(ksy2(1:100)-fliplr(ksy2(102:end))));
    out.ratmean_n = mean(angles{2}(angles{2}>0))/mean(angles{2}(angles{2}<0));
else
    out.symks_n = NaN; out.ratmean_n = NaN;
end


%% z-score
zangles = cell(2,1);
zangles{1} = zscore(angles{1}); zangles{2} = zscore(angles{2}); zallangles = zscore(allangles);

%% how stationary are the angle sets?
% Do simple StatAvs for mean and std:
% ap_buff_2 = buffer(angles{1},floor(n(1)/2));
% if size(ap_buff_2,2)>2,ap_buff_2 = ap_buff_2(:,1:2);end % lose last point

% There are positive angles
if ~isempty(zangles{1});
    % StatAv2
    [statav_m statav_s] = sub_statav(zangles{1},2);
    out.statav2_p_m = statav_m;
    out.statav2_p_s = statav_s;
    
    % StatAv3
    [statav_m statav_s] = sub_statav(zangles{1},3);
    out.statav3_p_m = statav_m;
    out.statav3_p_s = statav_s;
    
    % StatAv 4
    [statav_m statav_s] = sub_statav(zangles{1},4);
    out.statav4_p_m = statav_m;
    out.statav4_p_s = statav_s;
    
    % StatAv 5
    [statav_m statav_s] = sub_statav(zangles{1},5);
    out.statav5_p_m = statav_m;
    out.statav5_p_s = statav_s;
    
else
    out.statav2_p_m = NaN; out.statav2_p_s = NaN;
    out.statav3_p_m = NaN; out.statav3_p_s = NaN;
    out.statav4_p_m = NaN; out.statav4_p_s = NaN;
    out.statav5_p_m = NaN; out.statav5_p_s = NaN;
end

% There are negative angles
if ~isempty(zangles{2});
    % StatAv2
    [statav_m statav_s] = sub_statav(zangles{2},2);
    out.statav2_n_m = statav_m;
    out.statav2_n_s = statav_s;
    
    % StatAv3
    [statav_m statav_s] = sub_statav(zangles{2},3);
    out.statav3_n_m = statav_m;
    out.statav3_n_s = statav_s;
    
    % StatAv4
    [statav_m statav_s] = sub_statav(zangles{2},4);
    out.statav4_n_m = statav_m;
    out.statav4_n_s = statav_s;
    
    % StatAv5
    [statav_m statav_s] = sub_statav(zangles{2},5);
    out.statav5_n_m = statav_m;
    out.statav5_n_s = statav_s;
else
    out.statav2_n_m = NaN; out.statav2_n_s = NaN;
    out.statav3_n_m = NaN; out.statav3_n_s = NaN;
    out.statav4_n_m = NaN; out.statav4_n_s = NaN;
    out.statav5_n_m = NaN; out.statav5_n_s = NaN;
end

% All angles:

% StatAv2
[statav_m, statav_s] = sub_statav(zallangles,2);
out.statav2_all_m = statav_m;
out.statav2_all_s = statav_s;

% StatAv3
[statav_m, statav_s] = sub_statav(zallangles,3);
out.statav3_all_m = statav_m;
out.statav3_all_s = statav_s;

% StatAv4
[statav_m, statav_s] = sub_statav(zallangles,4);
out.statav4_all_m = statav_m;
out.statav4_all_s = statav_s;

% StatAv5
[statav_m, statav_s] = sub_statav(zallangles,5);
out.statav5_all_m = statav_m;
out.statav5_all_s = statav_s;

%% correlations?
if ~isempty(zangles{1});
    out.tau_p = CO_FirstZero(zangles{1},'ac');
    out.ac1_p = CO_AutoCorr(zangles{1},1);
    out.ac2_p = CO_AutoCorr(zangles{1},2);
else
    out.tau_p = NaN; out.ac1_p = NaN; out.ac2_p = NaN;
end

if ~isempty(zangles{2});
    out.tau_n = CO_FirstZero(zangles{2},'ac');
    out.ac1_n = CO_AutoCorr(zangles{2},1);
    out.ac2_n = CO_AutoCorr(zangles{2},2);
else
    out.tau_n=NaN; out.ac1_n = NaN; out.ac2_n = NaN;
end

out.tau_all = CO_FirstZero(zallangles,'ac');
out.ac1_all = CO_AutoCorr(zallangles,1);
out.ac2_all = CO_AutoCorr(zallangles,2);

%% What does the distribution look like?
% Some quantiles and moments
if ~isempty(zangles{1});
    out.q1_p = quantile(zangles{1},0.01);
    out.q10_p = quantile(zangles{1},0.1);
    out.q90_p = quantile(zangles{1},0.9);
    out.q99_p = quantile(zangles{1},0.99);
    out.skewness_p = skewness(angles{1});
    out.kurtosis_p = kurtosis(angles{1});
else
    out.q1_p = NaN; out.q10_p = NaN;
    out.q90_p = NaN; out.q99_p = NaN;
    out.skewness_p = NaN; out.kurtosis_p = NaN;
end

if ~isempty(zangles{2});
    out.q1_n = quantile(zangles{2},0.01);
    out.q10_n = quantile(zangles{2},0.1);
    out.q90_n = quantile(zangles{2},0.9);
    out.q99_n = quantile(zangles{2},0.99);
    out.skewness_n = skewness(angles{2});
    out.kurtosis_n = kurtosis(angles{2});
else
    out.q1_n = NaN; out.q10_n = NaN;
    out.q90_n = NaN; out.q99_n = NaN;
    out.skewness_n = NaN; out.kurtosis_n = NaN;
end

F_quantz = @(x) quantile(zallangles,x);
out.q1_all = F_quantz(0.01);
out.q10_all = F_quantz(0.1);
out.q90_all = F_quantz(0.9);
out.q99_all = F_quantz(0.99);
out.skewness_all = skewness(allangles);
out.kurtosis_all = kurtosis(allangles);

%% Outliers?
% forget about this, I think.

function [statavmean, statavstd] = sub_statav(x,n)
    % Does a n-partition statav.
    % Require 2*n points (i.e., minimum of 2 in each partition) to do a
    % statav that even has a chance of being meaningful.
    NN = length(x);
    if NN < 2*n % not long enough
        statavmean = NaN; statavstd = NaN;
        return
    end
    
    x_buff = buffer(x,floor(NN/n));
    if size(x_buff,2) > n, x_buff = x_buff(:,1:n); end % lose last point
    statavmean = std(mean(x_buff))/std(x);
    statavstd = std(std(x_buff))/std(x);
end


end