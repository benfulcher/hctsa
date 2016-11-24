function out = CO_Embed2(y,tau)
% CO_Embed2     2-dimensional embedding space
%
% Embeds the z-scored time series in a two-dimensional time-delay
% embedding space with a given time-delay, tau, and outputs a set of
% statistics about the structure in this space, including angular
% distribution, etc.
%
%---INPUTS:
% y, the column-vector time series
% tau, the time-delay (can be 'tau' for first zero-crossing of ACF)
%
%---OUTPUTS: include the distribution of angles between successive points in the
% space, stationarity of this angular distribution, euclidean distances from the
% origin, and statistics on outliers.

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


doPlot = 0; % can set to 1 to plot some outputs

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 'tau';
end

% Set tau to the first zero-crossing of the autocorrelation function, with the 'tau' input
if strcmp(tau,'tau'),
    tau = CO_FirstZero(y,'ac');
    if tau > length(y)/10
        tau = floor(length(y)/10);
    end
end

% Ensure that y is a column vector
if size(y,2) > size(y,1);
    y = y';
end

% Construct the two-dimensional recurrence space
m = [y(1:end-tau), y(1+tau:end)];
N = size(m,1); % number of points in the recurrence space

if doPlot
    figure('color','w');
    plot(m(:,1),m(:,2),'.');
end

% 1) Distribution of angles time series; angles between successive points in
% 	 this space

theta = diff(m(:,2))./diff(m(:,1));
theta = atan(theta); % measured as deviation from the horizontal


if doPlot, ksdensity(theta); end % can plot distribution of angles
out.theta_ac1 = CO_AutoCorr(theta,1,'Fourier');
out.theta_ac2 = CO_AutoCorr(theta,2,'Fourier');
out.theta_ac3 = CO_AutoCorr(theta,3,'Fourier');

out.theta_mean = mean(theta);
out.theta_std = std(theta);

binEdges = linspace(-pi/2,pi/2,11); % 10 bins in the histogram
[px,binEdges] = histcounts(theta,binEdges,'Normalization','probability');
binWidths = diff(binEdges);
out.hist10std = std(px);
out.histent = -sum(px(px>0).*log(px(px>0)./binWidths(px>0)));

%-------------------------------------------------------------------------------
% Stationarity in fifths of the time series
%-------------------------------------------------------------------------------
% Use histograms with 4 bins
x = linspace(-pi/2,pi/2,5); % 4 bins
afifth = floor((N-1)/5); % -1 because angles are correlations *between* points
n = zeros(length(x),5);
for i = 1:5
	n(:,i) = histc(theta(afifth*(i-1)+1:afifth*i),x);
end
n = n/afifth;
n(4,:) = n(4,:) + n(5,:);
n(5,:) = [];

% Output the standard deviation in each bin:
for i = 1:4
    out.(sprintf('stdb%u',i)) = std(n(:,i));
end

%-------------------------------------------------------------------------------
% STATIONARITY of points in the space (do they move around in the space)
%-------------------------------------------------------------------------------

% (1) in terms of distance from origin
afifth = floor(N/5);
buffer_m = cell(5,1); % stores a fifth of the time series (embedding vector) in each entry
for i = 1:5
    buffer_m{i} = m(afifth*(i-1)+1:afifth*i,:);
end

% Mean euclidean distance in each segment
eucdm = cellfun(@(x)mean(sqrt(x(:,1).^2 + x(:,2).^2)),buffer_m);
for i = 1:5
    out.(sprintf('eucdm%u',i)) = eucdm(i);
end
out.std_eucdm = std(eucdm);
out.mean_eucdm = mean(eucdm);

% Standard deviation of Euclidean distances in each segment
eucds = cellfun(@(x)std(sqrt(x(:,1).^2 + x(:,2).^2)),buffer_m);
for i = 1:5
    out.(sprintf('eucds%u',i)) = eucds(i);
end
out.std_eucds = std(eucds);
out.mean_eucds = mean(eucds);

% Maximum volume in each segment
% (defined as area of rectangle of max span in each direction)
maxspanx = cellfun(@(x)range(x(:,1)),buffer_m);
maxspany = cellfun(@(x)range(x(:,2)),buffer_m);
spanareas = maxspanx.*maxspany;
out.stdspana = std(spanareas);
out.meanspana = mean(spanareas);

% ------------------------------------------------------------------------------
% Outliers in the embedding space
% ------------------------------------------------------------------------------
% area of max span of all points; versus area of max span of 50% of points closest to origin
d = sqrt(m(:,1).^2 + m(:,2).^2);
[~, ix] = sort(d,'ascend');

out.areas_all = range(m(:,1))*range(m(:,2));
r50 = ix(1:round(end/2)); % 50% of point closest to origin
out.areas_50 = range(m(r50,1))*range(m(r50,2));
out.arearat = out.areas_50 / out.areas_all;

end
