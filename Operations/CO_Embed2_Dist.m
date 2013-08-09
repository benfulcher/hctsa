% CO_Embed2_Dist
% 
% Returns statistics on the sequence of successive Euclidean distances between
% points in a two-dimensional time-delay embedding space with a given
% time-delay, tau.
% 
% Outputs include the autocorrelation of distances, the mean distance, the
% spread of distances, and statistics from an exponential fit to the
% distribution of distances.
% 
% INPUTS:
% y, a z-scored column vector representing the input time series
% tau, the time delay.
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

function out = CO_Embed2_Dist(y,tau)
% Ben Fulcher, September 2009

doplot = 0; % plot results

%% Check inputs:
if nargin < 2 || isempty(tau)
    tau = 'tau'; % set to the first minimum of autocorrelation function
end

N = length(y); % time-series length

if strcmp(tau,'tau'),
    tau = CO_FirstZero(y,'ac');
    if tau > N/10
        tau = floor(N/10);
    end
end

% Make sure the time series is a column vector
if size(y,2) > size(y,1);
    y = y';
end

m = [y(1:end-tau), y(1+tau:end)];
% m2 = y(1+tau:end);


% Plot some outputs
if doplot
    figure('color','w'); box('on');
    plot(m(:,1),m(:,2),'.');
end

% Calculate Euclidean distances between successive points in this space, d:

d = sum(abs(diff(m(:,1)).^2 + diff(m(:,2)).^2),2); % Correct form, updated 9/7/2011

% Outputs statistics obtained from ordered set of distances between successive points in the recurrence space
out.d_ac1 = CO_AutoCorr(d,1); % Autocorrelation at lag 1
out.d_ac2 = CO_AutoCorr(d,2); % Autocorrelation at lag 2
out.d_ac3 = CO_AutoCorr(d,3); % Autocorrelation at lag 3

out.d_mean = mean(d); % Mean distance
out.d_median = median(d); % Median distance
out.d_std = std(d); % standard deviation of distances
out.d_iqr = iqr(d); % Interquartile range of distances
out.d_max = max(d); % Maximum distance
out.d_min = min(d); % Minimum distance
out.d_cv = mean(d)/std(d); % coefficient of variation of distances

% Empirical distance distribution often fits Exponential distribution quite well
% Fit to all values (often some extreme outliers, but oh well)
% Use a histogram with fixed bins
[n, x] = hist(d,20);
n = n/(sum(n)*(x(2)-x(1))); % normalize to proportional bin counts
l = expfit(d);
nlogL = explike(l,d);
expf = exppdf(x,l);
out.d_expfit_l = l;
out.d_expfit_nlogL = nlogL;
% Sum of abs differences between exp fit and observed:
out.d_expfit_sumdiff = sum(abs(n - expf));


end