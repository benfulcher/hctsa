function out = EX_MovingThreshold(y,a,b)
% EX_MovingThreshold    Moving threshold model for extreme events in a
%                       time series
%
% Inspired by an idea contained in:
% "Reactions to extreme events: Moving threshold model"
% Altmann et al., Physica A 364, 435--444 (2006)
%
% This algorithm is based on this idea: it uses the occurrence of extreme events
% to modify a hypothetical 'barrier' that classes new points as 'extreme' or not.
% The barrier begins at sigma, and if the absolute value of the next data point
% is greater than the barrier, the barrier is increased by a proportion 'a',
% otherwise the position of the barrier is decreased by a proportion 'b'.
%
%---INPUTS:
% y, the input (z-scored) time series
% a, the barrier jump parameter (in extreme event)
% b, the barrier decay proportion (in absence of extreme event)
%
%---OUTPUTS: the mean, spread, maximum, and minimum of the time series for the
% barrier, the mean of the difference between the barrier and the time series
% values, and statistics on the occurrence of 'kicks' (times at which the
% threshold is modified), and by how much the threshold changes on average.
%
% In future could make a variant operation that optimizes a and b to minimize the
% quantity meanqover/pkick (hugged the shape as close as possible with the
% minimum number of kicks), and returns a and b...?

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

doPlot = 0; % can set to 1 to plot outputs

% ------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
% Check that the time series is z-scored (just a warning)
if ~BF_iszscored(y)
    warning('The input time series should be z-scored')
end

if nargin < 2 || isempty(a)
    a = 1; % set default
end

if nargin < 2 || isempty(b)
    b = 0.1; % set default
end
if (b < 0) || (b > 1)
    error('The decay proportion, b, should be between 0 and 1');
end

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % time-series length
y = abs(y); % extreme events defined in terms of absolute deviation from mean
q = zeros(N,1); % the barrier
kicks = zeros(N,1);

% Treat the barrier as knowing nothing about the time series, until it encounters it
% (except for the std! -- starts at 1)

% Initial condition of barrier q:
% The barrier will get smarter about the distribution
% but will decay to simulate 'forgetfulness' in the original model(!)

q(1) = 1; % begin at sigma

for i = 2:N
	if y(i) > q(i-1) % Extreme event -- time series value more extreme than the barrier
		% q(i) = (1+a)*q(i-1); % increase barrier by proportion a
		q(i) = (1+a)*y(i); % increase barrier above the new observation by a factor a
		kicks(i) = q(i) - q(i-1); % The size of the increase
	else
		q(i) = (1-b)*q(i-1); % Decrease barrier by proportion b
	end
end

if doPlot
    figure('color','w'); box('on')
    hold on
    plot(y,'.-k')
    plot(q,'--r')
    hold off
end

% ------------------------------------------------------------------------------
%% Outputs
% ------------------------------------------------------------------------------

% Basic statistics on the barrier dynamics, q
out.meanq = mean(q);
out.medianq = median(q);
out.iqrq = iqr(q);
out.maxq = max(q);
out.minq = min(q);
out.stdq = std(q);
out.meanqover = mean(q - y);

% Kicks (when the barrier is changed due to extreme event)
out.pkick = sum(kicks)/(N-1); % probability of a kick
fkicks = find(kicks);
I_kick = diff(fkicks); % time intervals between successive kicks
out.stdkickf = std(I_kick);
out.meankickf = mean(I_kick);
out.mediankickf = median(I_kick);

end
