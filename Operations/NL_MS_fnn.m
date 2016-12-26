function out = NL_MS_fnn(y,de,tau,th,kth,justBest,bestp)
% NL_MS_fnn     False nearest neighbors of a time series.
%
% Determines the number of false nearest neighbors for the embedded time series
% using Michael Small's false nearest neighbor code, fnn (renamed MS_fnn here).
%
% False nearest neighbors are judged using a ratio of the distances between the
% next k points and the neighboring points of a given datapoint.
%
%---INPUTS:
% y, the input time series
%
% de, the embedding dimensions to compare across (a vector)
%
% tau, the time-delay (can be 'ac' or 'mi' to be the first zero-crossing of ACF,
%                       or first minimum of AMI, respectively)
%
% th, the distance threshold for neighbours
%
% kth, the the distance to next points
%
% [opt] justBest, can be set to 1 to just return the best embedding dimension, m_{best}
%
% [opt] bestp, if justBest = 1, can set bestp as the proportion of false nearest
%              neighbours at which the optimal embedding dimension is selected.
%
% This function returns statistics on the proportion of false nearest neighbors
% as a function of the embedding dimension m = m_{min}, m_{min}+1, ..., m_{max}
% for a given time lag tau, and distance threshold for neighbors, d_{th}.
%
%---OUTPUTS: include the proportion of false nearest neighbors at each m, the mean
% and spread, and the smallest m at which the proportion of false nearest
% neighbors drops below each of a set of fixed thresholds.

% cf. M. Small, Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance (book) World Scientific, Nonlinear Science Series A,
% Vol. 52 (2005)
% Code available at http://small.eie.polyu.edu.hk/matlab/
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
%% CHECK INPUTS:
% ------------------------------------------------------------------------------
% Embedding dimension(s), de
if nargin < 2 || isempty(de)
    de = (1:10);
end

% Time delay, tau
if nargin < 3 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac'); % first zero-crossing of autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi'); % first minimum of automutual information function
end

% A distance threshold for neighbours
if nargin < 4
    th = 5;
end

% Distance to next points
if nargin < 5
    kth = 1;
end

% (Actually better to use MS_unfolding now -- does a near-identical thing
% to this...)
if nargin < 6 || isempty(justBest)
    justBest = 0;
end

if nargin < 7 || isempty(bestp)
    bestp = 0.1; % first time under 10% of neighest neighbours
end

% ------------------------------------------------------------------------------
%% Run Michael Small's false nearest neighbour code, MS_fnn:
% ------------------------------------------------------------------------------
p = MS_fnn(y,de,tau,th,kth);

% ------------------------------------------------------------------------------
%% Now make output
% ------------------------------------------------------------------------------
% Assuming we've set tau, and m is a vector, we should have p (the
% proportion of false neighbours) and de (the corresonding embedding
% dimensions) as vectors

if justBest
    % We just want a scalar to choose the embedding with
    out = firstunderf(bestp);
    return
else
    % Output all of them
    for i = 1:length(de)
        out.(sprintf('pfnn_%u',de(i))) = p(i);
    end

    % Output mean
    out.meanpfnn = mean(p);

    % Standard deviation
    out.stdpfnn = std(p);

    % Find embedding dimension for the first time p goes under x%
    out.firstunder02 = firstunderf(0.2); % 20%
    out.firstunder01 = firstunderf(0.1); % 10%
    out.firstunder005 = firstunderf(0.05); % 5%
    out.firstunder002 = firstunderf(0.02); % 2%
    out.firstunder001 = firstunderf(0.01); % 1%

    % Maximum step-wise change across p
    out.max1stepchange = max(abs(diff(p)));
end

% ------------------------------------------------------------------------------
function firsti = firstunderf(x)
    %% Find de for the first time p goes under x%
    firsti = de(find(p < x,1,'first'));
    if isempty(firsti)
        firsti = de(end) + 1;
    end
end

end
