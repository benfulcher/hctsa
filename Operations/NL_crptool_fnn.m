function out = NL_crptool_fnn(y,maxm,r,taum,th,randomSeed)
% NL_crptool_fnn    Analyzes the false-nearest neighbours statistic.
%
%---INPUTS:
% y, the input time series
% maxm, the maximum embedding dimension to consider
% r, the threshold; neighbourhood criterion
% taum, the method of determining the time delay, 'corr' for first zero-crossing
%       of autocorrelation function, or 'mi' for the first minimum of the mutual
%       information
%
% th [opt], returns the first time the number of false nearest neighbours drops
%           under this threshold
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed

% Computation uses N. Marwan's code from the CRP Toolbox:
% http://tocsy.pik-potsdam.de/CRPtoolbox/
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
%% Preliminaries, input checking
% ------------------------------------------------------------------------------
doPlot = 0; % plot outputs to figure
N = length(y); % length of the input time series

% 1) maxm: the maximum embedding dimension
if nargin < 2 || isempty(maxm)
    maxm = 10; % default maximum embedding dimension
end

% 2) r, the neighbourhood criterion
if nargin < 3 || isempty(r)
    r = 2; % neighbourhood criterion
end

% 3) determine the time delay
if nargin < 4 || isempty(taum)
    taum = 'mi'; % by default determine time delay by first minimum of AMI
end
if ischar(taum)
    if strcmp(taum,'mi')
        tau = CO_FirstMin(y,'mi'); % time-delay
    elseif strcmp(taum,'ac')
        tau = CO_FirstZero(y,'ac'); % time-delay
    else
        error('Invalid time-delay method ''%s''.',taum)
    end
else % give a numeric answer
    tau = taum;
end
% Don't want tau to be too large
if tau > N/10;
    tau = floor(N/10);
end

% 4) Just output a scalar embedding dimension rather than statistics on the
% method?
if nargin < 5
    th = []; % default is to return statistics
end

% 5) randomSeed: how to treat the randomization
if nargin < 6
    randomSeed = []; % default
end

% ------------------------------------------------------------------------------
%% Here's where the action happens:
% ------------------------------------------------------------------------------
if ~exist(fullfile('Marwan_crptool','crptool_fnn'),'file')
    error('Error -- the CRP Toolbox functions for calculating nearest neighbours can not be found');
end

% Control the random seed (for reproducibility):
BF_ResetSeed(randomSeed);

% Run Marwan's CRPToolbox false nearest neighbors code:
nn = crptool_fnn(y,maxm,tau,r,'silent');

if isnan(nn);
    error('Error running the function ''fnn'' from Marwan''s CRPToolbox')
end

if doPlot
    figure('color','w')
    plot(1:maxm,nn,'o-k');
end

if isempty(th) % output summary statistics

    % nn drops
    dnn = diff(nn);
    out.mdrop = mean(dnn); % same information as in fnn(maxm)
    out.pdrop = -sum(sign(dnn))/(maxm-1);

    % fnn
    for i = 2:maxm
        out.(sprintf('fnn%u',i)) = nn(i);
    end

    % first time NN error goes below a set of thresholds
    % firstunderfn = @(x) find(nn < x,1,'first');
    out.firstunder08 = firstunderf(0.8,1:maxm,nn);
    out.firstunder07 = firstunderf(0.7,1:maxm,nn);
    out.firstunder05 = firstunderf(0.5,1:maxm,nn);
    out.firstunder02 = firstunderf(0.2,1:maxm,nn);
    out.firstunder01 = firstunderf(0.1,1:maxm,nn);
    out.firstunder005 = firstunderf(0.05,1:maxm,nn);

else % in this case return a scalar of embedding dimension as output
    out = firstunderf(th,1:maxm,nn);
end


% ------------------------------------------------------------------------------
function firsti = firstunderf(x,m,p)
    %% Find m for the first time p goes under x%
    firsti = m(find(p < x,1,'first'));
    if isempty(firsti)
        firsti = m(end) + 1;
    end
end

end
