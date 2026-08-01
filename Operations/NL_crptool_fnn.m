function out = NL_crptool_fnn(y,maxm,r,taum,th,randomSeed)
% NL_crptool_fnn    Analyzes the false-nearest neighbors statistic.
%
%---INPUTS:
% y, the input time series
% maxm, the maximum embedding dimension to consider
% r, the threshold; neighbourhood criterion
% taum, the method of determining the time delay, 'ac' for first zero-crossing
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
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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
% Resolve default inputs (must happen here, using the caller's nargin, before
% the cache key below is built from these variables -- NL_crptool_fnn_Compute
% below always receives all 6 arguments explicitly, so nargin there would
% otherwise always read as 6).
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(maxm)
    maxm = 10; % default maximum embedding dimension
end
if nargin < 3 || isempty(r)
    r = 2; % neighbourhood criterion
end
if nargin < 4 || isempty(taum)
    taum = 'mi'; % by default determine time delay by first minimum of AMI
end
if nargin < 5
    th = []; % default is to return statistics
end
if nargin < 6
    randomSeed = []; % default
end

% ------------------------------------------------------------------------------
% Cache the resolved output (keyed on an exact, NaN-safe match of all inputs
% that affect the result). This computation resets the random seed to a fixed
% value before running (BF_ResetSeed(randomSeed), which defaults to
% rng(0,'twister')), so it is fully deterministic given its inputs -- a repeat
% call with the same series and parameters always gives the same output,
% making it safe to cache. Several operations share this function's default
% parameters (via BF_Embed's 'fnnmar' embedding-dimension method), each
% independently re-running Marwan's CRPToolbox false-nearest-neighbors code
% for the same series. A miss falls through to exactly the same computation
% as before (see CO_AutoCorr.m for the same pattern and rationale, including
% why this is safe under parfor).
% ------------------------------------------------------------------------------
persistent cacheY cacheMaxm cacheR cacheTaum cacheTh cacheRandomSeed cacheOut
maxCacheEntries = 4;
if isempty(cacheY)
    cacheY = {};
    cacheMaxm = {};
    cacheR = {};
    cacheTaum = {};
    cacheTh = {};
    cacheRandomSeed = {};
    cacheOut = {};
end

for ci = 1:numel(cacheY)
    if isequaln(y,cacheY{ci}) && isequaln(maxm,cacheMaxm{ci}) && isequaln(r,cacheR{ci}) && ...
            isequaln(taum,cacheTaum{ci}) && isequaln(th,cacheTh{ci}) && isequaln(randomSeed,cacheRandomSeed{ci})
        out = cacheOut{ci};
        return
    end
end

out = NL_crptool_fnn_Compute(y,maxm,r,taum,th,randomSeed);

cacheY{end+1} = y;
cacheMaxm{end+1} = maxm;
cacheR{end+1} = r;
cacheTaum{end+1} = taum;
cacheTh{end+1} = th;
cacheRandomSeed{end+1} = randomSeed;
cacheOut{end+1} = out;
if numel(cacheY) > maxCacheEntries
    cacheY(1) = [];
    cacheMaxm(1) = [];
    cacheR(1) = [];
    cacheTaum(1) = [];
    cacheTh(1) = [];
    cacheRandomSeed(1) = [];
    cacheOut(1) = [];
end

end

% ------------------------------------------------------------------------------
function out = NL_crptool_fnn_Compute(y,maxm,r,taum,th,randomSeed)
% The original computation, unchanged -- moved into its own local function so
% the cache wrapper above can store its result. All inputs are already
% resolved (non-default) by the caller above.
% ------------------------------------------------------------------------------
doPlot = 0; % plot outputs to figure
N = length(y); % length of the input time series

% determine the time delay
if ischar(taum)
    % time delay
    if strcmp(taum,'mi')
        tau = CO_FirstMin(y,'mi');
    elseif strcmp(taum,'ac')
        tau = CO_FirstCrossing(y,'ac',0,'discrete');
    else
        error('Invalid time-delay method ''%s''.',taum)
    end
else % give a numeric answer
    tau = taum;
end
% Don't want tau to be too large
if tau > N/10
    tau = floor(N/10);
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

% Output some summary statistics
if isempty(th)

    % nn drops
    dnn = diff(nn);
    out.mdrop = mean(dnn); % same information as in fnn(maxm)
    out.pdrop = -mean(sign(dnn)); % proportion of m -> m+1 for which fnn decreased

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
