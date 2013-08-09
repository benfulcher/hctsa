% NL_crptool_fnn
% 
% Computes and analyzes the false-nearest neighbours statistic.
% 
% Computation is done by referencing N. Marwan's code from the CRP Toolbox:
% http://tocsy.pik-potsdam.de/CRPtoolbox/
% 
% INPUTS:
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

function out = NL_crptool_fnn(y,maxm,r,taum,th)
% Ben Fulcher, October 2009

%% Preliminaries
doplot = 0; % plot outputs to figure
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

% Here's where the action happens:
if ~exist('crptool_fnn')
    error('Error -- the CRP Toolbox functions for calculating nearest neighbours can not be found');
end
nn = crptool_fnn(y,maxm,tau,r,'silent'); % run Marwan's CRPToolbox false nearest neighbors code

if isnan(nn);
    error('Error running the function ''fnn'' from Marwan''s CRPToolbox')
end

if doplot
    figure('color','w')
    plot(1:maxm,nn,'o-k');
end

if isempty(th) % output summary statistics

    % nn drops
    dnn = diff(nn);
    out.mdrop = mean(dnn);
    out.pdrop = -sum(sign(dnn))/(maxm-1);
    
    % fnn
    for i = 2:maxm
        eval(sprintf('out.fnn%u = nn(%u);',i,i));
    end

    % first time NN error goes below a set of thresholds
    firstunderfn = @(x) find(nn < x,1,'first');
    out.m005 = firstunderfn(0.05);
    if isempty(out.m005), out.m005 = maxm + 1; end
    
    out.m01 = firstunderfn(0.1);
    if isempty(out.m01), out.m01 = maxm + 1; end
    
    out.m02 = firstunderfn(0.2);
    if isempty(out.m02), out.m02 = maxm + 1; end
    
    out.m05 = firstunderfn(0.5);
    if isempty(out.m05), out.m05 = maxm + 1; end


else % just want a scalar of embedding dimension as output
    out = find(nn < th,1,'first');
    if isempty(out)
        out = maxm + 1;
    end
end


end