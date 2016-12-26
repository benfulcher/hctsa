function out = NL_MS_nlpe(y,de,tau,maxN)
% NL_MS_nlpe   Normalized drop-one-out constant interpolation nonlinear prediction error.
%
% Computes the nlpe for a time-delay embedded time series using Michael Small's
% code, nlpe (renamed MS_nlpe here):
%
%---INPUTS:
% y, the input time series
% de, the embedding dimension (can be an integer, or 'fnn' to select as the
%       point where the proportion of false nearest neighbors falls below 5%
%       using NL_MS_fnn)
% tau, the time-delay (can be an integer or 'ac' to be the first zero-crossing
%       of the ACF or 'mi' to be the first minimum of the automutual information
%       function)
%
%---OUTPUTS: include measures of the meanerror of the nonlinear predictor, and a
% set of measures on the correlation, Gaussianity, etc. of the residuals.
%
% cf. M. Small, Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance (book) World Scientific, Nonlinear Science Series A,
% Vol. 52 (2005)
%
% Michael Small's Matlab code is available at http://small.eie.polyu.edu.hk/matlab/

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

if ~BF_iszscored(y)
    warning('The input time series should be z-scored')
end

N = length(y); % time-series length

% ------------------------------------------------------------------------------
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------
% Embedding dimension:
if nargin < 2 || isempty(de)
    de = 3;
end

% Time-delay, tau
if nargin < 3 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac');
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi');
end

if nargin < 4 || isempty(maxN)
    % Set a maximum time-series length to analyze, due to memory constraints
    % with longer time series
    maxN = 5000;
end

% nlpe can cause memory pains for long time series
% Let's do this dirty cheat
if N > maxN
    % Crop the time series to the first maxN samples:
    y = y(1:maxN);
    fprintf(1,['Michael Small''s ''nlpe'' code is only being evaluated ' ...
                    'on the first %u (/%u) time-series samples...\n'],maxN,N);
    N = maxN;
end
if N < 20 % Short time series cause problems:
    fprintf(1,'Time series (N = %u) is too short\n',length(y));
    out = NaN; return
end

% Do false nearest neighbours to compute an appropriate embedding dimension, if needed
if strcmp(de,'fnn')
    de = NL_MS_fnn(y,1:10,tau,5,1,1,0.05);
end

% Run Michael Small's nonlinear prediction error code:
res = MS_nlpe(y,de,tau); % residuals

% ------------------------------------------------------------------------------
%% Compute outputs
% ------------------------------------------------------------------------------
out.msqerr = mean(res.^2);
% out.rmserr = sqrt(mean(e.^2));
% out.mabserr = mean(abs(e));
% out.meanres = mean(e);

% Use MF_ResidualAnalysis on the residuals
% 1) Get statistics on residuals
residstats = MF_ResidualAnalysis(res);

% Convert output of residual analysis to local outputs in quick loop
fields = fieldnames(residstats);
for k = 1:length(fields);
    out.(fields{k}) = residstats.(fields{k});
end

end
