function out = CO_PartialAutoCorr(y,maxTau,whatMethod)
% CO_PartialAutoCorr   Compute the partial autocorrelation of an input time series
%
%---INPUTS:
% y, a scalar time series column vector.
%
% maxTau, the maximum time-delay. Returns for lags up to this maximum.
%
% whatMethod, the method used to compute: 'ols' or 'yule_walker'
%
%---OUTPUT: the partial autocorrelations across the set of time lags.
%

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
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
if nargin < 2
    % Use a maximum lag of 10 by default
    maxTau = 10;
end

if nargin < 3 || isempty(whatMethod)
    % ordinary least square by default
    whatMethod = 'ols';
end

%-------------------------------------------------------------------------------
%% Initial checks on maxTau
%-------------------------------------------------------------------------------
N = length(y); % time-series length

assert(maxTau > 0)

if maxTau < 0
    error('Negative time lags not applicable')
end

% ------------------------------------------------------------------------------
%% Do the computation
% ------------------------------------------------------------------------------

pacf = parcorr(y,'NumLags',maxTau,'Method',whatMethod);

% Zero lag is the first entry in the PACF (and should always be 1)

for i = 1:maxTau
    out.(sprintf('pac_%u',i)) = pacf(i+1);
end

end
