function y_embed = BF_Embed(y,tau,m,justGiveMeParams,randomSeed,beVocal)
% BF_Embed  Time-delay embedding
%
% Returns a time-delay embedding of the input time series into an m dimensional
% space at a time delay tau.
%
%---INPUTS:
% y, univariate scalar time series
%
% tau, time-delay. Can be a string, 'ac', 'mi', ...
%
% m, the embedding dimension. Must be a cell specifying method and parameters,
%    e.g., {'fnn',0.1} does fnn method using a threshold of 0.1...
%
% justGiveMeParams [opt], if true returns a vector of [tau,m] rather than doing any
%           actual embedding (default = 0, i.e., do the embedding).
%
% randomSeed [currently unused], the sole remaining embedding-dimension
%           method ('tisean') is deterministic given its inputs, so this
%           parameter has no effect; kept only so existing positional calls
%           don't need updating.
%
%---OUTPUT:
% A matrix of width m containing the vectors in the new embedding space...

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

N = length(y); % length of the input time series, y

% ---[[NOTE: other input checking done later, below]]---
% randomSeed: how to treat the randomization
if nargin < 5
    randomSeed = []; % default
end
if nargin < 6
    beVocal = false; % by default, do not display information about the embedding
end

% ------------------------------------------------------------------------------
%% (1) Time-delay, tau
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 1; % default time delay is 1
    sstau = 'to default of 1';
else
    if ischar(tau) % use a routine to inform tau
        switch tau
            case 'mi' % first minimum of mutual information function
                tau = CO_FirstMin(y,'mi');
                if isnan(tau)
                    % Could not get time delay by mutual information (time series too short?)
                    y_embed = NaN; return
                end
                sstau = sprintf('by first minimum of mutual information to tau = ');
            case 'ac' % first zero-crossing of ACF
                tau = CO_FirstCrossing(y,'ac',0,'discrete');
                sstau = sprintf('by first zero crossing of autocorrelation function to tau = ');
                if isnan(tau)
                    % Could not get time delay by ACF (time series too short?)
                    y_embed = NaN; return
                end
            otherwise
                error('Invalid time-delay method ''%s''.',tau)
        end
    else
        sstau = sprintf('by user to %u',tau);
    end
end
% we now have an integer time delay tau
% explanation stored in string sstau for printing later

% ------------------------------------------------------------------------------
%% Determine the embedding dimension, m
% ------------------------------------------------------------------------------
if nargin < 3 || isempty(m) % set to default value
    m = 2; % Embed in 2-dimensional space by default
    ssm = sprintf('to (strange) default of %u',m);
else % use a routine to inform m
    if ~iscell(m), m = {m}; end
    if ischar(m{1})
        switch m{1}
            case 'fnn'
                % Selects an embedding dimension via false nearest neighbors,
                % using TISEAN's false_nearest code. escapeFactor=5 matches
                % the false-neighbor escape-ratio threshold conventionally
                % used for this kind of test (TISEAN's own internal default
                % of 2.0 is considerably stricter and was found, in practice,
                % to make the test look like it never converges even for
                % textbook low-dimensional systems).
                if length(m) == 1
                    th = 0.4;
                else
                    th = m{2};
                end
                escapeFactor = 5;
                m = NL_FNN(y,tau,10,0.05,1,th,escapeFactor);
                ssm = sprintf('by TISEAN false_nearest code with 5%% theiler window and threshold %f to m = %u',th,m);

            case 'fnnsmall'
                % Uses Michael Small's fnn code. Not used by any operation in
                % the current hctsa feature set (standardized on 'fnn'
                % above), but kept available as an alternative.
                if length(m) == 1
                    th = 0.01;
                else
                    th = m{2};
                end
                m = MS_unfolding(y,th,1:10,tau);
                ssm = sprintf('by Michael Small''s FNN code with threshold %f to m = %u',th,m);

            otherwise
                error('Embedding dimension, m, incorrectly specified.')
        end
    else
        m = m{1};
        ssm = sprintf('by user to %u',m);
    end
end
if isnan(m)
    % Could not determine an embedding dimension (time series too short/degenerate?)
    y_embed = NaN; return
end
% we now have an integral embedding dimension, m

% ------------------------------------------------------------------------------
%% Do the embedding
% ------------------------------------------------------------------------------
if nargin < 4
    justGiveMeParams = false; % Don't return a signal object, return a matrix
end

if justGiveMeParams % Just return the embedding parameters
    y_embed = [tau, m];
    return
end

% Make sure it's a column vector:
if size(y,2) > size(y,1)
    y = y';
end

% Matlab-based matrix embedding:
N_embed = N-(m-1)*tau;
if N_embed <=0
    % Time series too short to embed with these embedding parameters
    y_embed = NaN; return
end

% Each embedding vector is a row (of length m columns)
% Number of embedding vectors is N_embed = N - (m-1)*tau
y_embed = zeros(N_embed,m);

for i = 1:m
   y_embed(:,i) = y(1+(i-1)*tau:N_embed+(i-1)*tau);
end

% Tell me about it:
if beVocal
    fprintf(1,['Time series embedded successfully:\n--Time delay %s%u\n' ...
                        '--Embedding dimension m %s\n'],sstau,tau,ssm);
end


end
