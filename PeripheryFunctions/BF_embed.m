function y_embed = BF_embed(y,tau,m,makeSignal,randomSeed,beVocal)
% BF_embed  Time-delay embedding
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
% makeSignal [opt], if 1, uses TSTOOL to embed and returns a signal object.
%           (default = 0, i.e., not to do this and instead return matrix).
%           If 2, returns a vector of [tau m] rather than any explicit embedding
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%
%---OUTPUT:
% A matrix of width m containing the vectors in the new embedding space...
%
% The makeSignal option uses the TSTOOL code 'embed'
% TSTOOL: http://www.physik3.gwdg.de/tstool/

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

N = length(y); % length of the input time series, y

% randomSeed: how to treat the randomization
if nargin < 5
    randomSeed = []; % default
end
if nargin < 6
    beVocal = 0; % by default, do not display information about the embedding
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
                sstau = sprintf('by first minimum of mutual information to tau = ');
            case 'ac' % first zero-crossing of ACF
                tau = CO_FirstZero(y,'ac');
                sstau = sprintf('by first zero crossing of autocorrelation function to tau = ');
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
            case 'tisean'
                % Ben Fulcher, 2015-03-21
                % Uses TISEAN false_nearest code
                if length(m) == 1
                    th = 0.4;
                else
                    th = m{2};
                end
                m = NL_TISEAN_fnn(y,tau,10,0.05,1,th);
                ssm = sprintf('by TISEAN false_nearest code with 5% theiler window and threshold %f to m = %u',th,m);

            case 'fnnsmall'
                % uses Michael Small's fnn code
                if length(m) == 1
                    th = 0.01;
                else
                    th = m{2};
                end
                m = MS_unfolding(y,th,1:10,tau);
                ssm = sprintf('by Michael Small''s FNN code with threshold %f to m = %u',th,m);

            case 'fnnmar'
                % uses Marwin's fnn code in CRPToolbox
                % should specify threshold for proportion of fnn
                % default is 0.1
                % e.g., {'fnnmar',0.2} would use a threshold of 0.2
                % uses time delay determined above

                if length(m) == 1 % no threshold specified
                    th = 0.4; % set default threshold 0.4
                else
                    th = m{2};
                end
                try
                    m = NL_crptool_fnn(y,10,2,tau,th,randomSeed);
                catch
                    fprintf(1,'Error with FNN code');
                    y_embed = NaN;
                    return
                end
                ssm = sprintf('by N. Marwan''s CRPtoolbox (GPL) ''fnn'' code with threshold %f to m = %u',th,m);

            case 'cao'
                % Uses TSTOOL code for cao method to determine optimal
                % embedding dimension
                % max embedding dimension of 10
                % time delay determined by above method
                % 3 nearest neighbours
                % 20% of time series length as reference points
                if length(m) == 1
                    th = 10;
                end
                try
                    m = NL_CaosMethod(y,10,tau,3,0.2,{'mmthresh',th});
                catch
                    fprintf(1,'Call to TSTOOL function ''cao'' failed');
                    y_embed = NaN; return
                end
                ssm = sprintf('by TSTOOL function ''cao'' using ''mmthresh'' with threshold %f to m = %u',th,m);

            otherwise
                error('Embedding dimension, m, incorrectly specified.')
        end
    else
        m = m{1};
        ssm = sprintf('by user to %u',m);
    end
end
% we now have an integral embedding dimension, m

% ------------------------------------------------------------------------------
%% Do the embedding
% ------------------------------------------------------------------------------
if nargin < 4
    makeSignal = 0; % Don't return a signal object, return a matrix
end

if makeSignal == 2 % Just return the embedding parameters
    y_embed = [tau, m];
    return
end

% Make sure it's a column vector:
if size(y,2) > size(y,1)
    y = y';
end

if makeSignal
    % Use the TSTOOL embed function:
    try
        y_embed = embed(signal(y),m,tau);
    catch me
        if strcmp(me.message,'Time series to short for chosen embedding parameters')
            fprintf(1,'Time series (N = %u) too short to embed\n',N);
            y_embed = NaN; return
        else
            % Could always try optimizing my own routine (below) so TSTOOL is not required for this step...
            error('Embedding time series using TSTOOL function ''embed'' failed: %s',me.message)
        end
    end

    % if ~makeSignal
    %    y_embed = data(y_embed);
    %    % this is actually faster than my implementation, which is commented out below
    % end

else
    % Use a Matlab-based implementation:
    N_embed = N-(m-1)*tau;
    if N_embed <=0
        error('Time Series (N = %u) too short to embed with these embedding parameters',N);
    end

    % Each embedding vector is a row (of length m columns)
    % Number of ebmedding vectors is N_embed = N - (m-1)*tau
    y_embed = zeros(N_embed,m);

    for i = 1:m
       y_embed(:,i) = y(1+(i-1)*tau:N_embed+(i-1)*tau);
    end
end

% Tell me about it:
if beVocal
    fprintf(1,['Time series embedded successfully:\n--Time delay %s%u\n' ...
                        '--Embedding dimension m %s\n'],sstau,tau,ssm);
end


end
