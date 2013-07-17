function y_embed = BF_embed(y,tau,m,sig)
% Embeds the time series y into an m dimensional space at a time delay tau
% Uses the TSTOOL code 'embed'
% INPUTS:
% y: univariate scalar time series
% tau: time-delay. Can be a string, 'ac', 'mi', ...
% m: embedding dimension. Must be a cell specifying method and parameters,
% e.g., {'fnn',0.1} does fnn method using a threshold of 0.1...
% sig [opt]: if 1, uses TSTOOL to embed and returns a signal object.
%           (default = 0, i.e., not to do this and instead return matrix).
%           If 2, returns a vector of [tau m] rather than any explicit embedding
% OUTPUT:
% a matrix of width m containing the vectors in the new embedding space...
% Ben Fulcher October 2009

bevocal = 0; % display information about embedding
N = length(y); % length of the input time series, y

%% (1) Time-delay, tau
if nargin < 2 || isempty(tau)
    tau = 1; % default time delay is 1
    sstau = 'to default of 1';
else
    if ischar(tau) % use a routine to inform tau
        switch tau
            case 'mi' % first minimum of mutual information function
                tau = CO_firstmin(y,'mi');
                sstau = sprintf('by first minimum of mutual information to tau = %u');
            case 'ac' % first zero-crossing of ACF
                tau = CO_fzcac(y);
                sstau = sprintf('by first zero crossing of autocorrelation function to tau = %u',tau);
            otherwise
                error('Invalid time-delay method ''%s''.',tau)
        end
    else
        sstau = sprintf('by user to %u',tau);
    end
end
% we now have an integer time delay tau
% explanation stored in string sstau for printing later

%% (2) Set embedding dimension, m
if nargin < 3 || isempty(m) % set to default value
    m = 2; % embed in 2-dimensional space by default! Probably not a great default!
    ssm = sprintf('to (strange) default of %u',m);
else % use a routine to inform m
    if ~iscell(m), m = {m}; end
    if ischar(m{1})
        switch m{1}
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
                if length(m) == 1 % set default threshold 0.1
                    th = 0.1;
                else
                    th = m{2};
                end
                try
                    m = NL_crptool_fnn(y,10,2,tau,th);
                catch
                    fprintf(1,'Error with FNN code')
                    y_embed = NaN;
                    return
                end
                ssm = sprintf('by N. Marwan''s CRPtoolbox ''fnn'' code with threshold %f to m = %u',th,m);
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
                    m = TSTL_cao(y,10,tau,3,0.2,{'mmthresh',th});
                catch
                    fprintf(1,'Call to TSTOOL function ''cao'' failed')
                    y_embed = NaN; return
                end
                ssm = sprintf('by TSTOOL function ''cao'' using ''mmthresh'' with threshold %f to m = %u',th,m);
            otherwise
                fprintf(1,'embedding dimension, m, incorrectly specified. Exiting.')
                y_embed = NaN; return
        end
    else
        m = m{1};
        ssm = sprintf('by user to %u',m);
    end
end
% we now have an integral embedding dimension, m

%% Now do the embedding
if nargin < 4
    sig = 0; % don't return a signal object, return a matrix
end

if sig == 2 % just return the embedding parameters
    y_embed = [tau, m];
    return
end

% Use the TSTOOL embed function.
if size(y,2) > size(y,1) % make sure it's a column vector
    y = y';
end
try
    y_embed = embed(signal(y),m,tau);
catch me
    if strcmp(me.message,'time series to short for chosen embedding parameters')
        fprintf(1,'Time series (N = %u) too short to embed\n',N);
        y_embed = NaN; return
    else
        % could always try optimizing my own routine (below) so TSTOOL is not required for this step...
        error('Embedding time series using TSTOOL function ''embed'' failed')
    end
end

if ~sig
   y_embed = data(y_embed);
   % this is actually faster than my implementation, which is commented out below
end

if bevocal
    fprintf(1,'Time series embedded using time delay, tau = %s and embedding dimension m = %s',sstau,ssm);
end


% if sig
%     y_embed = embed(signal(y),m,tau);
% else
% %     % my own routine
% %     % create matrix of element indicies
% %     % m wide (each embedding vector)
% %     % N-m*tau long (number of embedding vectors)
% %     y_embed = zeros(N-(m-1)*tau,m);
% %     for i=1:m
% %        y_embed(:,i) = (1+(i-1)*tau:1+(i-1)*tau+N-(m-1)*tau-1)';
% %     end
% % %     keyboard
% %     y_embed = arrayfun(@(x)y(x),y_embed); % take these elements of y
%     
%     % Shit, it's actually faster for large time series to use the TSTOOL version!!:
%     y_embed = data(embed(signal(y),m,tau));
%     
% end


end