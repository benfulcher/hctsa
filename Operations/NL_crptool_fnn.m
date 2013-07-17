function out = NL_crptool_fnn(y,maxm,r,taum,th)
% Uses Marwan's code from the CRP Toolbox
% taum is the method of determining the time delay, 'corr' for first
% zero-crossing of autocorrelation function, or 'mi' for the first minimum
% of the mutual information
% ** If th is set, then just outputs a scalar, the first time the number of
% false neighbours goes under this value
% Ben Fulcher, October 2009

%% Preliminaries
doplot = 1; % plot outputs
N = length(y); % length of the input time series

% 1) maxm: the maximum embedding dimension
if nargin < 2 || isempty(maxm)
    maxm = 10; % default maximum embedding dimension
end

% 2) r, the neubourhood criterion
if nargin < 3 || isempty(r)
    r = 2; % neighbourhood criterion
end

% 3) determine the time delay
if nargin < 4 || isempty(taum)
    taum = 'mi'; % by default determine time delay by first minimum of AMI
end
if ischar(taum)
    if strcmp(taum,'mi')
        tau = CO_firstmin(y,'mi'); % time-delay
    elseif strcmp(taum,'ac')
        tau = CO_fzcac(y); % time-delay
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

% HERE'S WHERE THE ACTION HAPPENS:
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