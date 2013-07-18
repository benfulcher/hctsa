function out = NL_MS_fnn(y,de,tau,th,kth,justbest,bestp)
% Wrapper for Michael Small's False Nearest Neighbour Code:
% http://small.eie.polyu.edu.hk/matlab/
% Ben Fulcher 19/2/2010

%% INPUTS

% embedding dimension(s), de
if nargin < 2 || isempty(de)
    de = (1:10);
end

% Time delay, tau
if nargin < 3 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_fzcac(y); % first zero-crossing of autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_firstmin(y,'mi'); % first minimum of automutual information function
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
if nargin < 6 || isempty(justbest)
    justbest = 0;
end
if nargin < 7 || isempty(bestp)
    bestp = 0.05; % first time under 5% of neighest neighbours
end

% Run Michael Small's false nearest neighbour code:
p = MS_fnn(y,de,tau,th,kth);

%% Now make output
% Assuming we've set tau, and m is a vector, we should have p (the
% proportion of false neighbours) and de (the corresonding embedding
% dimensions) as vectors

if justbest
    % We just want a scalar to choose the embedding with
    out = firstunderf(bestp);
    return
else
    % Output all of them
    for i = 1:length(de)
        eval(sprintf('out.pfnn_%u = %f;',de(i),p(i)));
    end

    % Output mean
    out.meanpfnn = mean(p);
    % Standard deviation
    out.stdpfnn = std(p);
        
    %% Find de for the first time p goes under x%
    out.firstunder02 = firstunderf(0.2); % 20%
    out.firstunder01 = firstunderf(0.1); % 10%
    out.firstunder005 = firstunderf(0.05); % 5%
    out.firstunder002 = firstunderf(0.02); % 2%
    out.firstunder001 = firstunderf(0.01); % 1%

    % Maximum step-wise change across p
    out.max1stepchange = max(abs(diff(p)));
end

function firsti = firstunderf(x)
    %% Find de for the first time p goes under x%
    firsti = de(find(p < x,1,'first'));
    if isempty(firsti)
        firsti = de(end) + 1;
    end
end

end