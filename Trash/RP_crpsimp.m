% RP_crpsimp
% 
% REDUNDANT
% 

function out = RP_crpsimp(y,m,tau,eps,dist,tw)
% My attempt to do a bunch of recurrence plot measures using the CRP
% Toolbox by Marwan, Version 5.13, Release 26
% This function isn't actually used -- quite buggy software
% y should be z-scored
% Ben Fulcher, October 2009

%% How to deal with long time series
% Now, because this only really works well for short time series, which is a bit
% of a pity, we have a few options:
% 1) restrict use only to short time series (l<2000)
% 2) take the first 2000 samples for longer time series
% 3) downsample longer time series and compute on that
% Only taking the first 2000 samples I think is best, the downsampling
% could produce some spurious recurrences... I'll do that for now...

N = length(y); % time-series length

maxN = 2000; % really time consuming if input time series is longer than this threshold...
if N > maxN
    y = y(1:maxN);
    N = maxN;
    fprintf(1,'Analyzing just the first %u points of the input time series\n', maxN)
end

%% 1) Set parameters
% Set time-delay, tau
if nargin < 3 || isempty(tau)
    tau = 'mi'; %use first minimum of the automutual information function
end
if ischar(tau)
    switch tau
    case 'ac'
        tau = CO_FirstZero(y,'ac');
    case 'mi'
        tau = CO_FirstMin(y,'mi');
    otherwise
        error('Unknown tau setting: ''%s''',tau);
    end
end
if tau > N/10
    tau = max(1,floor(N/10));
end

% Set embedding dimension, m -- maybe use false nearest neighbours
if nargin < 2 || isempty(m)
    m = 3; % use an embedding dimension of 3 by default
end
if strcmp(m,'fnn')
    % first time fnn goes below 0.1 (up to maximum of 10)
    m = NL_crptool_fnn(y,10,2,tau,0.1);
end

% Set the neighbourhood size, eps
if nargin < 4 || isempty(eps)
    eps = 1;
end


% Set the distance metric, dist
if nargin < 5 || isempty(dist)
    dist = 'euclidean';
else
    dists = {'maxnorm','euclidean','minnorm','nrmnorm','rr','fan',...
        'inter','omatrix','opattern','distance'};
    if ~ismember(dists,dist)
        error('Invalid distance measure')
    end
end

% Set the Theiler Window, tw
if nargin < 6 || isempty(tw)
    tw = 1;
end

% Set the minimum length of diagonal/vertical structure, lmin and vmin
lmin = 2;
vmin = 2;

%% Compute CRQA measures

rqa = crqa(y,m,tau,eps,[],[],lmin,vmin,tw,dist,'nonormalize','silent');

out.recur_rate = rqa(1); % recurrence rate
out.det = rqa(2); % determinism measure
out.meandiagl = rqa(3); % mean diagonal length
out.longdiag = rqa(4); % length of longest diagonal line
out.entdiagl = rqa(5); % entropy of diagonal line lengths
out.laminarity = rqa(6); % laminarity
out.trapt = rqa(7); % trapping time
out.longvert = rqa(8); % length of longest vertical line
out.rec1 = rqa(9); % recurrence time of the 1st type
out.rec2 = rqa(10); % recurrence time of the 2nd type

%% Compute CRP
% In princpile could extract measures myself off this, but also quite
% computationally expensive; will just stick to the statistics above from
% the Marwan implementation. Could also use built in trapping time (tt) and
% diagonal line lengths (dl) routines...
% Would be good to see how statistics on this varied with the parameters;
% particularly epsilon, for an informed choice of m and tau.

%% Compute Diagonalwise CRQA measures

% rqad = crqad(y,m,tau,eps,[],'nonormalize','silent');

% out.recur_rate_xy = rqad.RRp;
% out.recur_rate_xny = 


% Not sure how these scale with sample size...


% %% Compute CRP
% % Using CRP Toolbox code
% 
% rp = crp(y,m,tau,eps,dist,'nonormalize','silent');
% spy(rp);
% % keyboard
% 
% rpl = size(rp,1);
% 
% %% Get some measures out
% 
% % percentage of points that are recurrence points
% out.precurr = length(find(rp))/rpl^2;
% 
% % entropy of distribution in triangles in the two-dimensional space
% % 5 triangles in upper diagonal with horizontal bases
% % Can't be bothered, actually...!
% % x = rpl*(1-sqrt((1:4)/5));
% % bins = zeros(5,1);
% % bins(1) = length(find());



end