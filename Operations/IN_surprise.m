function out = IN_surprise(y,whatinf,memory,ng,cgmeth,nits)
% How surprised you are at this data point given the information recorded
% by whatinf (the information recorded) at a given memory in the past
% ng is number of groups to coarse-grain into
% cgmeth is the method for coarse-graining: 'quantile' or 'updown'
% Ben Fulcher, September 2009

%% Check inputs and set defaults
if nargin < 2 || isempty(whatinf)
    whatinf = 'dist'; % expect probabilities based on prior observed distribution
end
% Memory: how far into the past to base your priors on
if nargin < 3 || isempty(memory)
    memory = 0.2; % set it as 20% of the time-series length
end
if memory > 0 && memory < 1 % specify memory as a proportion of the time-series length
    memory = round(memory*length(y));
end
% ng -- number of groups for the time-series coarse-graining/symbolization
if nargin < 4 || isempty(ng)
    ng = 3; % use three symbols to approximate the time-series values
end
% cgmeth: the coase-graining method to use
if nargin < 5 || isempty(cgmeth)
    cgmeth = 'quantile'; % symbolize time series by their values (quantile)
end
% nits: number of iterations
if nargin < 6 || isempty(nits)
    nits = 500;
    % number of iterations of the procedure to perform (does it with random samples)
    % could also imagine doing it exhaustively...?!
end

%% Course Grain
yth = SUB_coursegrain(y,ng,cgmeth); % a coarse-grained time series using the numbers 1:ng

N = length(yth); % will be the same as y, for 'quantile', and 'updown'

%% get prior information
rs = randperm(N-memory) + memory;
rs = rs(1:min(nits,end));

store = zeros(nits,1);
for i = 1:length(rs)
    switch whatinf
        case 'dist'
            % Uses the distribution up to memory to inform the next point
            
            % Calculate probability of this given past memory
            p = sum(yth(rs(i)-memory:rs(i)-1) == yth(rs(i)))/memory;
            store(i) = p;
            
        case 'T1'
            % Uses one-point correlations in memory to inform the next point
            
            % Estimate transition probabilities from data in memory
            % Find where in memory this has been observed before, and what
            % preceeded it:
            memorydata = yth(rs(i)-memory:rs(i)-1);
            % Previous value observed in memory here:
            inmem = find(memorydata(1:end-1) == yth(rs(i)-1));
            if isempty(inmem)
                p = 0;
            else
                p = sum(memorydata(inmem+1) == yth(rs(i)))/length(inmem);
            end
            store(i) = p;
            
        case 'T2'
            % Uses two-point correlations in memory to inform the next point
            
            memorydata = yth(rs(i)-memory:rs(i)-1);
            % Previous value observed in memory here:
            inmem1 = find(memorydata(2:end-1) == yth(rs(i)-1)); % the 2:end makes the next line ok...?
            inmem2 = find(memorydata(inmem1) == yth(rs(i)-2));
            if isempty(inmem2)
                p = 0;
            else
                p = length(find(memorydata(inmem2+2) == yth(rs(i))))/length(inmem2);
            end
            store(i) = p;
    end
end

% information gained from next observation is log(1/p) = -log(p)
iz = find(store == 0);
store(iz) = 1; % to avoid log(0) error in next line
store = -log(store); % transform to surprises/information gains
store(iz) = 0; % so that log(0) = 0
% May be strange? maybe remove these points rather than setting to zero?
% plot(store)

out.min = min(store); % minimum amount of information you can gain in this way
out.max = max(store); % maximum amount of information you can gain in this way
out.mean = mean(store); % mean
out.median = median(store); % median
out.lq = quantile(store,0.25); % lower quartile
out.uq = quantile(store,0.75); % upper quartile
out.std = std(store); % standard deviation

% t-statistic to information gain of 1
out.tstat = abs((out.mean-1)/(out.std/sqrt(nits)));


end