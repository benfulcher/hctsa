function out = IN_surprise(y,annax,memory,ng,cgmeth)
% How surprised you are at this data point given the information recorded
% by annax at a given memory
% ng is number of groups to coarse-grain into
% cgmeth is the method for coarse-graining: 'quantile' or 'updown'
% Ben Fulcher, September 2009

%% Course Grain
yth = SUB_coursegrain(y,ng,cgmeth); % a coarse-grained time series using the numbers 1:ng

N = length(yth); % will be the same as y, for 'quantile', and 'updown'

%% get prior information
nits = 500; % how many iterations of the procedure to perform
rs = randperm(N-memory) + memory;
rs = rs(1:min(nits,end));

store = zeros(nits,1);
for i = 1:length(rs)
    switch annax
        case 'dist'
            % Uses the distribution up to memory to inform the next point
            
            % Calculate probability of this given past memory
            p = sum(yth(rs(i)-memory:rs(i)-1)==yth(rs(i)))/memory;
            store(i) = p;
            
        case 'T1'
            % Uses one-point correlations in memory to inform the next point
            
            % Estimate transition probabilities from data in memory
            % Find where in memory this has been observed before, and what
            % preceeded it:
            memorydata = yth(rs(i)-memory:rs(i)-1);
            % Previous value observed in memory here:
            inmem = find(memorydata(1:end-1)==yth(rs(i)-1));
            if isempty(inmem)
                p = 0;
            else
                p = sum(memorydata(inmem+1)==yth(rs(i)))/length(inmem);
            end
            store(i) = p;
            
        case 'T2'
            % Uses two-point correlations in memory to inform the next point
            
            memorydata = yth(rs(i)-memory:rs(i)-1);
            % Previous value observed in memory here:
            inmem1 = find(memorydata(2:end-1)==yth(rs(i)-1)); % the 2:end makes the next line ok...?
            inmem2 = find(memorydata(inmem1)==yth(rs(i)-2));
            if isempty(inmem2)
                p = 0;
            else
                p = length(find(memorydata(inmem2+2)==yth(rs(i))))/length(inmem2);
            end
            store(i) = p;
    end
end

% information gained from next observation is log(1/p) = -log(p)
iz = find(store==0);
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