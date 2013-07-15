function out = MS_shannon(y,bin,depth)
% Wrapped up from Michael Small's code 'shannon.m', available at
% http://small.eie.polyu.edu.hk/matlab/
% Uses a number 'depth' symbols from a 'bin'-symbol alphabet size

if nargin < 2 || isempty(bin)
    bin = 2; % two bins to discretize the time series, y
end
if nargin < 3 || isempty(depth)
    depth = 3; % three-long strings
end

%% (*)
if length(bin) == 1 && length(depth) == 1
    % Run the code, just return a number
    out = shannon(y,bin,depth)/depth;
end

%% (*) Return statistics over depths (constant number of bins)
% Somewhat strange behaviour -- very variable
if length(bin) == 1 && length(depth) > 1
    % range over depths specified in the vector
    % return statistics on results
    ndepths = length(depth);
    ents = zeros(ndepths,1);
    for i = 1:ndepths
        ents(i) = shannon(y,bin,depth(i));
    end
    % should scale with depth: normalize by this:
    ents = ents./depth';
    out.maxent = max(ents);
    out.minent = min(ents);
    out.medent = median(ents);
    out.meanent = mean(ents);
    out.stdent = std(ents);
end

%% (*) Statistics over different bin numbers
if length(bin) > 1 && length(depth) == 1
    % range over depths specified in the vector
    % return statistics on results
    nbins = length(bin);
    ents = zeros(nbins,1);
    for i = 1:nbins
        ents(i) = shannon(y,bin(i),depth);
    end
    ents = ents/depth; % should scale with depth: normalize by this:
    out.maxent = max(ents);
    out.minent = min(ents);
    out.medent = median(ents);
    out.meanent = mean(ents);
    out.stdent = std(ents);
end

%% (*) statistics over both bins and depths
if length(bin) > 1 && length(depth) > 1
    nbins = length(bin);
    ndepths = length(depth);
    
    ents = zeros(nbins,ndepths);
    for i = 1:nbins
        for j = 1:ndepths
            ents(i,j) = shannon(y,bin(i),depth(j))/depth(j);
        end
    end
    % Don't know what quite to do -- I think stick to above, where only one
    % input is a vector at a time.
end


end