function out = FC_Surprise(y,whatPrior,memory,numGroups,cgmeth,numIters,randomSeed)
% FC_Surprise   How surprised you would be of the next data point given recent memory.
%
% Coarse-grains the time series, turning it into a sequence of symbols of a
% given alphabet size, numGroups, and quantifies measures of surprise of a
% process with local memory of the past memory values of the symbolic string.
%
% We then consider a memory length, memory, of the time series, and
% use the data in the proceeding memory samples to inform our expectations of
% the following sample.
%
% The 'information gained', log(1/p), at each sample using expectations
% calculated from the previous memory samples, is estimated.
%
%---INPUTS:
% y, the input time series
%
% whatPrior, the type of information to store in memory:
%           (i) 'dist': the values of the time series in the previous memory
%                       samples,
%           (ii) 'T1': the one-point transition probabilities in the previous
%                       memory samples, and
%           (iii) 'T2': the two-point transition probabilities in the previous
%                       memory samples.
%
% memory, the memory length (either number of samples, or a proportion of the
%           time-series length, if between 0 and 1)
%
% numGroups, the number of groups to coarse-grain the time series into
%
% cgmeth, the coarse-graining, or symbolization method:
%          (i) 'quantile': an equiprobable alphabet by the value of each
%                          time-series datapoint,
%          (ii) 'updown': an equiprobable alphabet by the value of incremental
%                         changes in the time-series values, and
%          (iii) 'embed2quadrants': 4-letter alphabet of the quadrant each data
%                            point resides in a two-dimensional embedding space.
%
% numIters, the number of iterations to repeat the procedure for.
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%
%---OUTPUTS: summaries of this series of information gains, including the
%            minimum, maximum, mean, median, lower and upper quartiles, and
%            standard deviation.

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

% ------------------------------------------------------------------------------
%% Check inputs and set defaults
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(whatPrior)
    whatPrior = 'dist'; % expect probabilities based on prior observed distribution
end

% memory: how far into the past to base your priors on
if nargin < 3 || isempty(memory)
    memory = 0.2; % set it as 20% of the time-series length
end
if (memory > 0) && (memory < 1) % specify memory as a proportion of the time-series length
    memory = round(memory*length(y));
end

% ng -- number of groups for the time-series coarse-graining/symbolization
if nargin < 4 || isempty(numGroups)
    numGroups = 3; % use three symbols to approximate the time-series values
end

% cgmeth: the coase-graining method to use
if nargin < 5 || isempty(cgmeth)
    cgmeth = 'quantile'; % symbolize time series by their values (quantile)
end

% numIters: number of iterations
if nargin < 6 || isempty(numIters)
    numIters = 500;
    % number of iterations of the procedure to perform (does it with random samples)
    % could also imagine doing it exhaustively...?!
end

% randomSeed: how to treat the randomization
if nargin < 7
    randomSeed = []; % default for BF_ResetSeed
end

% ------------------------------------------------------------------------------
%% Course Grain
% ------------------------------------------------------------------------------
yth = SB_CoarseGrain(y,cgmeth,numGroups); % a coarse-grained time series using the numbers 1:numGroups

N = length(yth); % will be the same as y, for 'quantile', and 'updown'

% Select random samples to test:
BF_ResetSeed(randomSeed); % control random seed (for reproducibility)
rs = randperm(N-memory) + memory; % Can't do beginning of time series, up to memory
rs = sort(rs(1:min(numIters,end))); % Just use a random sample of numIters points to test

%-------------------------------------------------------------------------------
% Compute empirical probabilities from time series
%-------------------------------------------------------------------------------
store = zeros(numIters,1); % store probabilities
for i = 1:length(rs)
    switch whatPrior
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
            memoryData = yth(rs(i)-memory:rs(i)-1);
            % Previous value observed in memory here:
            inmem = find(memoryData(1:end-1) == yth(rs(i)-1));
            if isempty(inmem)
                p = 0;
            else
                p = mean(memoryData(inmem+1) == yth(rs(i)));
            end
            store(i) = p;

        case 'T2'
            % Uses two-point correlations in memory to inform the next point

            memoryData = yth(rs(i)-memory:rs(i)-1);
            % Previous value observed in memory here:
            inmem1 = find(memoryData(2:end-1) == yth(rs(i)-1)); % the 2:end makes the next line ok...?
            inmem2 = find(memoryData(inmem1) == yth(rs(i)-2));
            if isempty(inmem2)
                p = 0;
            else
                p = sum(memoryData(inmem2+2) == yth(rs(i)))/length(inmem2);
            end
            store(i) = p;

        otherwise
            error('Unknown method ''%s''',whatPrior);
    end
end

%-------------------------------------------------------------------------------
% Information gained from next observation is log(1/p) = -log(p)
%-------------------------------------------------------------------------------
store(store==0) = 1; % so that we set log(0)==0
store = -log(store); % transform to surprises/information gains
% histogram(store)

if any(store > 0)
    out.min = min(store(store > 0)); % Minimum amount of information you can gain in this way
else
    out.min = NaN;
end
out.max = max(store); % Maximum amount of information you can gain in this way
out.mean = mean(store); % mean
out.sum = sum(store); % sum (contains same information as mean since length(store) the same)
out.median = median(store); % median
out.lq = quantile(store,0.25); % lower quartile
out.uq = quantile(store,0.75); % upper quartile
out.std = std(store); % standard deviation

% t-statistic to information gain of 1
if out.std==0
    out.tstat = NaN; % can't compute this if there is no variation
else
    out.tstat = abs((out.mean-1)/(out.std/sqrt(numIters)));
end

end
