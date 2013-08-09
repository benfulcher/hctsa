% FC_Surprise
% 
% How surprised you might be by the next recorded data points given the data recorded
% in recent memory.
% 
% Coarse-grains the time series, turning it into a sequence of symbols of a
% given alphabet size, ng, and quantifies measures of surprise of a
% process with local memory of the past memory values of the symbolic string.
% 
% We then consider a memory length, memory, of the time series, and
% use the data in the proceeding memory samples to inform our expectations of
% the following sample.
% 
% The 'information gained', log(1/p), at each sample using expectations
% calculated from the previous memory samples, is estimated.
% 
% INPUTS:
% y, the input time series
% 
% whatinf, the type of information to store in memory:
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
% ng, the number of groups to coarse-grain the time series into
% 
% cgmeth, the coarse-graining, or symbolization method:
%          (i) 'quantile': an equiprobable alphabet by the value of each
%                          time-series datapoint,
%          (ii) 'updown': an equiprobable alphabet by the value of incremental
%                         changes in the time-series values, and
%          (iii) 'embed2quadrants': by the quadrant each data point resides in
%                          in a two-dimensional embedding space.
% 
% nits, the number of iterations to repeat the procedure for.
% 
% Outputs of this operation are summaries of this series of information gains,
% including the minimum, maximum, mean, median, lower and upper quartiles, and
% standard deviation.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = FC_Surprise(y,whatinf,memory,ng,cgmeth,nits)
% Ben Fulcher, September 2009

%% Check inputs and set defaults
if nargin < 2 || isempty(whatinf)
    whatinf = 'dist'; % expect probabilities based on prior observed distribution
end
% Memory: how far into the past to base your priors on
if nargin < 3 || isempty(memory)
    memory = 0.2; % set it as 20% of the time-series length
end
if (memory > 0) && (memory < 1) % specify memory as a proportion of the time-series length
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
yth = SB_coarsegrain(y,cgmeth,ng); % a coarse-grained time series using the numbers 1:ng

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
                p = sum(memorydata(inmem2+2) == yth(rs(i)))/length(inmem2);
            end
            store(i) = p;
            
        otherwise
            error('Unknwon method ''%s''',whatinf);
    end
end

% information gained from next observation is log(1/p) = -log(p)
iz = (store == 0);
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