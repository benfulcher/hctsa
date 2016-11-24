function out = EN_MS_shannon(y,nbin,depth)
% EN_MS_shannon     Approximate Shannon entropy of a time series.
%
% Uses an nbin-bin encoding and depth-symbol sequences.
% Uniform population binning is used, and the implementation uses Michael Small's code
% MS_shannon.m (renamed from the original, simply shannon.m)
%
% cf. M. Small, Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance (book) World Scientific, Nonlinear Science Series A,
% Vol. 52 (2005)
% Michael Small's code is available at available at http://small.eie.polyu.edu.hk/matlab/
%
% In this wrapper function, you can evaluate the code at a given n and d, and
% also across a range of depth and nbin to return statistics on how the obtained
% entropies change.
%
%---INPUTS:
% y, the input time series
% nbin, the number of bins to discretize the time series into (i.e., alphabet size)
% depth, the length of strings to analyze

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

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(nbin)
    nbin = 2; % two bins to discretize the time series, y
end
if nargin < 3 || isempty(depth)
    depth = 3; % three-long strings
end

% ------------------------------------------------------------------------------
%% Evaluate the shannon entropy for a given set of parameters
% ------------------------------------------------------------------------------
if (length(nbin) == 1) && (length(depth) == 1)
    % Run the code, just return a number
    % scales with depth, so it's nice to normalize by this factor:
    out = MS_shannon(y,nbin,depth) / depth;
end

% ------------------------------------------------------------------------------
%% (*) Return statistics over depths (constant number of bins)
% ------------------------------------------------------------------------------
% Somewhat strange behaviour -- very variable
if (length(nbin) == 1) && (length(depth) > 1)
    % Range over depths specified in the vector and return statistics on results
    numDepths = length(depth);
    ents = zeros(numDepths,1);
    for i = 1:numDepths
        ents(i) = MS_shannon(y,nbin,depth(i));
    end
    % Should scale with depth: normalize by this:
    ents = ents./depth';
    out.maxent = max(ents);
    out.minent = min(ents);
    out.medent = median(ents);
    out.meanent = mean(ents);
    out.stdent = std(ents);
end

% ------------------------------------------------------------------------------
%% (*) Statistics over different bin numbers
% ------------------------------------------------------------------------------
if (length(nbin) > 1) && (length(depth) == 1)
    % Range over bins specified in the vector nbin; return statistics on results
    nbins = length(nbin);
    ents = zeros(nbins,1);
    for i = 1:nbins
        ents(i) = MS_shannon(y,nbin(i),depth);
    end
    out.maxent = max(ents);
    out.minent = min(ents);
    out.medent = median(ents);
    out.meanent = mean(ents);
    out.stdent = std(ents);
end

% ------------------------------------------------------------------------------
%% (*) statistics over both nbins and depths
% ------------------------------------------------------------------------------
if (length(nbin) > 1) && (length(depth) > 1)
    nbins = length(nbin);
    numDepths = length(depth);

    ents = zeros(nbins,numDepths);
    for i = 1:nbins
        for j = 1:numDepths
            ents(i,j) = MS_shannon(y,nbin(i),depth(j))/depth(j);
        end
    end
    % Don't know what quite to do -- I think stick to above, where only one
    % input is a vector at a time.
    % ***INCOMPLETE*** don't do this.
    error('Comparing both bins and depth not implemented')
end


end
