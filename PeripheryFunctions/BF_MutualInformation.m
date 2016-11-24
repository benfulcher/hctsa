function mi = BF_MutualInformation(v1,v2,r1,r2,numBins)
% BF_MutualInformation  Mutual information between two data vectors using bin counting.
%
% Mutual information computed using a histogram-based, bin-counting method.
%
%---INPUTS:
% v1, the first input vector
% v2, the second input vector
% r1, the bin-partitioning method for the first input vector, v1
% r2, the bin-partitioning method for the second input vector, v2
% numBins, the number of bins to partition each vector into.
%
% NB: r1 and r2 can also be two-component vectors, that specify a custom range
%     for binning
%
%---OUTPUT:
% mi, the mutual information computed between v1 and v2

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
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
% by default, take a range equal to the range of the vectors
if nargin < 3 || isempty(r1)
    r1 = 'range';
end
if nargin < 4 || isempty(r2)
    r2 = 'range';
end
% use 10 bins for the histograms (=100 bins in 2d)
if nargin < 5 || isempty(numBins)
    numBins = 10;
end
% ------------------------------------------------------------------------------

if length(v1) ~= length(v2)
    error('input vectors must be the same length');
end

N = length(v1); % length of vectors

% Make sure both column vectors
if size(v1,2) > size(v1,1), v1 = v1'; end
if size(v2,2) > size(v2,1), v2 = v2'; end

% Create histograms:
% (i) in v1
edgesi = SUB_GiveMeEdges(r1,v1,numBins);
[ni, bini] = histc(v1, edgesi);

% (ii) in v2
edgesj = SUB_GiveMeEdges(r2,v2,numBins);
[nj, binj] = histc(v2, edgesj);

% ------------------------------------------------------------------------------
%% Create a joint histogram:
% ------------------------------------------------------------------------------

% We have the edges in each dimension: edgesi, and edgesj
histxy = zeros(numBins);
for i2 = 1:numBins
    for j2 = 1:numBins
        histxy(i2,j2) = sum(bini==i2 & binj==j2);
    end
end

% Normalize counts to probabilities
p_i = ni(1:numBins)/N;
p_j = nj(1:numBins)/N;
p_ij = histxy/N;
p_ixp_j = p_i*p_j';
summe = (p_ixp_j > 0 & p_ij > 0);

% Do a matrix-sum mutual information calculation:
if any(summe(:) == 1)
    mi = sum(p_ij(summe).*log(p_ij(summe)./p_ixp_j(summe)));
else
    fprintf(1,['The histograms aren''t catching any points?? ' ...
            'Perhaps due to an inappropriate custom range for binning the data...\n']);
    mi = NaN; return
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function edges = SUB_GiveMeEdges(r,v,nbins)
    EE = 1E-6; % this small addition gets lost in the last bin
    if strcmp(r,'range')
        edges = linspace(min(v),max(v)+EE,nbins+1);

    elseif strcmp(r,'quantile') % bin edges based on quantiles
        edges = quantile(v,linspace(0,1,nbins+1));
%             edges(1) = edges(1) - 0.1;
        edges(end) = edges(end) + EE;
%             edges = sort(unique(edges)); % in case you have many repeated values -- will bias MI calculation
    elseif length(r)==2 % a two-component vector
        edges = linspace(r(1),r(2)+EE,nbins+1);
    else
        error('Unknown partitioning method ''%s''',r);
    end
end
% ------------------------------------------------------------------------------
end
