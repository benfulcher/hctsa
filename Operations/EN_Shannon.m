function out = EN_Shannon(y, numBin, depth)
% EN_Shannon   Approximate Shannon entropy of a time series.
%
% Uses a numBin-bin encoding and depth-symbol sequences. The series is coarse-grained
% into numBin symbols using uniform population binning (thresholds at equally
% spaced quantiles, so each symbol is used about equally often), every overlapping
% word of depth successive symbols is counted, and the Shannon entropy of the word
% distribution, -sum(p.*log(p)) in nats, is divided by depth to give the entropy
% per symbol (the entropy otherwise scales with depth).
%
% In this wrapper function, you can evaluate the code at a given numBin and depth,
% or across a range of depths (or of numbers of bins) to return statistics on how
% the obtained entropies change.
%
% The implementation uses Michael Small's code MS_shannon.m (renamed from the
% original, simply shannon.m), available at http://small.eie.polyu.edu.hk/matlab/
%
% ---INPUTS:
% y, the input time series
% numBin, the number of bins to discretize the time series into (i.e., alphabet
%    size) (default: 2). Can be a vector (e.g., 2:10) if depth is a single number.
% depth, the length of strings (words) to analyze (default: 3). Can be a vector
%    (e.g., 1:10) if numBin is a single number.
%
% ---OUTPUTS:
% If numBin and depth are both single numbers, a scalar: the Shannon entropy per
% symbol. If one of them is a vector, a structure with fields summarizing the
% entropy per symbol across the range tested:
% maxent, the maximum
% minent, the minimum
% medent, the median
% meanent, the mean
% stdent, the standard deviation
% (Both numBin and depth being vectors is not implemented and gives an error.)
%
% ---REFERENCES:
% M. Small, "Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance", World Scientific, Nonlinear Science Series A, Vol. 52
% (2005).

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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

% -------------------------------------------------------------------------------
% Check inputs:
% -------------------------------------------------------------------------------
if nargin < 2 || isempty(numBin)
	numBin = 2; % two bins to discretize the time series, y
end
if nargin < 3 || isempty(depth)
	depth = 3; % three-long strings
end
binRangeSize = length(numBin);
depthRangeSize = length(depth);

% ------------------------------------------------------------------------------
if binRangeSize == 1
	if depthRangeSize == 1
		%% Evaluate the shannon entropy of discretization
		% Run the code, just return a number
		% This scales with depth, so it's nice to normalize by this factor:
		out = MS_shannon(y, numBin, depth) / depth;
	elseif depthRangeSize > 1
		% Range over depths specified in the vector and return statistics on results
		% (constant number of bins)
		% Somewhat strange behaviour -- very variable
		numDepths = length(depth);
		ents = zeros(numDepths, 1);
		for i = 1:numDepths
			ents(i) = MS_shannon(y, numBin, depth(i)) / depth(i);
		end
		% Output statistics on variation across the range tested:
		out.maxent = max(ents);
		out.minent = min(ents);
		out.medent = median(ents);
		out.meanent = mean(ents);
		out.stdent = std(ents);
	end
elseif binRangeSize > 1
	if depthRangeSize == 1
		%% (*) Statistics over different bin numbers (constant depth)
		% Range over bins specified in the vector numBin; return statistics on results
		ents = zeros(binRangeSize, 1);
		for i = 1:binRangeSize
			ents(i) = MS_shannon(y, numBin(i), depth) / depth;
		end
		out.maxent = max(ents);
		out.minent = min(ents);
		out.medent = median(ents);
		out.meanent = mean(ents);
		out.stdent = std(ents);
	elseif depthRangeSize > 1
		% Don't know what quite to do -- I think stick to above, where only one
		% input is a vector at a time.
		% ***INCOMPLETE*** don't do this.
		error('Comparing both bins and depth not implemented')
		%% (*) stats over numBins and depths
		% ents = zeros(binRangeSize,depthRangeSize);
		% for i = 1:numBins
		%     for j = 1:numDepths
		%         ents(i,j) = MS_shannon(y,numBin(i),depth(j))/depth(j);
		%     end
		% end
	end
end

end
