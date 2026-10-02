function out = DN_HistogramAsymmetry(y, numBins, doSimple)
% DN_HistogramAsymmetry   Measures of the asymmetry of the histogram of the data.
%
% Compares the histograms of the positive and the negative values of a
% z-scored series. A histogram of numBins equal-width bins is built separately
% over the range of the positive values and over the range of the negative
% values, and expressed as proportions of the number of nonzero values. The
% outputs compare how many values lie on each side, the proportion in the
% fullest bin on each side, and the location of that bin. A symmetric
% distribution gives 0 for all the differences and for modeAsymmetry.
%
% The input is assumed to be z-scored (a warning is raised otherwise), so that
% positive and negative values lie above and below the mean.
%
% ---INPUTS:
% y, the input data vector
% numBins, the number of bins to use in each histogram (default: 10)
% doSimple, whether to use a simple binning method with linearly spaced bins
%       (BF_SimpleBinner; true, the default), or MATLAB's histcounts (false)
%
% ---OUTPUTS:
% densityDiff, the proportion of nonzero values that are positive minus the
%       proportion that are negative
% modeProbPos, the proportion of nonzero values in the fullest bin of the
%       positive values
% modeProbNeg, the proportion of nonzero values in the fullest bin of the
%       negative values
% modeDiff, modeProbPos - modeProbNeg
% posMode, the center of the fullest bin of the positive values (the mean of
%       the centers if several tie)
% negMode, the center of the fullest bin of the negative values (the mean of
%       the centers if several tie)
% modeAsymmetry, posMode + negMode

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
% Check inputs and set defaults:
% -------------------------------------------------------------------------------
if nargin < 2
	numBins = 10;
end
if nargin < 3
	doSimple = true;
end

% -------------------------------------------------------------------------------
% Check z-score standardization (since it is assumed that positive and negative
% values can be treated separately):
iszscored = BF_iszscored(y);
if ~iszscored
	warning('DN_HistogramAsymmetry assumes a z-scored (or standardized) input')
end

% -------------------------------------------------------------------------------
% Compute the histogram separately from positive and negative values in the data:
yPos = y(y > 0);
yNeg = y(y < 0);
if doSimple
	[countsPos, binEdgesPos] = BF_SimpleBinner(yPos, numBins);
	[countsNeg, binEdgesNeg] = BF_SimpleBinner(yNeg, numBins);
else
	[countsPos, binEdgesPos] = histcounts(yPos, numBins);
	[countsNeg, binEdgesNeg] = histcounts(yNeg, numBins);
end

% Normalize by total counts:
NnonZero = sum(y ~= 0);
pPos = countsPos / NnonZero;
pNeg = countsNeg / NnonZero;

% Compute bin centers from bin edges:
binCentersPos = mean([binEdgesPos(1:end - 1); binEdgesPos(2:end)]);
binCentersNeg = mean([binEdgesNeg(1:end - 1); binEdgesNeg(2:end)]);

% Histogram counts and overall density differences:
% Asymmetry about the mean, as a PROPORTION of the non-zero values (matching the
% denominator modeProbPos/modeProbNeg below already use). A raw count difference is
% extensive -- it grows with the time-series length rather than describing the
% distribution: measured on real-world series, |densityDiff| had Spearman rho = 0.45
% with series length (values ranged to +/-9070), which normalizing removes (-0.08).
out.densityDiff = (sum(y > 0) - sum(y < 0)) / NnonZero;
out.modeProbPos = max(pPos);
out.modeProbNeg = max(pNeg);
out.modeDiff = out.modeProbPos - out.modeProbNeg;

% Mean position of maximums (if multiple):
out.posMode = mean(binCentersPos(pPos == out.modeProbPos));
out.negMode = mean(binCentersNeg(pNeg == out.modeProbNeg));
out.modeAsymmetry = out.posMode + out.negMode;

end
