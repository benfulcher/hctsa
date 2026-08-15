function out = PP_SchreiberDenoise(y, m, d, numIter, neighborhoodStd)
% PP_SchreiberDenoise   Nonlinear noise reduction, and how it changes the series.
%
% Applies Schreiber's simple nonlinear noise-reduction method -- replace each
% embedded point by the average of its state-space neighbors -- via TISEAN's
% 'nrlazy' (the whole-vector-correcting C reimplementation of the original
% single-component-correcting Fortran 'lazy'; TISEAN's own documentation
% reports nrlazy tends to do better than lazy on flow-like data, which is
% most of what hctsa sees):
%
% Schreiber, T. "Extremely simple nonlinear noise reduction method",
% Phys. Rev. E 47, 2401 (1993).
%
% This is the denoising stage of the same Chaos Decision Tree Algorithm
% pipeline as CO_Oversampling (its oversampling-diagnosis stage) and
% NL_ZeroOneTest (its chaos-classification stage):
%
% Toker, D. et al. "A simple method for detecting chaos in nature",
% Commun. Biol. 3, 11 (2020). DOI: 10.1038/s42003-019-0715-9
%
% Denoising is not itself a scalar feature, so -- following the same idiom
% as PP_Compare/PP_ModelFit -- this reports how several time-series
% properties change as a result of applying it:
% - rmsCorrection, corrOrigDenoised: how much/little the series changed.
% - fracVarRemoved: the fraction of variance attributed to noise (removed).
% - ac1Change: change in lag-1 autocorrelation (denoising should smooth the
%   series, increasing it).
% - meanNeighbors, fracNoCorrection: nrlazy's own per-point neighbor counts
%   (a neighbor count of 1 means no correction was possible/applied at that
%   point -- diagnoses whether neighborhoodStd was too small for this data).
% - KDenoised, KChange: this pipeline's actual point in denoising first is
%   to make NL_ZeroOneTest's 0-1 test for chaos more reliable on noisy data
%   (measurement noise inflates the apparent diffusion of the test's
%   (p,q) trajectory) -- KDenoised is NL_ZeroOneTest's K statistic computed
%   on the denoised series, and KChange = KDenoised - K(y) directly answers
%   "does removing noise change the chaos verdict for this series?" (K(y)
%   itself isn't reported here since it's already registered directly via
%   NL_ZeroOneTest -- reporting it again here would be a pure duplicate).
%
% ---INPUTS:
% y, the input time series
% m, the embedding dimension (nrlazy '-m', default 5)
% d, the embedding delay (nrlazy '-d', default 1)
% numIter, the number of correction passes (nrlazy '-i'; more iterations
%          denoise more aggressively but risk distorting real dynamics --
%          TISEAN's own default, and most published use, is 1)
% neighborhoodStd, neighborhood radius in units of the data's standard
%          deviation (nrlazy '-v'; scale-invariant, unlike the raw '-r'
%          alternative which is a fixed data-interval fraction)
%
% ---OUTPUTS:
% rmsCorrection, fracVarRemoved, corrOrigDenoised, ac1Change, meanNeighbors,
% fracNoCorrection, KDenoised, KChange (see above)
%
% Redundancy check (Empirical1000, m=5,d=1,v=0.3, numIter in {1,4}): all
% fields stay well clear of the usual r>=0.9 threshold except two borderline,
% inconsistent-across-numIter cases -- fracNoCorrection (r=0.90 at numIter=1
% vs NL_TISEAN_fnn's neighborhood-size fields, but only 0.89 at numIter=4)
% and KDenoised (r=0.90 at numIter=1 vs SP_Summaries' spectral log-log slope
% fields, but only 0.89 at numIter=4) -- both sensible mechanistic overlaps
% (local embedding-space density; broadband-vs-periodic structure) rather
% than duplicates, and too marginal/inconsistent to drop.

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

% ------------------------------------------------------------------------------
%% Check inputs & set defaults:
% ------------------------------------------------------------------------------
y = y(:);
N = length(y);
if nargin < 2 || isempty(m)
	m = 5;
end
if nargin < 3 || isempty(d)
	d = 1;
end
if nargin < 4 || isempty(numIter)
	numIter = 1;
end
if nargin < 5 || isempty(neighborhoodStd)
	neighborhoodStd = 0.5;
end

% Need enough points for (m-1)*d + 1 to fit an embedding vector, with margin
% for the local-neighborhood statistics nrlazy computes to be meaningful:
if N < max(50, 10 * ((m - 1) * d + 1))
	warning('Time series too short to denoise at m=%u, d=%u', m, d);
	out.rmsCorrection = NaN;
	out.fracVarRemoved = NaN;
	out.corrOrigDenoised = NaN;
	out.ac1Change = NaN;
	out.meanNeighbors = NaN;
	out.fracNoCorrection = NaN;
	out.KDenoised = NaN;
	out.KChange = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Run nrlazy:
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);

% -V4: stdout carries both the denoised value and the number of neighbors
% found for each point (two columns), rather than just the value alone.
% -m takes two numbers, 'components,embeddingDim' -- always 1 component here
% since y is a scalar (univariate) time series.
tiseanCommand = sprintf('nrlazy -m1,%u -d%u -i%u -v%g -V4 %s', m, d, numIter, neighborhoodStd, filePath);
[~, res] = BF_TiseanSystem(tiseanCommand);

if isempty(res)
	warning('No output from TISEAN routine nrlazy on the data');
	out.rmsCorrection = NaN;
	out.fracVarRemoved = NaN;
	out.corrOrigDenoised = NaN;
	out.ac1Change = NaN;
	out.meanNeighbors = NaN;
	out.fracNoCorrection = NaN;
	out.KDenoised = NaN;
	out.KChange = NaN;
	return
end
data = textscan(res, '%f%f');
yDenoised = data{1};
numNeighbors = data{2};

if length(yDenoised) ~= N || any(isnan(yDenoised))
	warning('nrlazy produced unusable output for this data (m=%u, d=%u too large for the series length?)', m, d);
	out.rmsCorrection = NaN;
	out.fracVarRemoved = NaN;
	out.corrOrigDenoised = NaN;
	out.ac1Change = NaN;
	out.meanNeighbors = NaN;
	out.fracNoCorrection = NaN;
	out.KDenoised = NaN;
	out.KChange = NaN;
	return
end

% ------------------------------------------------------------------------------
%% How much/how did the series change?
% ------------------------------------------------------------------------------
correction = y - yDenoised;
out.rmsCorrection = sqrt(mean(correction.^2));
out.fracVarRemoved = 1 - var(yDenoised) / var(y);
cc = corrcoef(y, yDenoised);
out.corrOrigDenoised = cc(1, 2);
out.ac1Change = CO_AutoCorr(yDenoised, 1, 'Fourier') - CO_AutoCorr(y, 1, 'Fourier');

out.meanNeighbors = mean(numNeighbors);
out.fracNoCorrection = mean(numNeighbors == 1);

% ------------------------------------------------------------------------------
%% Does denoising change the 0-1-test chaos verdict? (the actual point of
%% denoising first in the source pipeline)
% ------------------------------------------------------------------------------
zeroOneOrig = NL_ZeroOneTest(y, 20);
zeroOneDenoised = NL_ZeroOneTest(yDenoised, 20);
out.KDenoised = zeroOneDenoised.K;
out.KChange = zeroOneDenoised.K - zeroOneOrig.K;

end
