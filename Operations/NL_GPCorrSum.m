function out = NL_GPCorrSum(y, Nref, r, thwin, nbins, embedParams, doTwo)
% NL_GPCorrSum correlation sum scaling by Grassberger-Proccacia algorithm
%
% ---INPUTS:
% y, column vector of time-series data
% Nref, number of (randomly-chosen) reference points (-1: use all points,
%       if a decimal, then use this fraction of the time series length)
% r, maximum search radius relative to attractor size, 0 < r < 1
% thwin, number of samples to exclude before and after each reference index
%        (~ Theiler window)
% nbins, number of partitioned bins
% embedParams, embedding parameters to feed BF_Embed.m for embedding the
%               signal in the form {tau,m}
% doTwo, if this is set to 1, will use corrsum, if set to 2, will use corrsum2.
%           For corrsum2, n specifies the number of pairs per bin. Default is 1,
%           to use corrsum.
%
% ---OUTPUTS: basic statistics on the correlation-sum scaling, including
% iteratively re-weighted least squares linear fits to log-log plots using
% the robustfit function in Matlab's Statistics Toolbox.
%
% Uses TISEAN's d2 (rather than TSTOOL's corrsum/corrsum2, which this
% operation used previously) to compute the correlation sum for a
% time-delay-embedded time series by the Grassberger-Procaccia algorithm.
% d2's raw correlation-sum output (its .c2 file, also used by
% NL_TISEAN_d2.m) gives (ln r, C(r)) pairs directly; only doTwo = 1 (the
% corrsum-equivalent binning, log-spaced radii) is supported -- doTwo = 2
% (corrsum2's fixed-pairs-per-bin binning) has no TISEAN equivalent and was
% never used by any of this operation's own mop-file entries.
%
% cf. "Characterization of Strange Attractors", P. Grassberger and I. Procaccia,
% Phys. Rev. Lett. 50(5) 346 (1983)

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Preliminaries
% ------------------------------------------------------------------------------
doPlot = 0; % whether to plot outputs to figure
N = length(y); % length of time series

% (1) Number of reference points, Nref
if nargin < 2 || isempty(Nref)
	Nref = 500; % 500 points
end
if (Nref > 0) && (Nref < 1)
	Nref = round(N * Nref); % specify a proportion of time series length
end
if Nref >= N
	Nref = -1; % Number of reference points capped at time series length
end

% (2) Maximum relative search radius, r
if nargin < 3 || isempty(r)
	r = 0.05; % 5% of attractor radius
end

% (3) Remove spurious correlations of adjacent points, thwin
if nargin < 4 || isempty(thwin)
	thwin = 10; % default window length
end

% (4) Number of bins, nbins
if nargin < 5 || isempty(nbins)
	nbins = 20; % defulat number of bins
end

% (5) Set embedding parameters to defaults
if nargin < 6 || isempty(embedParams)
	embedParams = {'ac', 'fnnmar'};
else
	if length(embedParams) ~= 2
		error('Embedding parameters are formatted incorrectly -- need {tau,m}')
	end
end

if nargin < 7 || isempty(doTwo)
	doTwo = 1; % use corrsum rather than corrsum2
end

if doTwo == 2
	error(['NL_GPCorrSum: doTwo = 2 (corrsum2''s fixed-pairs-per-bin ' ...
		   'binning) has no TISEAN equivalent and is not supported.']);
end

% ------------------------------------------------------------------------------
%% Resolve the embedding parameters (tau, m)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, embedParams{1}, embedParams{2}, true);
tau = tm(1);
m = tm(2);

if (N - (m - 1) * tau) < thwin
	warning(['Embedded time series (N = %u, m = %u, tau = %u) too short' ...
			 ' to do a correlation sum\n'], N, m, tau);
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% Write the file for TISEAN to work with
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);

% Map TSTOOL's Nref convention (-1 = use all points) onto TISEAN's d2
% (-N 0 = use all pairs):
if Nref == -1
	NrefTISEAN = 0;
else
	NrefTISEAN = Nref;
end

% Maximum search radius, in standard deviations of y (matching the
% convention already established for NL_TakensEstimator.m and
% NL_TISEAN_d2.m's takens05, citing Kantz & Schreiber -- TSTOOL's own r was
% instead a proportion of "attractor size", equivalent in spirit but not
% numerically identical), scaled by sqrt(m): a pairwise Euclidean distance
% in an m-dimensional embedding of (roughly independent) per-coordinate
% scale std(y) itself scales as std(y)*sqrt(m), not std(y) alone. Without
% this correction, small r requests at higher embedding dimensions can ask
% for a radius smaller than the closest pairs of points in the actual
% reconstructed phase space, leaving every bin with zero pairs (confirmed
% empirically: for m=7 noise, the smallest radius with any pairs at all
% was ~2.6% away from r*std(y)*sqrt(m), vs. 2.6x too large for r*std(y)
% alone). d2's minimum epsilon is left at its own default (max data
% interval/1000): tightening it in proportion to maxEps was tried and
% found to shrink the useful (non-zero, non-saturated) part of the range
% rather than help. A handful of bins with too few points to fit is
% handled by the "not enough points" check below, same as this operation
% has always done for degenerate cases -- much rarer now, but not
% entirely eliminated: a small requested r can still legitimately fall
% below the closest-pair distance for a genuinely high-dimensional,
% unstructured series (verified: this remains possible for structureless
% noise embedded at a high FNN-selected m, but not for real
% low-dimensional structure, e.g. a chaotic logistic map continues to
% give a sensible estimate at the same small r). In that regime a NaN
% output is an honest "no reliable low-dimensional scaling could be
% detected at this length scale", not a computation failure to paper over.
maxEps = r * std(y) * sqrt(m);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, d2
% ------------------------------------------------------------------------------
% d2 over embedding dimensions 1:m (only the m-th is used below). Note:
% "-M<m>,<m>" (i.e. asking for a single, fixed embedding dimension) triggers
% a bug in this TISEAN build where it reports "0 lines read" and produces no
% output at all; "-M1,<m>" (a genuine range, as NL_TISEAN_d2.m and
% NL_TakensEstimator.m already use) works correctly, so that's used
% here too and the dimension-m block is picked out afterwards.
[~, res] = system(sprintf('d2 -d%u -M1,%u -t%u -R%g -N%u -#%u %s', ...
						  tau, m, thwin, maxEps, NrefTISEAN, nbins, filePath));
if exist([filePath '.stat'], 'file'), delete([filePath '.stat']); end
if exist([filePath '.d2'], 'file'), delete([filePath '.d2']); end
if exist([filePath '.h2'], 'file'), delete([filePath '.h2']); end

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	if exist([filePath '.c2'], 'file'), delete([filePath '.c2']); end
	error('Call to TISEAN function ''d2'' failed.');
end

% ------------------------------------------------------------------------------
%% Parse the raw correlation-sum data (.c2 file)
% ------------------------------------------------------------------------------
if ~exist([filePath '.c2'], 'file')
	error('TISEAN function ''d2'' did not produce a .c2 output file.');
end
fid = fopen([filePath '.c2']);
fileLines = textscan(fid, '%[^\n]');
fclose(fid);
delete([filePath '.c2']);
fileLines = fileLines{1};

% First line is a '#center=' header (not per-dimension); the rest is one
% '#dim=' block per embedding dimension 1:m -- take the last one (m):
w = strmatch('#dim=', fileLines);
if length(w) ~= m
	error('TISEAN function ''d2'' returned an unexpected number of data blocks.');
end
w(end + 1) = length(fileLines) + 1;
ss = fileLines(w(m) + 1:w(m + 1) - 1);

rc = zeros(length(ss), 2); % [ln(r), C(r)]
nn = 0;
for jj = 1:length(ss)
	tmp = textscan(ss{jj}, '%f%f');
	if all(cellfun(@isempty, tmp))
		break % a trailing comment line
	end
	nn = nn + 1;
	rc(nn, :) = horzcat(tmp{:});
end
rc = rc(1:nn, :);

if isempty(rc)
	fprintf(1, 'No output obtained from d2''s correlation sum.\n');
	out = NaN; return
end

% The .c2 file stores raw r (log-spaced, not pre-logged -- confirmed by the
% consecutive ratios of column 1 being ~constant, i.e. geometric, whereas
% ln(r) would be arithmetic/evenly-spaced) and raw C(r) (not ln(C(r))):
lnr = log(rc(:, 1));
lnCr = log(rc(:, 2));

if doPlot
	figure('color', 'w'); box('on');
	plot(lnr, lnCr, '.-k');
end

% Contains ln(r) in rows and values are ln(C(r));

% ------------------------------------------------------------------------------
%% Remove any Infs in lnCr
% ------------------------------------------------------------------------------
rGood = (isfinite(lnCr));
if ~any(rGood)
	fprintf(1, 'No good outputs obtained from the correlation sum.\n');
	out = NaN; return
end

% Only keep finite values:
lnCr = lnCr(rGood);
lnr = lnr(rGood);

% ------------------------------------------------------------------------------
%% Output Statistics
% ------------------------------------------------------------------------------
% Basic statistics:
out.minlnr = min(lnr);
out.maxlnr = max(lnr);
out.minlnCr = min(lnCr);
out.maxlnCr = max(lnCr);
out.rangelnCr = range(lnCr);
out.meanlnCr = mean(lnCr);

% Fit linear to log-log plot (full range)
enoughpoints = true;
try
	[a, stats] = robustfit(lnr, lnCr);
catch me
	if strcmp(me.message, 'Not enough points to perform robust estimation.')
		enoughpoints = false;
	end
end

if enoughpoints
	out.robfit_a1 = a(1);
	out.robfit_a2 = a(2);
	out.robfit_sigrat = stats.ols_s / stats.robust_s;
	out.robfit_s = stats.s;
	out.robfit_sea1 = stats.se(1);
	out.robfit_sea2 = stats.se(2);

	fit_lnCr = a(2) * lnr + a(1);
	if doPlot, hold on; plot(lnr, fit_lnCr, 'r'); hold off; end

	% Compute residuals (lnCr and fit_lnCr are both column vectors here; a
	% stray transpose on this line -- present in the original TSTOOL-based
	% version of this operation too -- turned this into an N-by-N matrix
	% via implicit broadcasting instead of an N-by-1 residual vector,
	% making robfitresmeanabs/robfitresmeansq vectors instead of scalars):
	res = lnCr - fit_lnCr;

	out.robfitresmeanabs = mean(abs(res));
	out.robfitresmeansq = mean(res.^2);
	out.robfitresac1 = CO_AutoCorr(res, 1, 'Fourier');
else
	out.robfit_a1 = NaN;
	out.robfit_a2 = NaN;
	out.robfit_sigrat = NaN;
	out.robfit_s = NaN;
	out.robfit_sea1 = NaN;
	out.robfit_sea2 = NaN;

	out.robfitresmeanabs = NaN;
	out.robfitresmeansq = NaN;
	out.robfitresac1 = NaN;
end

% now non-robust linear fit
% [p, S] = polyfit(lnr',lnCr,1);

end
