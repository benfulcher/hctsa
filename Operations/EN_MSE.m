function out = EN_MSE(y, scaleRange, m, r, preProcessHow, whatEntropy, numClasses)
% EN_MSE   Multiscale entropy of a time series.
%
% At each scale s the time series is coarse-grained by averaging over
% non-overlapping windows of s samples (scale 1 is the original series), and the
% entropy of the coarse-grained series is computed: by default the sample entropy,
% SampEn(m,r) (EN_SampEn), as in the multiscale entropy of Costa et al. All s
% possible starting offsets of the windows are used, instead of just offset 0 (see
% NOTES). The sample entropy at scale s is the refined composite multiscale
% entropy (RCMSE) of Wu et al. (2014): the template matches of length m and m+1 are
% counted for each of the s coarse-grained series (as in EN_SampEn), and the value
% is -log of the summed matches of length m+1 over the summed matches of length m.
% Offsets for which the coarse-grained series has fewer than 20 samples are
% skipped, and a scale with no remaining offsets, or with no matches, is NaN. The
% entropy is also summarized across scales (extremes and where they occur, mean,
% spread, trend).
%
% ---INPUTS:
% y, the input time series
% scaleRange, a vector of scales (default: 1:10)
% m, the embedding dimension (length of sequence to match) (default: 2)
% r, the similarity threshold for matching, an absolute value, so it is not
%    rescaled with the scale (default: 0.15). It is a fraction
%    of the standard deviation of the input if y is z-scored (as in hctsa).
% preProcessHow, how to preprocess the data before coarse-graining (default: do
%    not). The preprocessed series (BF_PreProcess) is z-scored. Options include
%    'diff1' (incremental differences) and 'rescale_tau' (first coarse-grain at
%    the 1/e autocorrelation time, BF_GetTau(y,'ac1e')).
% whatEntropy, which entropy to evaluate at each scale:
%    'sampen' (default): sample entropy, the classical multiscale entropy of Costa
%        et al.
%    'dispen': normalized dispersion entropy (EN_DispEn, with tau = 1), i.e.,
%        multiscale dispersion entropy (MDE); cf. H. Azami, M. Rostaghi,
%        D. Abasolo, J. Escudero, "Refined Composite Multiscale Dispersion
%        Entropy and its Application to Biomedical Signals", IEEE Trans.
%        Biomed. Eng. 64(12) 2872 (2017).
%    'fdispen': the fluctuation-based variant of the same.
%    r is unused for the dispersion-based settings (they partition amplitude into
%    classes rather than applying a distance tolerance), and output field names
%    carry the corresponding suffix (dispen_s1, meanDispEn, ...).
% numClasses, the number of amplitude classes for the dispersion-based settings
%    (default: 6, EN_DispEn's own default); unused for 'sampen'.
%
% ---OUTPUTS:
% A structure with fields (for the sample entropy setting; for 'dispen' and
% 'fdispen' the entropy-specific names change as given in brackets):
% sampen_s1, sampen_s2, sampen_s3, sampen_s4, sampen_s5, sampen_s6, sampen_s7,
%     sampen_s8, sampen_s9, sampen_s10 [dispen_s1, ... or fdispen_s1, ...], the
%     entropy at each scale in scaleRange (named by the scale; these names are for
%     the default scaleRange)
% The following summaries across scales are NaN if the entropy at any scale is NaN:
% maxSampEn [maxDispEn, maxFDispEn], the maximum entropy across scales
% maxScale, the scale at which the maximum occurs
% minSampEn [minDispEn, minFDispEn], the minimum entropy across scales
% minScale, the scale at which the minimum occurs
% meanSampEn [meanDispEn, meanFDispEn], the mean entropy across scales
% stdSampEn [stdDispEn, stdFDispEn], the standard deviation of the entropy across
%     scales
% cvSampEn [cvDispEn, cvFDispEn], the coefficient of variation of the entropy
%     across scales (stdSampEn/meanSampEn)
% meanch, the mean change in entropy from one scale to the next
% slope, the slope of a robust linear fit (robustfit) of entropy against scale
%     (NaN for fewer than 3 scales)
% slopeSE, the standard error of that slope
% A scalar NaN is returned instead if no scale has enough samples.
%
% ---REFERENCES:
% Costa, Goldberger and Peng, "Multiscale entropy analysis of biological
% signals", Phys. Rev. E 71, 021906 (2005).
% http://physionet.comp.nus.edu.sg/physiotools/mse/papers/pre-2005.pdf
%
% Original C implementation and docs here:
% http://physionet.org/physiotools/mse/tutorial/node3.html
%
% S.-D. Wu, C.-W. Wu, S.-G. Lin, C.-C. Wang and K.-Y. Lee, "Time Series Analysis
% Using Composite Multiscale Entropy", Entropy 15(3), 1069-1084 (2013).
% DOI: 10.3390/e15031069
%
% S.-D. Wu, C.-W. Wu, S.-G. Lin, K.-Y. Lee and C.-K. Peng, "Analysis of complex
% time series using refined composite multiscale entropy", Physics Letters A 378,
% 1369-1374 (2014), Eq. (8).
% DOI: 10.1016/j.physleta.2014.03.034
%
% ---NOTES:
% The dispersion setting estimates a similar quantity to SampEn, but more precisely
% and much faster. The advantage reverses only under extreme quantization (very few
% distinct values), where the amplitude classes become degenerate and SampEn is
% better.
%
% Composite coarse-graining: the classic Costa coarse-graining at scale s uses a
% single, arbitrary non-overlapping partition of the series starting at the first
% sample, but s-1 other equally valid partitions exist (starting at offsets
% 1,...,s-1), and which one is used has no scientific meaning. For short series,
% the single-offset SampEn at large scales can depend strongly on this arbitrary
% offset. All s offsets are therefore used. For the sample entropy, the match
% counts are pooled across offsets before the log (refined composite MSE, Wu et
% al. 2014), so the value is defined even if some offsets have no matches of length
% m+1 (for which their own SampEn is undefined). For the dispersion
% settings, the value at scale s is the mean over the offsets (composite MSE, Wu et
% al. 2013).
%
% Registered summaries: for the dispersion setting, only meanDispEn, minDispEn and
% stdDispEn are registered as features; the remaining summaries were less reliable
% across independent realizations of the same process. For the sample entropy setting,
% maxSampEn is not registered (largely redundant with meanSampEn, which is essentially
% Costa's "Complexity Index"), whereas maxScale (the location of the peak, distinct
% from its height) is. meanch is also not registered: slope and slopeSE summarize the
% same trend more completely, using all scales at once.

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
% Check inputs, set defaults
% -------------------------------------------------------------------------------
if nargin < 2 || isempty(scaleRange)
	scaleRange = 1:10;
end
if nargin < 3 || isempty(m)
	m = 2;
end
if nargin < 4 || isempty(r)
	r = 0.15;
end
if nargin < 5
	preProcessHow = '';
end

if nargin < 6 || isempty(whatEntropy)
	whatEntropy = 'sampen'; % the classical Costa et al. multiscale entropy
end
if ~ismember(whatEntropy, {'sampen', 'dispen', 'fdispen'})
	error('Unknown entropy ''%s'' (expected ''sampen'', ''dispen'' or ''fdispen'')', whatEntropy);
end

if nargin < 7 || isempty(numClasses)
	numClasses = 6; % EN_DispEn's own default; unused for 'sampen'
end

% Field-name suffix (the 'sampen' outputs use the plain SampEn/sampen names):
switch whatEntropy
	case 'sampen',  enName = 'SampEn'; enPrefix = 'sampen';
	case 'dispen',  enName = 'DispEn'; enPrefix = 'dispen';
	case 'fdispen', enName = 'FDispEn'; enPrefix = 'fdispen';
end

% -------------------------------------------------------------------------------
% Impose a minimum time-series length of 20 samples to perform a SampEn
% (should probably be even higher...?)
minTSLength = 20;

doPlot = false; % whether to plot outputs
numScales = length(scaleRange);

% -------------------------------------------------------------------------------
% Preprocess
% -------------------------------------------------------------------------------
% Do the specified pre-processing BEFORE applying the coarse-graining
if ~isempty(preProcessHow)
	y = zscore(BF_PreProcess(y, preProcessHow));
end

% -------------------------------------------------------------------------------
% Coarse-graining and entropy across scales:
% -------------------------------------------------------------------------------
% cf. Eq. (16) in Costa et al. (2005) for the coarse-graining at offset 0; here all
% `scale` possible non-overlapping starting offsets are used (see NOTES above).
sampEns = nan(numScales, 1);
for si = 1:numScales
	scale = scaleRange(si);
	switch whatEntropy
		case 'sampen'
			% Refined composite multiscale entropy (Wu et al., 2014): pool the
			% template-match counts over the offsets before taking the log:
			numMatchesM1 = 0; % matches at length m+1
			numMatchesM = 0; % matches at length m
			for off = 0:scale - 1
				y_cg = mean(BF_MakeBuffer(y(off + 1:end), scale), 2);
				if length(y_cg) >= minTSLength
					[~, countA, countB] = sampen_mex(y_cg', m + 1, r);
					numMatchesM1 = numMatchesM1 + countA(m + 1);
					numMatchesM = numMatchesM + countB(m + 1);
				end
			end
			if numMatchesM1 > 0 && numMatchesM > 0
				sampEns(si) = -log(numMatchesM1 / numMatchesM);
			end
		otherwise % 'dispen' / 'fdispen': mean over offsets (composite)
			offsetEns = nan(scale, 1);
			for off = 0:scale - 1
				y_cg = mean(BF_MakeBuffer(y(off + 1:end), scale), 2);
				if length(y_cg) >= minTSLength
					dispEnStruct = EN_DispEn(y_cg, m, numClasses, 1);
					if isstruct(dispEnStruct)
						if strcmp(whatEntropy, 'dispen')
							offsetEns(off + 1) = dispEnStruct.normDispEn;
						else
							offsetEns(off + 1) = dispEnStruct.normFDispEn;
						end
					end
				end
			end
			sampEns(si) = mean(offsetEns, 'omitnan');
	end
end

% -------------------------------------------------------------------------------
% Outputs: multiscale entropy
% -------------------------------------------------------------------------------
if all(isnan(sampEns))
	if ~isempty(preProcessHow)
		ppText = sprintf('after %s pre-processing', preProcessHow);
	else
		ppText = '';
	end
	warning('Not enough samples (%u %s) to compute %s at multiple scales', ...
			length(y), ppText, enName)
	out = NaN;
	return
end

if doPlot
	figure('color', 'w') %#ok<UNRCH>
	subplot(2, 1, 1);
	plot(y);
	subplot(2, 1, 2);
	plot(sampEns, 'o-k')
end

% Output raw values
for i = 1:numScales
	out.(sprintf('%s_s%u', enPrefix, scaleRange(i))) = sampEns(i);
end

% -------------------------------------------------------------------------------
% Summary statistics of the variation across scales (NaN if any scale is NaN):
% -------------------------------------------------------------------------------
summaryFields = {['max' enName], 'maxScale', ['min' enName], 'minScale', ...
				['mean' enName], ['std' enName], ['cv' enName], 'meanch', 'slope', 'slopeSE'};
if any(isnan(sampEns))
	for i = 1:length(summaryFields)
		out.(summaryFields{i}) = NaN;
	end
	return
end

% Maximum, and where it occurred
[out.(['max' enName]), maxInd] = max(sampEns);
out.maxScale = scaleRange(maxInd);
% Minimum, and where it occurred
[out.(['min' enName]), minInd] = min(sampEns);
out.minScale = scaleRange(minInd);
% Mean, std, coefficient of variation:
out.(['mean' enName]) = mean(sampEns);
out.(['std' enName]) = std(sampEns);
out.(['cv' enName]) = out.(['std' enName]) / out.(['mean' enName]);
% Mean change across the range of scales:
out.meanch = mean(diff(sampEns));

% Trend across scales: robust linear fit of the entropy vs. scale (a more
% principled alternative to meanch's raw adjacent-scale averaging -- uses
% all scales at once, and comes with its own standard error). A robust fit needs
% more points than its two coefficients:
if numScales > 2
	[linfit, linstats] = robustfit(scaleRange(:), sampEns);
	out.slope = linfit(2);
	out.slopeSE = linstats.se(2);
else
	out.slope = NaN;
	out.slopeSE = NaN;
end

end
