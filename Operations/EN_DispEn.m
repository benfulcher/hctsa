function out = EN_DispEn(y, m, c, tau, mappingHow)
% EN_DispEn     Dispersion entropy of a time series.
%
% Maps the time series onto c amplitude classes, symbolizes each embedding
% vector by the sequence of classes it visits (a 'dispersion pattern'), and
% returns the Shannon entropy of the resulting pattern distribution.
%
% Unlike permutation entropy (EN_PermEn), which records only the rank
% ordering within each embedding vector and so discards amplitude
% information entirely ([1,2,3] and [1,2,300] are the same pattern),
% dispersion entropy assigns each point to an amplitude class first, so the
% size of an excursion, not just its direction, shapes the symbol sequence.
% It is also markedly cheaper than sample entropy (EN_SampEn) and degrades
% more gracefully on short, noisy series.
%
% cf. M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for
% Time-Series Analysis", IEEE Signal Processing Letters 23(5) 610 (2016).
% DOI: 10.1109/LSP.2016.2542881
%
% Also returns the fluctuation-based variant, which symbolizes the
% *differences* between successive classes rather than the classes
% themselves, and so responds to the size of class-to-class changes rather
% than to absolute amplitude level (i.e., it is blind to a local trend that
% shifts every point into a higher class together):
%
% cf. H. Azami and J. Escudero, "Amplitude- and Fluctuation-Based
% Dispersion Entropy", Entropy 20(3) 210 (2018). DOI: 10.3390/e20030210
%
% ---INPUTS:
% y, the input time series
%
% m, the embedding dimension (default 2, following the source papers; the
%    number of possible patterns grows as c^m, so m must stay small for the
%    pattern frequencies to be estimable -- see the reliability check below)
%
% c, the number of amplitude classes (default 6, the value used throughout
%    the source papers; c > 1 is required, since c = 1 puts every point in
%    the same class and yields the trivial single-pattern case)
%
% tau, the time delay (default 1; can also be 'ac' or 'mi', resolved as in
%      the rest of the library)
%
% mappingHow, how to map the time series onto (0,1) before classifying:
%       (i) 'ncdf' (default), the normal cumulative distribution function
%           with the series' own mean and standard deviation. This is the
%           mapping the method was introduced with: a linear mapping assigns
%           the majority of points to only a few classes whenever the
%           maximum or minimum is far from the median, so a single outlier
%           can collapse the symbolization.
%       (ii) 'linear', a min-max rescaling onto [0,1]. Retained mainly
%            because the worked example in the source papers uses it (and
%            this operation's unit test reproduces that example), but it is
%            outlier-sensitive for the reason above.
%
% ---OUTPUTS:
% dispEn, normDispEn: the dispersion entropy (in nats) and the same
%       normalized by its maximum possible value, log(c^m).
% fDispEn, normFDispEn: the fluctuation-based dispersion entropy and the
%       same normalized by log((2c-1)^(m-1)).
%
% Only the normalized outputs are registered as hctsa features: for a fixed
% (m,c) the raw and normalized versions differ by a constant factor and are
% therefore the same feature.

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
%% Check inputs and set defaults
% ------------------------------------------------------------------------------
y = y(:);
N = length(y);

if nargin < 2 || isempty(m)
	m = 2;
end
if nargin < 3 || isempty(c)
	c = 6;
end
if nargin < 4 || isempty(tau)
	tau = 1;
end
if nargin < 5 || isempty(mappingHow)
	mappingHow = 'ncdf';
end

if c < 2
	error('Need at least two amplitude classes (c = %g given)', c);
end
if m < 1
	error('Embedding dimension must be at least 1 (m = %g given)', m);
end

% Resolve a character-specified time delay, as elsewhere in the library:
if ischar(tau)
	switch tau
		case 'ac'
			tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
		case 'mi'
			tau = CO_FirstMin(y, 'mi');
		otherwise
			error('Unknown time delay ''%s''', tau);
	end
	if isnan(tau)
		out = NaN; return % data-dependent: no correlation length could be estimated
	end
end

numVectors = N - (m - 1) * tau; % number of embedding vectors
if numVectors < 5
	% Too few embedding vectors to estimate any pattern distribution
	% (matching the floor EN_PermEn applies for the same reason)
	warning('Time series (N = %u) too short for dispersion entropy at m = %u, tau = %u', N, m, tau);
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% Map the series onto (0,1), then onto c amplitude classes
% ------------------------------------------------------------------------------
switch mappingHow
	case 'ncdf'
		sigma = std(y);
		if sigma == 0
			% Constant series: the NCDF is degenerate (every point identical)
			warning('Constant time series has no dispersion structure');
			out = NaN; return
		end
		yMapped = normcdf(y, mean(y), sigma);
	case 'linear'
		yRange = max(y) - min(y);
		if yRange == 0
			warning('Constant time series has no dispersion structure');
			out = NaN; return
		end
		yMapped = (y - min(y)) / yRange;
	otherwise
		error('Unknown mapping ''%s'' (expected ''ncdf'' or ''linear'')', mappingHow);
end

% Assign to integer classes 1:c (Rostaghi & Azami's z = round(c*y + 0.5)).
% The rounding can reach c+1 at the very top of the range (exactly y = 1,
% which 'linear' always attains and 'ncdf' attains whenever normcdf
% saturates on an extreme value), and 0 is unreachable but guarded for
% symmetry, so clamp into 1:c:
z = round(c * yMapped + 0.5);
z = min(max(z, 1), c);

% ------------------------------------------------------------------------------
%% Form the embedding vectors of class indices
% ------------------------------------------------------------------------------
% Z(i,k) = z(i + (k-1)*tau), one row per embedding vector:
Z = z((1:numVectors)' + (0:m - 1) * tau);

% ------------------------------------------------------------------------------
%% Dispersion entropy: patterns are the class sequences themselves (c^m)
% ------------------------------------------------------------------------------
% Encode each row as a base-c integer in 1:c^m:
placeValues = c.^(m - 1:-1:0)';
patternIdx = (Z - 1) * placeValues + 1;
counts = accumarray(patternIdx, 1, [c^m, 1]);
p = counts / numVectors;
p = p(p > 0); % 0*log(0) = 0

out.dispEn = -sum(p .* log(p));
out.normDispEn = out.dispEn / log(c^m);

% ------------------------------------------------------------------------------
%% Fluctuation-based dispersion entropy: patterns are the successive
%% class DIFFERENCES, each in -(c-1):(c-1), giving (2c-1)^(m-1) patterns
% ------------------------------------------------------------------------------
if m < 2
	% A single-element vector has no differences to symbolize
	out.fDispEn = NaN;
	out.normFDispEn = NaN;
	return
end

dZ = diff(Z, 1, 2) + (c - 1); % shift -(c-1):(c-1) onto 0:(2c-2)
numFluctPatterns = (2 * c - 1)^(m - 1);
fPlaceValues = (2 * c - 1).^(m - 2:-1:0)';
fPatternIdx = dZ * fPlaceValues + 1;
fCounts = accumarray(fPatternIdx, 1, [numFluctPatterns, 1]);
pF = fCounts / numVectors;
pF = pF(pF > 0);

out.fDispEn = -sum(pF .* log(pF));
out.normFDispEn = out.fDispEn / log(numFluctPatterns);

end
