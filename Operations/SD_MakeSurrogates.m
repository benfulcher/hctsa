function out = SD_MakeSurrogates(x, surrMethod, numSurrs, extraParams, randomSeed)
% SD_MakeSurrogates   Generates surrogate time series.
%
% Method described relatively clearly in Guarin Lopez et al. (arXiv, 2010)
% Used bits of aaft code that references (and presumably was obtained from)
% "Surrogate data test for nonlinearity including monotonic
% transformations", D. Kugiumtzis, Phys. Rev. E, vol. 62, no. 1, 2000.
%
% Note that many other surrogate data methods exist that could later be
% implemented, cf. references in "Improvements to surrogate data methods for
% nonstationary time series", J. H. Lucio et al., Phys. Rev. E 85, 056202 (2012)
%
% ---INPUTS:
% x, the input time series
%
% surrMethod, the method for generating surrogates:
%             (i) 'RP' -- random phase surrogates
%             (ii) 'AAFT' -- amplitude adjusted Fourier transform
%             (iii) 'TFT' -- truncated Fourier transform
%             (iv) 'RandPerm' -- random permutation of the samples
%                    (destroys all temporal structure, not just nonlinear;
%                    equivalent to TSTOOL's surrogate method 3)
%
% numSurrs, the number of surrogates to generate
%
% extraParams, extra parameters required by the selected surrogate generation method
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed

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

beVocal = false; % Display text information/commentary to screen

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if size(x, 2) > size(x, 1)
	x = x'; % Time series must be a column vector -- SUB_RandomPhase's
			% Hermitian-symmetric phase construction (relied on by 'RP',
			% 'AAFT' and 'TFT' below) silently breaks for row-vector input
			% otherwise.
end

% Number of surrogates to generate:
if nargin < 3 || isempty(numSurrs)
	numSurrs = 1; % just create a single surrogate
end

% Any extra parameters (some methods require)
if nargin < 4
	extraParams = [];
end

% randomSeed: how to treat the randomization
if nargin < 5
	randomSeed = [];
end

% ------------------------------------------------------------------------------

N = length(x); % length of the time series
out = zeros(N, numSurrs); % each column is a new surrogate
tic % time it

% Control the random seed (for reproducibility):
BF_ResetSeed(randomSeed);

switch surrMethod
	case 'RP'
		% Random Phase Surrogates
		% Surrogates maintain linear correlations in the data, but any
		% nonlinear structure is destroyed by the phase randomization
		if beVocal
			fprintf(1, 'Constructing %u surrogates using the Random Phase Method\n', numSurrs)
			fprintf(1, ['Linear correlations are maintained but nonlinear structure will be ' ...
						'destroyed by the phase randomization\n'])
		end

		for surri = 1:numSurrs
			out(:, surri) = SUB_RandomPhase(x);
		end

	case 'AAFT'
		if beVocal
			fprintf(1, ['Constructing %u surrogates using the Amplitude Adjusted Fourier '...
						'Transform (AAFT) Method\n'], numSurrs)
			fprintf(1, ['Linear correlations are maintained but nonlinear structure will be destroyed ' ...
						'by the phase randomization. Amplitude Distribution is approximately maintained\n'])
		end

		% Sort and rank order the data
		[xSorted, ix] = sort(x);
		[~, xRO] = sort(ix); % rank ordered permutation

		for surri = 1:numSurrs
			% Rand order white Gaussian-distributed noise
			nSort = sort(randn(N, 1));
			y = nSort(xRO); % sorted Guassian white noise reordered as x

			% Random-phase surrogate of y (phase-randomized version of
			% random noise rank-ordered as x):
			yRP = SUB_RandomPhase(y);

			% Rank order x with respect to yRP:
			[~, ixyRP] = sort(yRP);
			[~, yRO] = sort(ixyRP);
			out(:, surri) = xSorted(yRO);
		end

	case 'TFT'
		if beVocal
			fprintf(1, ['Constructing %u surrogates using the Truncated Fourier '...
						'Transform (TFT) Method.\n'], numSurrs)
			fprintf(1, ['Low Frequency phases are preserved, and high frequency phases will be ' ...
						'randomized. A way of dealing with non-stationarity.\n'])
		end
		if isempty(extraParams)
			fprintf(1, 'You haven''t specified a cut-off frequency!! Setting N/8\n')
			fc = round(N / 8);
		else
			fc = extraParams; % extra input is the frequency cut-off
			if fc < 1
				fc = N * fc;
			end
		end

		for surri = 1:numSurrs
			out(:, surri) = SUB_RandomPhase(x, fc);
		end

	case 'RandPerm'
		% Random permutation surrogates: samples are shuffled, destroying
		% all temporal structure (not just nonlinear structure, unlike
		% RP/AAFT/TFT which preserve the linear/amplitude properties).
		if beVocal
			fprintf(1, 'Constructing %u surrogates using random permutation\n', numSurrs)
		end

		for surri = 1:numSurrs
			out(:, surri) = x(randperm(N));
		end

	otherwise
		error('Unknown surrogate generation method ''%s''', surrMethod)
end

% Cute farewell message
if beVocal
	fprintf(1, 'Generated %u %s surrogates in %s.\n', numSurrs, surrMethod, BF_TheTime(toc, 1))
end

end

% ------------------------------------------------------------------------------
function xNew = SUB_RandomPhase(x, fc)
	% Random-phase-Fourier-transform surrogate of column vector x.
	%
	% Preserves the magnitude spectrum exactly (reusing abs(fft(x)) as-is,
	% rather than manually re-assembling it from a half-spectrum slice --
	% that reassembly is what broke for row-vector input in an earlier
	% version of this code, since flipud() silently no-ops on a row
	% vector instead of reversing it) and randomizes every phase except
	% the DC and (for even N) Nyquist bins, which are preserved rather
	% than zeroed: for a real signal both are necessarily exactly 0 or pi,
	% and forcing them to 0 (as an earlier version did for the DC bin)
	% flips the sign of the surrogate's mean whenever x has a negative
	% mean. Randomized phases are mirrored (negated) onto the
	% conjugate-symmetric half of the spectrum so the result is real by
	% construction for any N, even or odd -- no truncation of x is needed
	% either way (an earlier version dropped the last sample for odd N,
	% then mismatched the FFT/IFFT lengths, distorting the magnitude
	% spectrum it was supposed to preserve).
	%
	% If fc (a bin count) is given, phases up to and including bin fc are
	% also preserved rather than randomized -- a truncated Fourier
	% transform (TFT) surrogate, for dealing with non-stationarity by
	% only randomizing high-frequency phases.
	if nargin < 2 || isempty(fc)
		fc = 0;
	end

	N = length(x);
	z = fft(x);
	zMag = abs(z);
	zPhase = angle(z);

	if mod(N, 2) == 0
		nFree = N / 2 - 1; % bins 2..N/2 (excludes DC at 1 and Nyquist at N/2+1)
	else
		nFree = (N - 1) / 2; % bins 2..(N+1)/2 (no separate Nyquist bin when N is odd)
	end
	freeBins = (2:nFree + 1)';

	randPhase = 2 * pi * rand(nFree, 1);
	nKeep = min(fc, nFree);
	randPhase(1:nKeep) = zPhase(freeBins(1:nKeep)); % preserve low-frequency phases (TFT)

	if mod(N, 2) == 0
		newPhase = [zPhase(1); randPhase; zPhase(N / 2 + 1); -flipud(randPhase)];
	else
		newPhase = [zPhase(1); randPhase; -flipud(randPhase)];
	end

	zNew = zMag .* exp(newPhase * 1i);
	xNew = real(ifft(zNew));
end
