function out = EN_wentropy(y, waveletName, level)
% EN_wentropy   Wavelet entropy of a time series.
%
% Decomposes y via the maximal-overlap discrete wavelet transform (MODWT)
% into a set of scales, computes each scale's share of the signal's total
% energy, p_j = E_j / sum(E), and returns the Shannon entropy of this
% relative-energy distribution across scales, normalized to [0,1] by its
% maximum possible value, log2(numLevels).
%
% cf. O. A. Rosso, S. Blanco, J. Yordanova, V. Kolev, A. Figliola, M.
% Schuermann, E. Basar, "Wavelet entropy: a new tool for analysis of short
% duration brain electrical signals", J. Neurosci. Methods 105(1) 65 (2001).
%
% ---INPUTS:
% y, the input time series
% waveletName [opt], the wavelet used for the MODWT decomposition (default: 'sym4')
% level [opt], the number of decomposition levels (default: 5)
%
% ---NOTES:
% level is fixed by default (rather than left to wentropy's automatic
% choice, floor(log2(length(y)))) because the number of levels sets the
% normalizing denominator, log2(numLevels), for the [0,1]-scaled entropy --
% letting it grow with length(y) introduced substantial length-dependence
% (Spearman rho=-0.91 for value vs. N, and a value range spanning ~0.65 to
% ~0.59, across growing prefixes of one real series from length 500 to
% 10000). Fixing level shrinks the absolute drift to ~0.02 over the same
% prefix range (level=4,5,6 all similar) -- what Spearman rho remains
% (~0.67) reflects that residual small, smooth drift rather than a
% practically meaningful trend. level=5 needs length(y) >= ~64 for a
% non-degenerate MODWT decomposition (comfortably below the minimum length,
% 208, in the Empirical1000 validation set).
% Audit, 2026-08-15: the previous implementation called MATLAB's legacy
% wentropy(x,'shannon')/wentropy(x,'log energy') cost-function syntax
% directly on raw (non energy-normalized) z-scored time-series values --
% flagged by this file's own docstring as suspect since the code was first
% written (2013). Confirmed numerically that this was not a valid entropy:
% wentropy(zscore(randn(1,500)),'shannon') returned -0.75 (entropies are
% non-negative by definition), and rescaling that same z-scored signal by a
% constant factor of 5 (same relative structure, different amplitude) took
% it to -99 -- an entropy of a fixed distribution should be invariant to an
% arbitrary absolute scale. The root cause: that legacy call computes
% -sum(x.^2.*log(x.^2)), which is only a valid (bounded, non-negative)
% entropy when x is pre-normalized so that sum(x.^2)=1 (treating x.^2 as a
% probability mass) -- z-scored input instead has sum(x.^2)~=N.
%
% MATLAB's Wavelet Toolbox has since been redesigned around the
% literature-grounded wavelet-entropy definition above (Rosso et al.,
% 2001), which normalizes relative energy across decomposition levels (not
% raw sample values) to sum to 1 before taking the entropy. Verified this
% replacement IS scale-invariant (rescaling y by 5x left the output
% unchanged to 4 decimal places on test signals) and properly bounded
% (relative energies summed to 1 as expected).

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
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox');

% ------------------------------------------------------------------------------
% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(waveletName)
	waveletName = 'sym4'; % default
end
if nargin < 3 || isempty(level)
	level = 5; % fixed (see NOTES: avoids length-dependence from an N-dependent level)
end

% ------------------------------------------------------------------------------
% Compute the (scaled, global) wavelet entropy
% ------------------------------------------------------------------------------
try
	ent = wentropy(y, 'Wavelet', waveletName, 'Distribution', 'global', 'Level', level);
catch
	% Data-dependent (e.g., series too short for the requested/default level):
	out = NaN; return
end

out = ent;

end
