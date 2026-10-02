function out = EN_wentropy(y, waveletName, level)
% EN_wentropy   Wavelet entropy of a time series.
%
% Decomposes y via the maximal-overlap discrete wavelet transform (MODWT) into
% level detail scales plus the remaining smooth (scaling) band, i.e., level + 1
% bands, computes each band's share of the signal's total energy,
% p_j = E_j / sum(E), and returns the Shannon entropy of this relative-energy
% distribution across bands, normalized to [0,1] by its maximum possible value,
% log2(level + 1). Low values mean the energy is concentrated in few bands; high
% values that it is spread evenly across bands. Uses MATLAB's wentropy (Wavelet
% Toolbox).
%
% ---INPUTS:
% y, the input time series
% waveletName [optional], the wavelet used for the MODWT decomposition (default:
%    'sym4')
% level [optional], the number of decomposition levels (default: 5)
%
% ---OUTPUTS:
% a scalar: the normalized wavelet entropy (NaN if the decomposition fails, e.g.,
% for a series too short for the requested number of levels).
%
% ---REFERENCES:
% O. A. Rosso, S. Blanco, J. Yordanova, V. Kolev, A. Figliola, M. Schuermann,
% E. Basar, "Wavelet entropy: a new tool for analysis of short duration brain
% electrical signals", J. Neurosci. Methods 105(1) 65 (2001).
%
% ---NOTES:
% The output is invariant to rescaling y, and is bounded in [0,1] (the value 1 is
% reached when the energy is equal in all level + 1 bands). The normalizing
% constant is log2(level + 1) (wentropy with 'Scaled' true, its default); this was
% checked numerically against the entropy of the relative band energies.
% level is fixed by default (rather than left to wentropy's automatic choice,
% floor(log2(length(y)))) because the number of levels sets the normalizing
% denominator, so letting it grow with length(y) introduced a strong length
% dependence. With level fixed, the value for white noise is independent of length
% (about 0.75, the entropy of the energy shares 1/2, 1/4, 1/8, 1/16, 1/32, 1/32
% of white noise, divided by log2(6)); the value changes with length only for
% series whose energy depends on length, e.g., a random walk, whose energy
% concentrates in the smooth band as the series grows. level = 5 needs
% length(y) >= ~64 for a non-degenerate MODWT decomposition.
% The earlier implementation used the legacy wentropy(x,'shannon') cost-function
% syntax on the raw values, which is only a valid entropy for x pre-normalized so
% that sum(x.^2) = 1 and gave negative, scale-dependent values for z-scored input.

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
	level = 5; % fixed (see NOTES: avoids length dependence from an N-dependent level)
end

% ------------------------------------------------------------------------------
% Compute the (scaled, global) wavelet entropy
% ------------------------------------------------------------------------------
try
	ent = wentropy(y, 'Wavelet', waveletName, 'Distribution', 'global', 'Level', level, 'Scaled', true);
catch
	% Data-dependent (e.g., series too short for the requested/default level):
	out = NaN; return
end

out = ent;

end
