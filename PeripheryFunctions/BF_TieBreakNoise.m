function y = BF_TieBreakNoise(y, seed)
% BF_TieBreakNoise  Adds tiny, reproducible jitter to break exact ties in y.
%
% Returns y unchanged unless it has a high proportion of repeated values
% (< 90% unique), in which case Gaussian noise with standard deviation
% 1e-10 x std(y) is added: small enough to leave any well-behaved continuous
% series' statistics untouched, but enough to break the exact ties that make
% nearest-neighbor (e.g., Kraskov/KSG mutual information) estimators
% degenerate on quantized or periodic-orbit data.
%
% The noise is drawn from a private, fixed-seed RandStream (rather than the
% global one, or JIDT's own internal Java RNG via its NOISE_LEVEL_TO_ADD
% property, which cannot be seeded in the bundled JIDT build), so that the
% same input always yields the same jittered output, and nothing about the
% caller's global random state is consumed or disturbed.
%
%---INPUTS:
% y, the input vector (or matrix; the repeat test and noise scale use all elements)
% seed, [opt, default 0] seed for the private stream (use different seeds for
%       different variables that will be jittered and then compared to one another)
%
%---OUTPUTS:
% y, the (possibly) jittered input, same size.

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

if nargin < 2 || isempty(seed)
    seed = 0;
end
if isempty(y) || numel(y) < 2
    return
end
uniqueFrac = numel(unique(y(:))) / numel(y);
sigma = std(y(:));
if uniqueFrac < 0.9 && sigma > 0
    rs = RandStream('mt19937ar', 'Seed', seed); % private stream: reproducible, leaves the global stream alone
    y = y + 1e-10 * sigma * randn(rs, size(y));
end

end
