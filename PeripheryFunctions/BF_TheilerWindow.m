function W = BF_TheilerWindow(y, spec, Nfrac)
% BF_TheilerWindow  Resolves a Theiler-window specification to a number of samples.
%
% Neighbor-based methods (recurrence, dimension, Lyapunov, and nonlinear
% prediction estimates, ...) exclude candidate neighbors that are close in
% time to a reference point, since those are close in state space only
% because successive values are correlated, not because the trajectory has
% returned (J. Theiler, Phys. Rev. A 34, 2427 (1986)). The window should
% therefore span the time over which values remain correlated, which differs
% from series to series: a fixed number of samples is too short for slowly
% varying series (and needlessly wasteful for fast ones), and a proportion of
% the series length also changes with the length of the series.
%
%---INPUTS:
% y, the time series (used to compute its autocorrelation time)
% spec, the Theiler window:
%       {'ac', k}: k times the first zero-crossing of the autocorrelation
%                  function (recommended);
%       an integer >= 0: a fixed number of samples;
%       a number in (0,1): a proportion of Nfrac (legacy; this scales with the
%                  length of the series).
% Nfrac, [opt, default numel(y)] the length that a proportional window refers to
%
%---OUTPUTS:
% W, the Theiler window (samples): neighbors j of a reference point i are
%       excluded when |i - j| <= W. NaN if the autocorrelation function never
%       crosses zero.

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

if nargin < 3 || isempty(Nfrac)
	Nfrac = numel(y);
end

if iscell(spec)
	if numel(spec) ~= 2 || ~strcmp(spec{1}, 'ac') || ~isnumeric(spec{2}) || spec{2} < 0
		error('Theiler window must be specified as {''ac'', k}, with k >= 0');
	end
	W = ceil(spec{2} * CO_FirstCrossing(y, 'ac', 0, 'discrete'));
elseif isnumeric(spec) && isscalar(spec) && spec >= 0
	if spec > 0 && spec < 1 % a proportion of the series length
		W = round(spec * Nfrac);
	else
		W = round(spec);
	end
else
	error('Unrecognized Theiler window specification');
end

end
