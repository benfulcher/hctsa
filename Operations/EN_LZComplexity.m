function out = EN_LZComplexity(y, n, preProc)
% EN_LZComplexity   Lempel-Ziv complexity of an n-symbol encoding of a time series.
%
% The time series is coarse-grained into n symbols using equiprobable bins (each
% symbol is used equally often), and the normalized Lempel-Ziv complexity of the
% resulting symbol string is computed: the number c of distinct symbol sequences
% found when reading the string, divided by the number expected for a random
% (noise) string, c*log(L)/(L*log(n)) for a string of length L. A value near 1
% indicates a random sequence; lower values indicate structure that can be
% compressed.
%
% ---INPUTS:
% y, the input time series
% n, the (integer) number of symbols to encode the data into (default: 2, a
%    binary encoding)
% preProc [optional], first apply a given preprocessing to the time series. For
%    now, just 'diff' is implemented, which z-scores the incremental differences
%    and then applies the complexity method. An empty input applies no
%    preprocessing (default).
%
% ---OUTPUTS:
% a scalar: the normalized Lempel-Ziv complexity.
%
% ---REFERENCES:
% M. Small, "Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance", World Scientific, Nonlinear Science Series A, Vol. 52
% (2005).
%
% ---NOTES:
% Uses Michael Small's code 'complexity' (renamed MS_complexity here), available
% at http://small.eie.polyu.edu.hk/matlab/. The code is a wrapper for Michael
% Small's original code and uses the associated mex file compiled from
% complexitybs.c (renamed MS_complexitybs.c here).

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

if nargin < 2 || isempty(n)
	n = 2; % n-bit encoding
end
if nargin < 3
	preProc = []; % no preprocessing
end

% Apply some pre-processing to the time series before performing the analysis
if ischar(preProc)
	switch preProc
		case 'diff'
			y = zscore(diff(y));
		otherwise
			error('Unknown preprocessing setting ''%s''', preProc);
	end
end

% Run Michael Small's (mexed) code for calcaulting the Lempel-Ziv complexity:
out = MS_complexity(y, n);

end
