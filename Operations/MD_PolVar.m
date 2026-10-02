function p = MD_PolVar(x, d, D)
% MD_PolVar   The POLVARd measure of a time series.
%
% Measures the probability of obtaining a sequence of D consecutive ones or zeros
% in the series of absolute differences, binarized at a threshold d. Each step
% (the absolute difference between consecutive values) is coded 1 if it is at
% least d and 0 otherwise. Scanning along the code from the start, windows of D
% consecutive symbols that are all ones or all zeros are counted, and after a
% match the scan moves on by D steps (so matches do not overlap). The output is
% the number of matches divided by the number of steps.
%
% The first mention may be in Wessel et al., PRE (2000), called Plvar
% cf. "Short-term forecasting of life-threatening cardiac arrhythmias based on
% symbolic dynamics and finite-time growth rates",
%       N. Wessel et al., Phys. Rev. E 61(1) 733 (2000)
%
% Although the original measure used raw thresholds, d, on RR interval sequences
% (measured in milliseconds), this code can be applied to general z-scored time
% series. So now d is not the time difference in milliseconds, but in units of
% std.
%
% The measure was originally applied to sequences of RR intervals and this code
% is heavily derived from that provided by Max A. Little, January 2009.
% cf. http://www.maxlittle.net/
%
% ---INPUTS:
% x, the input time series
% d, the symbolic coding (amplitude) difference: steps of at least this size are
%    coded 1 (default: 1)
% D, the word length (default: 6, the classical value)
%
% ---OUTPUTS:
% p, a scalar: the probability (proportion of steps) of obtaining a sequence of D
%    consecutive ones or zeros, as a count of non-overlapping such sequences
%    divided by the number of differences
%
% ---REFERENCES:
% N. Wessel et al., "Short-term forecasting of life-threatening cardiac arrhythmias
% based on symbolic dynamics and finite-time growth rates", Phys. Rev. E 61(1)
% 733 (2000).

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Max A. Little, <max.a.little@gmail.com>,
% <http://www.maxlittle.net/> and Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs:
% ------------------------------------------------------------------------------

% coding time difference, d
if nargin < 2 || isempty(d)
	d = 1;
end

% Word length, D
if nargin < 3 || isempty(D)
	D = 6;   % Look for symbol words of this default length
end

% ------------------------------------------------------------------------------

dx = abs(diff(x)); % absolute differences in consecutive values of the time series
N = length(dx); % number of differences in the input time series

% binary representation of time series based on consecutive changes being
% greater than d/1000...
xsym = (dx >= d); % consecutive differences exceed than some threshold, d
zseq = zeros(D, 1);
oseq = ones(D, 1);

% Search for D consecutive zeros/ones
i  = 1;
pc = 0;

% seqcnt = 0;
while (i <= (N - D + 1))
	xseq = xsym(i:(i + D - 1));
	%    seqcnt = seqcnt + 1;
	if (sum(xseq == zseq) == D) || (sum(xseq == oseq) == D)
		pc = pc + 1;
		i  = i + D;
	else
		i = i + 1;
	end
end

p = pc / N;
% p = pc / seqcnt;

end
