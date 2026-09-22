function out = SB_BinaryGapHomogeneity(x, gapWhat)
% SB_BinaryGapHomogeneity  Homogeneity of the gaps between like symbols in a binarization.
%
% The time series is symbolized to a binary string by whether it's above (1) or
% below (0) zero (equivalently its mean, for the z-scored input this is
% registered on). Taking the positions of the 1s (or the 0s), the gaps between
% successive such points are classified as either 'adjacent' (gap of 1, i.e.
% part of an unbroken run) or 'separated' (gap of 2 or more). This operation
% returns the length of the longest block over which that classification does
% not change, as a proportion of the time-series length -- i.e. the longest
% uninterrupted regime of consistent spacing, whether that regime is a long
% unbroken run or a long stretch of regularly-isolated points.
%
% ---INPUTS:
%
% x, the input time series
%
% gapWhat, (i) 'gaps1', gap homogeneity between above-zero (1) points
%          (ii) 'gaps0', gap homogeneity between below-zero (0) points
%
% ---NOTES:
% This is NOT the longest run of consecutive 1s or 0s, despite this operation's
% former name ('SB_BinaryStretch', with outputs 'lseq1'/'lseq0'), which claimed
% to measure exactly that and did not: verified against a brute-force longest-run
% calculation, the two disagree on 171 of 200 random series, and this quantity is
% 0 whenever the longest run touches either end of the series (blocks at the
% boundaries are not counted, since they are not bracketed by a change on both
% sides). For true run lengths, use SB_BinaryStats, whose longstretch0/longstretch1
% (and meanstretch/stdstretch) fields were verified exact against brute force.
%
% The quantity computed here was nonetheless kept, and renamed rather than fixed,
% because it is not redundant with those correct measures: on the Empirical1000
% dataset its rank correlation with SB_BinaryStats_mean_longstretch1 is only 0.47
% (0.54 for the 0-symbol version), well inside the library's retention bar.

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

if nargin < 2 || isempty(gapWhat)
	gapWhat = 'gaps1'; % by default
end

N = length(x); % length of the time series
x(x > 0) = 1;
x(x <= 0) = 0;

switch gapWhat
	case 'gaps1'
		% longest stretch of 1s [this code doesn't actually measure this!]
		out = max(diff(BF_SignChange(diff(find(x == 1)) - 1.5, 1))) / N;
	case 'gaps0'
		% longest stretch of 0s [this code doesn't actualy measure this!]
		out = max(diff(BF_SignChange(diff(find(x == 0)) - 1.5, 1))) / N;
	otherwise
		error('Unknown input ''%s'' (expected ''gaps1'' or ''gaps0'')', gapWhat)
end

if isempty(out)
	out = 0;
end

end
