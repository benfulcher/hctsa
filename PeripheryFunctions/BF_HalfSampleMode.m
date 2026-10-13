function m = BF_HalfSampleMode(y)
% BF_HalfSampleMode   Half-sample mode: a robust, bin-free estimate of the mode.
%
% Repeatedly keeps the half of the (sorted) data, ceil(n/2) consecutive values,
% that spans the shortest interval, until three or fewer values remain; the mode
% is then the mean of the two closest of them (the middle one if equally close).
% It has no bins or bandwidth to choose, is closed-form, and changes continuously
% with the data except where two windows have exactly the same width (the first,
% lowest, is taken).
%
% ---INPUTS:
% y, the data vector (NaNs are ignored)
%
% ---OUTPUTS:
% m, the estimated mode.
%
% ---REFERENCES:
% D. R. Bickel and R. Fruhwirth, "On a fast, robust estimator of the mode:
% Comparisons to other robust estimators with applications", Comput. Stat. Data
% Anal. 50(12), 3500 (2006).

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

y = sort(y(~isnan(y(:))));
while length(y) > 3
	h = ceil(length(y) / 2); % the number of values in a half-sample
	[~, i] = min(y(h:end) - y(1:end - h + 1)); % the shortest window of h consecutive values
	y = y(i:i + h - 1);
end

switch length(y)
	case 3
		if y(2) - y(1) < y(3) - y(2)
			m = mean(y(1:2));
		elseif y(2) - y(1) > y(3) - y(2)
			m = mean(y(2:3));
		else
			m = y(2);
		end
	otherwise % 1 or 2 values remain
		m = mean(y);
end

end
