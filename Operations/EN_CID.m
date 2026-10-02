function out = EN_CID(y)
% EN_CID   Complexity estimate of a time series from the length of its line graph.
%
% Estimates the 'complexity' of a time series as the stretched-out length of the
% line obtained by plotting the series as a line graph (the complexity estimate
% used in the complexity-invariant distance, CID). Two versions are computed:
% CE1, the root-mean-square of the increments, and CE2, the mean length of the
% line segments between successive points (unit time step, by Pythagoras). Each
% is also computed on the sorted time series, which gives the minimum value
% possible for these data, and the series' value is expressed as a ratio to it.
% Sums in the original definitions are replaced by means so that values scale
% properly with series length.
%
% ---INPUTS:
% y, the input time series
%
% ---OUTPUTS:
% A structure with fields:
% CE1, sqrt(mean(diff(y).^2))
% CE2, mean(sqrt(1 + diff(y).^2))
% minCE1, CE1 of the sorted time series
% minCE2, CE2 of the sorted time series
% CE1_norm, CE1/minCE1
% CE2_norm, CE2/minCE2
%
% ---REFERENCES:
% G. E. A. P. A. Batista, E. J. Keogh, O. M. Tataw, V. M. A. de Souza,
% "CID: an efficient complexity-invariant distance for time series",
% Data Min. Knowl. Disc. 28, 634-669 (2014). https://doi.org/10.1007/s10618-013-0312-3

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

% Original definition (in Table 2 of paper cited above)
% sum -> mean to deal with non-equal time-series lengths
% (now scales properly with length)
f_CE1 = @(x) sqrt(mean(diff(x).^2));

out.CE1 = f_CE1(y);

% Definition corresponding to the line segment example in Fig. 9 of the paper
% cited above (using Pythagoras's theorum):
f_CE2 = @(x) mean(sqrt(1 + diff(x).^2));

out.CE2 = f_CE2(y);

% Defined as a proportion of the minimum such value possible for this time series,
% this would be attained from putting close values close; i.e., sorting the time
% series

out.minCE1 = f_CE1(sort(y));
out.minCE2 = f_CE2(sort(y));

out.CE1_norm = out.CE1 / out.minCE1;
out.CE2_norm = out.CE2 / out.minCE2;

end
