function iszscored = BF_iszscored(x)
% BF_iszscored  Crude check for whether a data vector is z-scored.
%
% (~eps-close to being) z-scored.
% Used for displaying warning messages for functions that require z-scored inputs.
%
%---INPUT:
% x, the input time series (or any vector)
%
%---OUTPUT:
% iszscored, a logical with the verdict.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

% Give it a bit of numerical lee-way... Down in the 2e-14 region:
numericThreshold = 100*eps;

iszscored = ((abs(mean(x)) < numericThreshold) && (abs(std(x)-1) < numericThreshold));

end
