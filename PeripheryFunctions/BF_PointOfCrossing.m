function [firstCrossing, pointOfCrossing] = BF_PointOfCrossing(x,threshold)
% BF_PointOfCrossing  Linearly interpolate to the point of crossing a threshold
%
%---INPUT:
% x, a vector
% threshold, a threshold x crosses.
%
%---OUTPUTS:
% firstCrossing, the first discrete value after which a crossing event has occurred
% pointOfCrossing, the (linearly) interpolated point of crossing

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% Find index of x at which the first crossing event occurs:
if x(1) > threshold
    firstCrossing = find((x - threshold < 0),1,'first');
else
    firstCrossing = find((x - threshold > 0),1,'first');
end

if isempty(firstCrossing)
    % Never crosses
    N = length(x);
    firstCrossing = N;
    pointOfCrossing = N;
else
    % Continuous version---the point of crossing
    valueBeforeCrossing = x(firstCrossing - 1);
    valueAfterCrossing = x(firstCrossing);
    pointOfCrossing = firstCrossing - 1 + (threshold - valueBeforeCrossing)/(valueAfterCrossing - valueBeforeCrossing);
end

end
