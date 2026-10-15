function [timeSeriesGroup,classLabels,groupLabelsInteger,numGroups] = TS_ExtractGroupLabels(TimeSeries)
% TS_ExtractGroupLabels     Extract group labels from a TimeSeries table
%---INPUT:
% TimeSeries, the input TimeSeries table in HCTSA format.
%
%---OUTPUT:
% timeSeriesGroup, a categorical vector of group labels, one per time series
%                   (categories with no members are removed).
% classLabels, a cell of the names of the classes that have at least one member.
% groupLabelsInteger, a vector of integer class indices (indexing classLabels),
%                   with 0 for any time series that has no label.
% numGroups, the number of classes.
%
% If the time series have not been labeled (no Group variable, or every label is
% undefined), all time series are assigned to a single class, 'allData'.
% So a single class called 'allData' means no labels were assigned.

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

if ismember('Group',TimeSeries.Properties.VariableNames) && ~all(isundefined(TimeSeries.Group))
    timeSeriesGroup = TimeSeries.Group;
else
    % No group information assigned to time series
    timeSeriesGroup = cell(height(TimeSeries),1);
    timeSeriesGroup(:) = {'allData'};
    timeSeriesGroup = categorical(timeSeriesGroup);
end

% Only keep classes that have members, so that every caller sees occupied classes:
timeSeriesGroup = removecats(timeSeriesGroup);
classLabels = categories(timeSeriesGroup);
groupLabelsInteger = double(timeSeriesGroup); % NaN for undefined
groupLabelsInteger(isnan(groupLabelsInteger)) = 0;
numGroups = length(classLabels);

end
