function TellMeAboutLabeling(whatData)
% Tells the user information about the class-labeling of their data

%------------------------------------------------------------------------------
% If you use this code for your research, please cite these papers:
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
%------------------------------------------------------------------------------

% Can input TimeSeries table, or reference to where it is stored.
if istable(whatData)
    TimeSeries = whatData;
else
    [~,TimeSeries] = TS_LoadData(whatData);
end

% Check that the data are labeled:
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end

classLabels = categories(TimeSeries.Group);
numClasses = length(classLabels);
fprintf(1,'%u classes assigned to the time series in this dataset:\n',numClasses);
for i = 1:numClasses
    fprintf(1,'Class %u: %s (%u time series)\n',i,classLabels{i},...
            sum(TimeSeries.Group==classLabels{i}));
end

end
