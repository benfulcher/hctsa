function [groupLabels,classLabels,groupLabelsInteger,numGroups] = ExtractGroupLabels(TimeSeries)

if ismember('Group',TimeSeries.Properties.VariableNames) && ~all(isundefined(TimeSeries.Group))
    groupLabels = TimeSeries.Group;
else
    % No group information assigned to time series
    groupLabels = cell(height(TimeSeries),1);
    groupLabels(:) = {'allData'};
    groupLabels = categorical(groupLabels);
end

classLabels = categories(groupLabels);
groupLabelsInteger = arrayfun(@(x)find(classLabels==x),groupLabels);
numGroups = length(classLabels);

end
