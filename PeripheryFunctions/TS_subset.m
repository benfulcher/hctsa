function [TS_DataMat,TimeSeries,Operations] = TS_subset(whatData,ts_ids_keep,op_ids_keep,doSave)
% TS_subset save a given subset of data, based on time series and operation IDs

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
if nargin < 2
    ts_ids_keep = []; % all
end
if nargin < 3
    op_ids_keep = []; % all
end
if nargin < 4
    doSave = 1;
end

if isempty(ts_ids_keep) && isempty(op_ids_keep)
    error('Nothing to subset!');
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
numTimeSeries = length(TimeSeries);
numOperations = length(Operations);

%-------------------------------------------------------------------------------
% Do the subsetting:
%-------------------------------------------------------------------------------
i_keep = struct;
%--Match to TimeSeries IDs:
i_keep.TimeSeries = MatchMe([TimeSeries.ID],ts_ids_keep);
%--Match to Operation IDs:
i_keep.Operations = MatchMe([Operations.ID],op_ids_keep);
% Subset:
TS_DataMat = TS_DataMat(i_keep.TimeSeries,i_keep.Operations);
TimeSeries = TimeSeries(i_keep.TimeSeries);
Operations = Operations(i_keep.Operations);

fprintf('Data now contains %u -> %u time series and %u -> %u operations\n',...
            numTimeSeries,length(TimeSeries),numOperations,length(Operations))

if doSave
    % Save result to file

    % Remove group information because this will no longer be valid for sure
    if ~isempty(ts_ids_keep)
        rmfield(TimeSeries,'Group');
        fprintf('Warning: group information removed -- regenerate for subset data using TS_LabelGroups\n')
    end

    % Copy to a new subset .mat file, then save over with the new subset variables
    outputFileName = regexprep(whatDataFile,'.mat','_subset.mat');
    copyfile(whatDataFile,outputFileName);
    save(outputFileName,'TS_DataMat','TimeSeries','Operations','-append');

    varNames = whos('-file',whatDataFile);
    varNames = {varNames.name};
    if ismember('TS_Quality',varNames)
        load(whatDataFile,'TS_Quality')
        TS_Quality = TS_Quality(i_keep.TimeSeries,i_keep.Operations);
        save(outputFileName,'TS_Quality','-append');
    end
    if ismember('TS_CalcTime',varNames)
        load(whatDataFile,'TS_CalcTime')
        TS_CalcTime = TS_CalcTime(i_keep.TimeSeries,i_keep.Operations);
        save(outputFileName,'TS_CalcTime','-append');
    end

    % Possible inconsistency with grouping (e.g., if you remove a whole group of time series)
    % Matlab doesn't allow you to remove a variable from a .mat file easily
    if ~isempty(ts_ids_keep) && ismember('groupNames',varNames)
        groupNames = {};
        save(outputFileName,'groupNames','-append');
    end

    fprintf(1,'Data saved to %s!\n',outputFileName)
end

%-------------------------------------------------------------------------------
function ind = MatchMe(idsAll,idsMatch)
    % Find matches:
    if isempty(idsMatch)
        ind = logical(ones(size(idsAll)));
    else
        ind = ismember(idsAll,idsMatch);
    end
end

end
