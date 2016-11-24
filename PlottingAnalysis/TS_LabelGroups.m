function [groupLabels,newFileName] = TS_LabelGroups(keywordGroups,whatData,saveBack,filterMissing)
% TS_LabelGroups    Label groups of a time series using assigned keywords
%
% You provide a set of keyword options to store a specific grouping of time series.
% Useful when doing a classification task -- can store your classifications
% in the local structure arrays.
%
% Requires that each time series is labeled uniquely to a single group and that
% all groups contain members.
% (if needed, can trim your HCTSA file to generate a dataset that can meet these
% criteria)
%
%---EXAMPLE USAGE:
%
% Label all time series with keyword 'disease' as one group and all with the
% keyword 'healthy' in another group -- saving the group information back to
% HCTSA_N.mat:
%
% groupIndices = TS_LabelGroups({'disease','healthy'});
%
%---INPUTS:
% keywordGroups: The keyword groups, a cell of strings, as
%                   {'keyword_1', 'keyword2',...}
%                   Can also use an empty label, '', to select unique keywords
%                   automatically from the dataset.
%
% whatData: Where to retrive from (and write back to): 'HCTSA.mat' (default),
%           cf. TS_LoadData.
%
% saveBack: Can set to 0 to stop saving the grouping back to the input file.
%
% filterMissing: Set to 1 to remove data that don't match any keywords.
%
%---OUTPUTS:
% groupIndices: the indicies corresponding to each keyword in keywordGroups.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1
    keywordGroups = {};
    % Try to assign by unique keywords later
end
if ~isempty(keywordGroups) && ischar(keywordGroups);
    fprintf(1,'Grouping all items with ''%s''.\n',keywordGroups);
    keywordGroups = {keywordGroups};
end

if nargin < 2 || isempty(whatData)
    whatData = 'raw';
    fprintf(1,'Retrieving data from HCTSA.mat by default.\n');
end

if nargin < 3 || isempty(saveBack)
    saveBack = 1; % Saves the grouping back to the HCTSA_*.loc file
end

if nargin < 4
    filterMissing = 0; % don't filter out data with missing labels -- error if not all data are labeled
end

% ------------------------------------------------------------------------------
%% Load data from file
% ------------------------------------------------------------------------------
[~,TimeSeries,~,theFile] = TS_LoadData(whatData);
Keywords = SUB_cell2cellcell({TimeSeries.Keywords}); % Split into sub-cells using comma delimiter
numTimeSeries = length(TimeSeries);

% ------------------------------------------------------------------------------
% Set default keywords?
% ------------------------------------------------------------------------------
% Set group labels as each unique keyword in the data. Works only in simple cases.
if isempty(keywordGroups)
    fprintf(1,'No keywords assigned for labeling. Attempting to use unique keywords from data...\n');
    keywordsAll = [Keywords{:}]; % every keyword used across the dataset
    UKeywords = unique(keywordsAll);
    numUniqueKeywords = length(UKeywords);
    fprintf(1,'Shall I use the following %u keywords: %s?\n',numUniqueKeywords,...
                BF_cat(UKeywords,',',''''));
    reply = input('[y] for ''yes''...','s');
    if strcmp(reply,'y')
        keywordGroups = UKeywords;
    else
        fprintf(1,'No wuckers, thanks anyway ^_^\n'); return
    end
end

% ------------------------------------------------------------------------------
%% Label groups from keywords
% ------------------------------------------------------------------------------
numGroups = length(keywordGroups); % The number of groups

timer = tic;
groupIndices = false(numTimeSeries,numGroups);
for jo = 1:numGroups
    groupIndices(:,jo) = cellfun(@(x)any(ismember(keywordGroups{jo},x)),Keywords);
    if all(groupIndices(:,jo)==0)
        fprintf(1,'No matches found for ''%s''.\n',keywordGroups{jo});
    end
end
fprintf(1,'Group labeling complete in %s.\n',BF_thetime(toc(timer)));
clear timer % stop timing

%-------------------------------------------------------------------------------
% Checks:
%-------------------------------------------------------------------------------

% Check each group has some members:
emptyGroups = (sum(groupIndices,1)==0);
if any(emptyGroups)
    error('%u keywords have no matches to any time series in %s: %s',...
            sum(emptyGroups),theFile,BF_cat(keywordGroups(emptyGroups),',',''''));
end

% Check overlaps:
overlapping = (sum(groupIndices,2)>1);
if any(overlapping)
    error('%u time series have multiple group assignments: %s',sum(overlapping),...
                    BF_cat({TimeSeries(overlapping).Name},','));
end

% Check unlabeled:
newFileName = theFile; % by default you save back to the same file
unlabeled = (sum(groupIndices,2)==0);
if any(unlabeled)
    if ~filterMissing
        reply = input(sprintf('ERROR: %u/%u time series remain unlabeled (press enter to see them)',...
                                    sum(unlabeled),length(unlabeled)));
        isUnlabeled = find(unlabeled);
        for i = 1:length(isUnlabeled)
            fprintf('[%u] %s (%s)\n',TimeSeries(isUnlabeled(i)).ID, ...
                    TimeSeries(isUnlabeled(i)).Name,TimeSeries(isUnlabeled(i)).Keywords);
        end
        error('Unable to provide a unique label to all time series (NB: Can set filterMissing input to deal with this)');
    else
        reply = input(sprintf('%u/%u time series were unlabeled and WILL BE REMOVED from the dataset (press enter to see them)',...
                                    sum(unlabeled),length(unlabeled)));
        if ischar(whatData) && strcmp(theFile,'HCTSA_N.mat')
            warning('If data are already normalized, feature distributions will change after removing data; and normalizations may no longer be valid')
        end
        isUnlabeled = find(unlabeled);
        for i = 1:length(isUnlabeled)
            fprintf('[%u] %s (%s)\n',TimeSeries(isUnlabeled(i)).ID, ...
                    TimeSeries(isUnlabeled(i)).Name,TimeSeries(isUnlabeled(i)).Keywords);
        end
        input('Look alright? (cntrl-C to cancel)')

        % Remove the unlabled data from the TimeSeries structure and the data matrix
        % Filter data and save to new file:
        newFileName = TS_FilterData(theFile,[TimeSeries(~unlabeled).ID]);
        % Label this dataset with the same groups:
        groupLabels = TS_LabelGroups(keywordGroups,newFileName,saveBack,0);
        return
    end
end

%-------------------------------------------------------------------------------
% Everything checks out so now we can make group labels:
groupLabels = zeros(1,numTimeSeries);
for i = 1:numGroups
    groupLabels(groupIndices(:,i)) = i;
end

%-------------------------------------------------------------------------------
% User feedback:
%-------------------------------------------------------------------------------
fprintf(1,'We found:\n');
for i = 1:numGroups
    fprintf(1,'%s -- %u matches (/%u)\n',keywordGroups{i},sum(groupIndices(:,i)),numTimeSeries);
end

% ------------------------------------------------------------------------------
%% Save back to the input file?
% ------------------------------------------------------------------------------
if saveBack
    % You don't need to check variables, you can just append back to the input file:
    fprintf(1,'Saving group labels and information back to %s...',theFile);

    % First append/overwrite group names
    groupNames = keywordGroups;

    % Make a cell version of group indices (to use cell2struct)
    theGroupsCell = cell(size(groupLabels));

    % Cannot find an in-built function for this... :-/
    for i = 1:length(groupLabels), theGroupsCell{i} = groupLabels(i); end

    % First remove Group field if it exists
    if isfield(TimeSeries,'Group')
        TimeSeries = rmfield(TimeSeries,'Group');
    end

    % Add new field to the TimeSeries structure array
    newFieldNames = fieldnames(TimeSeries);
    newFieldNames{length(newFieldNames)+1} = 'Group';

    % Then append the new group information:
    % (some weird bug -- squeeze is sometimes needed here...:)
    TimeSeries = cell2struct([squeeze(struct2cell(TimeSeries));theGroupsCell],newFieldNames);

    % Save everything back to file:
    save(theFile,'TimeSeries','groupNames','-append')
    fprintf(1,' Saved.\n');
end

end
