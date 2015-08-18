function groupIndices = TS_LabelGroups(whatData,keywordGroups,saveBack)
% TS_LabelGroups    Label groups of a time series using assigned keywords
%
% You provide a set of keyword options to store a specific grouping of time series.
% Useful when doing a classification task -- can store your classifications
% in the local structure arrays.
%
%---EXAMPLE USAGE:
%
% Label all time series with 'disease' in one group and all with 'healthy' in
% another group:
%
% groupIndices = TS_LabelGroups('norm',{'disease',0;'healthy',0},'ts');
%
%---INPUTS:
% whatData: Where to retrive from (and write back to): 'orig', 'norm', or 'cl'
%
% keywordGroups: The keyword groups, a cell of strings, as
%                    {'Keyword_1',numberToRetrive; 'keyword2',numberToRetrive,...}
%                   Use '0' to retrieve all of a given class.
%                   Can also use an empty label, '', to select anything at
%                   random from all time series.
%
% TsorOps: Whether grouping is for operations ('ops') or time series ('ts')
%
% saveBack: Can set to 0 to stop saving the grouping back to the input file.
%
%---OUTPUTS:
% groupIndices: the indicies corresponding to each keyword in keywordGroups.

% --------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
    fprintf(1,'Retrieving data from HCTSA_N by default.\n');
end
if ~isstruct(whatData) && ~ismember(whatData,{'orig','norm'})
    error('When specifying data, we need ''orig'', ''norm'', or a custom filename.')
end

if nargin < 2
    keywordGroups = '';
    % Try to assign by unique keywords later
end
if ~isempty(keywordGroups) && ischar(keywordGroups);
    fprintf(1,'Grouping all items with ''%s''.\n',keywordGroups);
    keywordGroups = {keywordGroups,0};
end

if nargin < 4 || isempty(saveBack)
    saveBack = 1; % Saves the grouping back to the HCTSA_*.loc file
end

% ------------------------------------------------------------------------------
%% Load data from file
% ------------------------------------------------------------------------------
[~,TimeSeries,~,theFile] = TS_LoadData(whatData);
Keywords = {TimeSeries.Keywords};
IDs = [TimeSeries.ID];
numTimeSeries = length(TimeSeries);

% ------------------------------------------------------------------------------
% Set default keywords?
% ------------------------------------------------------------------------------
if isempty(keywordGroups)
    fprintf(1,'No keywords assigned for labeling. Attempting to use unique keywords from data...?\n');
    UKeywords = unique(Keywords);
    NumUniqueKeywords = length(UKeywords);
    fprintf(1,'Shall I use the following %u keywords?: %s\n',NumUniqueKeywords,BF_cat(UKeywords,',',''''));
    reply = input('[y] for ''yes''','s');
    if strcmp(reply,'y')
        keywordGroups = cell(NumUniqueKeywords,2);
        for i = 1:NumUniqueKeywords
            keywordGroups{i,1} = UKeywords{i};
            keywordGroups{i,2} = 0;
        end
    else
        fprintf(1,'Ok then, thanks anyway\n'); return
    end
end

% ------------------------------------------------------------------------------
%% Label groups from keywords
% ------------------------------------------------------------------------------
if ~all(cellfun(@ischar,keywordGroups(:))) % Hopefully, specified numbers of each keyword
    keywordNumbers = horzcat(keywordGroups{:,2}); % Just the number of each part
    keywordGroups = keywordGroups(:,1)'; % Just the keyword parts, a cell of strings
else
    keywordNumbers = zeros(length(keywordGroups),1); % include all of each keyword
end
numGroups = length(keywordGroups); % The number of groups
Keywords = SUB_cell2cellcell(Keywords);

timer = tic;
for jo = 1:numGroups
    if ~isempty(keywordGroups{jo}) % Collect time series with this keyword
        groupIndices{jo} = find(cellfun(@(x)any(ismember(keywordGroups{jo},x)),Keywords));
        if isempty(groupIndices{jo})
            fprintf(1,'No matches found for ''%s''.\n',keywordGroups{jo});
        end
        if (keywordNumbers(jo) ~= 0) && (keywordNumbers(jo) < length(groupIndices{jo})) % Take a random subset of matches
            rperm = randperm(length(groupIndices{jo}));
            groupIndices{jo} = groupIndices{jo}(rperm(1:keywordNumbers(jo)));
        end
    else % Take a certain number of random time series
         % integer: retrieve this many: in randomorder
        rperm = randperm(length(Keywords));
        keywordGroups{jo} = 'Others';
        groupIndices{jo} = [];
        notKeywordGroups = keywordGroups(setxor(1:numGroups,jo));
        for i = 1:length(Keywords)
            if all(~ismember(notKeywordGroups,Keywords{rperm(i)}))
                groupIndices{jo} = [groupIndices{jo}; rperm(i)];
                if (length(groupIndices{jo}) == keywordNumbers(jo))
                    break
                end
            end
        end
    end
end
fprintf(1,'Group labeling complete in %s.\n',BF_thetime(toc(timer)));
clear timer % stop timing

% More feedback
fprintf(1,'We found:\n');
for i = 1:numGroups
    fprintf(1,'%s -- %u matches (/%u)\n',keywordGroups{i},length(groupIndices{i}),numTimeSeries);
end

% ------------------------------------------------------------------------------
%% Save back to the input file?
% ------------------------------------------------------------------------------
if saveBack
    % You don't need to check variables, you can just append back to the input file:
    if ~all(cellfun(@isempty,groupIndices))
        fprintf(1,'Saving group labels and information back to %s...',theFile);

        % First append/overwrite group names
        groupNames = keywordGroups;

        % Then overwrite labels
        theGroups = BF_ToGroup(groupIndices,numTimeSeries)';
        % Now we need to make the cells
        theGroupsCell = cell(size(theGroups));
        % Cannot find an in-built for this... :-/
        for i = 1:length(theGroups), theGroupsCell{i} = theGroups(i); end

        % First remove Group field if it exists
        if isfield(TimeSeries,'Group')
            TimeSeries = rmfield(TimeSeries,'Group');
        end

        % Add fields to the TimeSeries structure array
        newFieldNames = fieldnames(TimeSeries);
        % Add two new fields:
        newFieldNames{length(newFieldNames)+1} = 'Group';

        % Then append the new group information:
        % (some weird bug -- squeeze is sometimes needed here...:)
        TimeSeries = cell2struct([squeeze(struct2cell(TimeSeries));theGroupsCell],newFieldNames);
        % {'ID','Name','Keywords','Length','Data','Group'}

        % Save everything back to file:
        save(theFile,'TimeSeries','-append')
        save(theFile,'groupNames','-append')
        fprintf(1,' Saved.\n');
    end
end

end
