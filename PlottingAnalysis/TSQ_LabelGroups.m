% --------------------------------------------------------------------------
% TSQ_LabelGroups
% --------------------------------------------------------------------------
% 
% You provide a set of keyword options to store a specific grouping of time series.
% Useful when doing a classification task -- can store your classifications
% in the local structure arrays.
% 
% Requires a very specific (and unfortunately ~unintuitive) structure:
%  {'Keyword_1',NumberToRetrive;'Keyword2',NumberToRetrive,...}
%  Use '0' to retrieve all of a given class.
%  Can also use an empty label, '', to select anything at random from all time
%  series.
% 
%---EXAMPLE USAGE:
% groupIndices =
% TSQ_LabelGroups('norm',{'space',100;'',200;'medical',0;...},'ts');
% 
%---INPUTS:
% keywordGroups: The keyword groups, a cell of strings: 
% TsorOps: Whether grouping is for operations ('ops') or time series ('ts')
% whatData: Where to retrive from (and write back to): 'orig', 'norm', or 'cl'
% saveBack: Can set to 0 to stop saving the grouping back to the input file.
%
%---OUTPUTS:
% groupIndices: the indicies corresponding to each keyword in keywordGroups
% 
%---HISTORY:
% Previously named 'SUB_autolabelQ'
% 
% --------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function groupIndices = TSQ_LabelGroups(whatData,keywordGroups,TsorOps,saveBack)

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
    fprintf(1,'Retrieving data from HCTSA_N by default.\n');
end
if ~isstruct(whatData) && ~ismember(whatData,{'orig','norm','cl'})
    error('When specifying data, we need ''orig'', ''norm'', or ''cl''.')
end

if nargin < 2
    keywordGroups = '';
    % Try to assign by unique keywords later
end
if ~isempty(keywordGroups) && ischar(keywordGroups);
    fprintf(1,'Grouping all items with ''%s''.\n',keywordGroups);
    keywordGroups = {keywordGroups,0};
end

if nargin < 3 || isempty(TsorOps)
    TsorOps = 'ts';
    fprintf(1,'Grouping time series.\n');
end
if ~ismember(TsorOps,{'ops','ts'})
    error('Specify either ''ops'' or ''ts''.')
end

if nargin < 4 || isempty(saveBack)
    saveBack = 1; % Saves the grouping back to the HCTSA_*.loc file
end

% ------------------------------------------------------------------------------
%% Load data from file
% ------------------------------------------------------------------------------
if isstruct(whatData)
    % Can make whatData a structure...? Some old functionality...//
    Keywords = whatData.Keywords;
    idsO = whatData.idsO;
else
    switch whatData
        case 'orig'
            TheFile = 'HCTSA_loc.mat';
        case 'norm'
            TheFile = 'HCTSA_N.mat';
        case 'cl'
            TheFile = 'HCTSA_cl.mat';
    end
    fprintf(1,'Loaded data from %s...',TheFile);
    if strcmp(TsorOps,'ts') % Time series
        load(TheFile,'TimeSeries')
        Keywords = {TimeSeries.Keywords};
        IDs = [TimeSeries.ID];
    else % Operations
        load(TheFile,'Operations')
        Keywords = {Operations.Keywords};
        IDs = [Operations.ID];
    end
    fprintf(1,' Loaded.\n');
end

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
    fprintf(1,'%s -- %u matches\n',keywordGroups{i},length(groupIndices{i}));
end

% ------------------------------------------------------------------------------
%% Save back to file?
% ------------------------------------------------------------------------------
if saveBack
    % You don't need to check variables, you can just append back to the input file:
    if ~all(cellfun(@isempty,groupIndices))
        fprintf(1,'Saving group labels and information back to %s...',TheFile);
        
        % First append/overwrite group names
        GroupNames = keywordGroups;
        
        % Then overwrite labels
        TheGroups = BF_ToGroup(groupIndices,length(TimeSeries))';
        % Now we need to make the cells
        TheGroupsCell = cell(size(TheGroups));
        % Cannot find an in-built for this... :-/
        for i = 1:length(TheGroups), TheGroupsCell{i} = TheGroups(i); end
        
        % First remove Group field if it exists
        if isfield(TimeSeries,'Group')
            TimeSeries = rmfield(TimeSeries,'Group');
        end
        
        % Add fields to the TimeSeries structure array
        NewFieldNames = fieldnames(TimeSeries);
        % Add two new fields:
        NewFieldNames{length(NewFieldNames)+1} = 'Group';

        % Then append the new group information:
        TimeSeries = cell2struct([struct2cell(TimeSeries);TheGroupsCell],NewFieldNames);
        % {'ID','FileName','Keywords','Length','Data','Group'}

        % Save everything back to file:
        save(TheFile,'TimeSeries','-append')
        save(TheFile,'GroupNames','-append')
        fprintf(1,' Saved.\n');
    end
end

end