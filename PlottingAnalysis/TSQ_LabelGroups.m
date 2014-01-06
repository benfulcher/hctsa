% TSQ_LabelGroups
% 
% You provide a set of keyword options to store a grouping of time series in the
% store.
% 
% Requires a very specific structure:
% {'Keyword_1',NumberToRetrive;'Keyword2',NumberToRetrive,...}
% Use '0' to retrieve all of a given class.
% Can also use an empty label, '', to select anything at random from all time series.
% 
% Example usage:
% KeywordGroups = {'space',100;'',200;'medical',0;...};
% 
%---INPUTS:
%--KeywordGroups: The keyword groups, a cell of strings: 
%--TsorOps: Whether grouping is for operations ('ops') or time series ('ts')
%--WhatData: Where to retrive from (and write back to): 'orig', 'norm', or 'cl'
%--SaveBack: Can set to 0 to stop saving the grouping back to the input file.
%
%---OUTPUTS:
%-- GroupIndices: the indicies corresponding to each keyword in KeywordGroups
% 
% [[Previously named 'SUB_autolabelQ']]
% 
% --------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function GroupIndices = TSQ_LabelGroups(KeywordGroups,TsorOps,WhatData,SaveBack)

%% Check inputs
if nargin < 1 || isempty(KeywordGroups)
    error('You must specify labels');
end
if ischar(KeywordGroups);
    fprintf(1,'Grouping all items with ''%s''.\n',KeywordGroups);
    KeywordGroups = {KeywordGroups,0};
end

if nargin < 2 || isempty(TsorOps)
    TsorOps = 'ts';
    fprintf(1,'Grouping time series.\n');
end
if ~ismember(TsorOps,{'ops','ts'})
    error('Specify either ''ops'' or ''ts''.')
end

if nargin < 3 || isempty(WhatData)
    WhatData = 'cl';
    fprintf(1,'Retrieving from HCTSA_cl by default.\n');
end
if ~isstruct(WhatData) && ~ismember(WhatData,{'orig','norm','cl'})
    error('Specify ''orig'', ''norm'', or ''cl''.')
end

if nargin < 4 || isempty(SaveBack)
    SaveBack = 1; % Saves the grouping back to the HCTSA_*.loc file
end

%% Load data from file
if isstruct(WhatData)
    % Can make WhatData a structure...? Some old functionality...//
    Keywords = WhatData.Keywords;
    idsO = WhatData.idsO;
else
    switch WhatData
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

if ~all(cellfun(@ischar,KeywordGroups(:))) % Have specified numbers of each
    KeywordNumbers = horzcat(KeywordGroups{:,2}); % Just the number of each part
    KeywordGroups = KeywordGroups(:,1)'; % Just the keyword parts, a cell of strings
else
    KeywordNumbers = zeros(length(KeywordGroups),1); % include all of each keyword
end
NumGroups = length(KeywordGroups); % The number of groups
Keywords = SUB_cell2cellcell(Keywords);

timer = tic;
for jo = 1:NumGroups
    if ~isempty(KeywordGroups{jo}) % Collect time series with this keyword
        GroupIndices{jo} = find(cellfun(@(x)any(ismember(KeywordGroups{jo},x)),Keywords));
        if isempty(GroupIndices{jo})
            fprintf(1,'No matches found for ''%s''.\n',KeywordGroups{jo});
        end
        if (KeywordNumbers(jo) ~= 0) && (KeywordNumbers(jo) < length(GroupIndices{jo})) % Take a random subset of matches
            rperm = randperm(length(GroupIndices{jo}));
            GroupIndices{jo} = GroupIndices{jo}(rperm(1:KeywordNumbers(jo)));
        end
    else % Take a certain number of random time series
         % integer: retrieve this many: in randomorder
        rperm = randperm(length(Keywords));
        KeywordGroups{jo} = 'Others';
        GroupIndices{jo} = [];
        notKeywordGroups = KeywordGroups(setxor(1:NumGroups,jo));
        for i = 1:length(Keywords)
            if all(~ismember(notKeywordGroups,Keywords{rperm(i)}))
                GroupIndices{jo} = [GroupIndices{jo}; rperm(i)];
                if (length(GroupIndices{jo}) == KeywordNumbers(jo))
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
for i = 1:NumGroups
    fprintf(1,'%s -- %u matches\n',KeywordGroups{i},length(GroupIndices{i}));
end

% Save the grouping back to file?
% Eventually would be nice to add a group index to each element, e.g., as [TimeSeries.Group]
% But now we'll just save GroupIndices
if SaveBack
    % Actually, you don't need to check variables, you can just append
    % Elements = who(matfile(TheFile));
    % Elements = whos('-file',TheFile);
    % if ismember({Elements.name},'GroupIndices')
    if ~all(cellfun(@isempty,GroupIndices))
        fprintf(1,'Saving group labels and information back to %s...',TheFile);
        
        % First append/overwrite group names
        GroupNames = KeywordGroups;
        
        % Then overwrite labels
        TheGroups = BF_ToGroup(GroupIndices,length(TimeSeries))';
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