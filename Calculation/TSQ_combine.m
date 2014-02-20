% ------------------------------------------------------------------------------
% TSQ_combine
% ------------------------------------------------------------------------------
% 
% This function joins two HCTSA_loc.mat files.
% Any data matrices are combined, and the guides are updated to reflect the
% concatenation.
% Note that in the case of duplicates, the first file will have precedence.
% Takes a union of time series, and an intersection of operations.
% 
%---INPUTS:
% HCTSA_loc_1: the path to the first HCTSA_loc.mat file
% HCTSA_loc_2: the path to the second HCTSA_loc.mat file
% 
%---OUTPUTS:
% Writes a new, combined HCTSA_loc.mat
%
% ------------------------------------------------------------------------------
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

function TSQ_combine(HCTSA_loc_1,HCTSA_loc_2)

if nargin < 2
    error('Must provide paths for two HCTSA_loc files')
end

% Combine the local files
HCTSA_locs = {HCTSA_loc_1, HCTSA_loc_2};

% Check paths point to valid files
% FullPath = cell(2,1);
for i = 1:2
    path = exist(HCTSA_locs{i});
    if (path==0)
        error('Could not find %s',HCTSA_locs{i});
    % else
    %     FullPath{i} = path;
    end
end

% Load the two local files
fprintf(1,'Loading data...');
LoadedData = cell(2,1);
for i = 1:2
    LoadedData{i} = load(HCTSA_locs{i});
end
fprintf(1,' Loaded.\n');

% Give some information
for i = 1:2
    fprintf(1,'%u: The file, %s, contains information for %u time series and %u operations.\n', ...
                    i,HCTSA_locs{i},length(LoadedData{i}.TimeSeries),length(LoadedData{i}.Operations));
end

% 1. Union of time series
% As a basic concatenation, then remove any duplicates
TimeSeries = cell2struct([struct2cell(LoadedData{1}.TimeSeries), ...
                          struct2cell(LoadedData{2}.TimeSeries)], ...
                                {'ID','FileName','Keywords','Length','Data'});
% Check for time series duplicates
[uniquetsids, ix] = unique(vertcat(TimeSeries.ID));
DidTrim = 0;
if length(uniquetsids) < length(TimeSeries)
    fprintf(1,'Trimming %u duplicate time series to a toal of %u\n', ...
                    length(TimeSeries)-length(uniquetsids),length(uniquetsids));
    TimeSeries = TimeSeries(ix);
    DidTrim = 1;
else
    fprintf(1,'All time series were distinct, we now have a total of %u.\n',length(TimeSeries));
end

% 2. Intersection of Operations
% Take intersection of operation ids, and take information from first input
[allopids,keepopi_1,keepopi_2] = intersect(vertcat(LoadedData{1}.Operations.ID),vertcat(LoadedData{2}.Operations.ID));
% keepopi_1 = ismember(LoadedData{1}.Operations.ID,allopids); % Indices of operations kept in Loaded Data (1)
% keepopi_2 = ismember(LoadedData{2}.Operations.ID,allopids); % Indices of operations kept in Loaded Data (2)
% Data from first file goes in (should be identical to in the second file)
Operations = LoadedData{1}.Operations(keepopi_1);

fprintf(1,'Keeping the %u overlapping operations.\n',length(allopids));

% 3. MasterOperations
% Take intersection, like operations -- those that are in both
[allmopids,keepmopi_1] = intersect(vertcat(LoadedData{1}.MasterOperations.ID),vertcat(LoadedData{2}.MasterOperations.ID));
MasterOperations = LoadedData{1}.MasterOperations(keepmopi_1);

% 3. Data:
GotData = 0;
if isfield(LoadedData{1},'TS_DataMat') && isfield(LoadedData{2},'TS_DataMat')
    % Both contain data matrices
    GotData = 1;
    fprintf(1,'Combining data matrices...');
    TS_DataMat = [LoadedData{1}.TS_DataMat(:,keepopi_1);LoadedData{2}.TS_DataMat(:,keepopi_2)];
    if DidTrim, TS_DataMat = TS_DataMat(ix,:); end
    fprintf(1,' Done.\n');
end

% 4. Quality
GotQuality = 0;
if isfield(LoadedData{1},'TS_Quality') && isfield(LoadedData{2},'TS_Quality')
    GotQuality = 1;
    % Both contain Quality matrices
    fprintf(1,'Combining quality label matrices...');
    TS_Quality = [LoadedData{1}.TS_Quality(:,keepopi_1);LoadedData{2}.TS_Quality(:,keepopi_2)];
    if DidTrim, TS_Quality = TS_Quality(ix,:); end
    fprintf(1,' Done.\n');
end

% 5. Calculation times
GotCalcTimes = 0;
if isfield(LoadedData{1},'TS_CalcTime') && isfield(LoadedData{2},'TS_CalcTime')
    GotCalcTimes = 1;
    % Both contain Calculation time matrices
    fprintf(1,'Combining calculation time matrices...');
    TS_CalcTime = [LoadedData{1}.TS_CalcTime(:,keepopi_1);LoadedData{2}.TS_CalcTime(:,keepopi_2)];
    if DidTrim, TS_CalcTime = TS_CalcTime(ix,:); end
    fprintf(1,' Done.\n');
end

% Now, save the results
% First check that HCTSA_loc.mat doesn't exist
fprintf(1,'A %u x %u matrix\n',size(TS_DataMat,1),size(TS_DataMat,2));
HereSheIs = which('HCTSA_loc.mat');
if isempty(HereSheIs)
    FileName = 'HCTSA_loc.mat';
        fprintf(1,['----------Saving to %s----------\n'],FileName)
else
    FileName = 'HCTSA_loc_combined.mat';
    fprintf(1,['----------Saving to %s----------\nYou''ll have to rename to HCTSA_loc.mat for' ...
    ' normal analysis routines like TSQ_normalize to work...\n'],FileName);
end

% Save
save(FileName,'TimeSeries','Operations','MasterOperations','-v7.3');
if GotData, save(FileName,'TS_DataMat','-append'); end % add data matrix
if GotQuality, save(FileName,'TS_Quality','-append'); end % add quality labels
if GotCalcTimes, save(FileName,'TS_CalcTime','-append'); end % add calculation times

fprintf(1,['Saved new Matlab file containing combined versions of %s ' ...
                    'and %s to %s\n'],HCTSA_locs{1},HCTSA_locs{2},FileName);
fprintf(1,'%s contains %u time series and %u operations\n',FileName, ...
                                length(TimeSeries),length(Operations));
    
end