function TS_combine(HCTSA_1,HCTSA_2,compare_tsids,merge_features,outputFileName)
% TS_combine   Combine two hctsa datasets (same features, different data)
%
% Takes a union of time series, and an intersection of operations from two hctsa
% datasets and writes the new combined dataset to a new .mat file
% Any data matrices are combined, and the structure arrays for TimeSeries and
% Operations are updated to reflect the concatenation.
%
% Note that in the case of duplicates, HCTSA_1 takes precedence over
% HCTSA_2.
%
% When using TS_init to generate datasets, be aware that the *same* set of
% operations and master operations must be used in both cases.
%
% NB: Use TS_merge if same data, different features
%
%---INPUTS:
% HCTSA_1: the first hctsa dataset (a .mat filename)
% HCTSA_2: the second hctsa dataset (a .mat filename)
% compare_tsids: (logical) whether to consider IDs in each file as the same.
%                If the two datasets to be joined are from different databases,
%                then a union of all time series results, regardless of the
%                uniqueness of their IDs (false, default).
%                However, if set to true (true, useful for different parts of a
%                dataset stored in the same mySQL database), IDs are matched so
%                that duplicate time series don't occur in the combined matrix.
% merge_features: (logical) whehter to merge distinct feature sets that occur
%                   between the two datasets. By default (false) assumes that
%                   both datasets were computed using the exact same feature set.
%                   Setting to true takes the union of the features present in
%                   both (assumes disjoint).
% outputFileName: output to a custom .mat file ('HCTSA.mat' by default).
%
%---OUTPUTS:
% Writes a new, combined .mat file (to the outputFileName)
%
%---USAGE:
% Combine two datasets computed using the same set of features, 'HCTSA_1.mat' and
% 'HCTSA_2.mat', into a new combined HCTSA file:
% TS_combine('HCTSA_1.mat','HCTSA_2.mat');

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check inputs:
% ------------------------------------------------------------------------------

if nargin < 2
    error('Must provide paths for two HCTSA*.mat files')
end

if nargin < 3
    % If compare_tsids = 1, we assume that both are from the same database and thus
    % filter out any intersection between ts_ids in the two datasets
    compare_tsids = false;
end
if compare_tsids
    fprintf(1,['Assuming both %s and %s came from the same database so that ' ...
            'time series IDs are comparable.\nAny intersection of IDs will be filtered out.\n'],...
            HCTSA_1,HCTSA_2);
else
    fprintf(1,['Assuming that %s and %s came different databases so' ...
            ' duplicate ts_ids can occur in the resulting matrix.\n'], ...
                                        HCTSA_1,HCTSA_2);
end
if nargin < 4
    % By default assumes identical feature sets between the two files:
    merge_features = false;
    % (If true, all time series must be identical)
end
if nargin < 5
    outputFileName = 'HCTSA.mat';
end
if ~strcmp(outputFileName(end-3:end),'.mat')
    error('Specify a .mat filename to output');
end


% ------------------------------------------------------------------------------
% Combine the local filenames
% ------------------------------------------------------------------------------
HCTSAs = {HCTSA_1, HCTSA_2};

% ------------------------------------------------------------------------------
% Check paths point to valid files
% ------------------------------------------------------------------------------
for i = 1:2
    if ~exist(HCTSAs{i},'file')
        error('Could not find %s',HCTSAs{i});
    end
end

% ------------------------------------------------------------------------------
% Load the two local files
% ------------------------------------------------------------------------------
fprintf(1,'Loading data...');
loadedData = cell(2,1);
for i = 1:2
    loadedData{i} = load(HCTSAs{i});
end
fprintf(1,' Loaded.\n');

% Give some information
for i = 1:2
    fprintf(1,'%u: The file, %s, contains information for %u time series and %u operations.\n', ...
                i,HCTSAs{i},length(loadedData{i}.TimeSeries),length(loadedData{i}.Operations));
end

%-------------------------------------------------------------------------------
% Check the fromDatabase flags
%-------------------------------------------------------------------------------
if loadedData{1}.fromDatabase ~= loadedData{2}.fromDatabase
    error('Weird that fromDatabase flags are inconsistent between the two HCTSA files.');
else
    fromDatabase = loadedData{1}.fromDatabase;
end

%-------------------------------------------------------------------------------
% Check the git data
%-------------------------------------------------------------------------------
if ~isfield(loadedData{1},'gitInfo') || ~isfield(loadedData{2},'gitInfo')
    % git info not present in at least one -- keep an empty structure
    gitInfo = struct();
elseif isempty(fieldnames(loadedData{1}.gitInfo)) && isempty(fieldnames(loadedData{2}.gitInfo))
    % git info not stored in either HCTSA file
    gitInfo = struct();
elseif ~strcmp(loadedData{1}.gitInfo.hash,loadedData{2}.gitInfo.hash)
    % Only check the hashes for consistency:
    error('Git versions are inconsistent between the two HCTSA files.');
else
    gitInfo = loadedData{1}.gitInfo;
end

%-------------------------------------------------------------------------------
% Remove any additional fields from the TimeSeries structure array:
%-------------------------------------------------------------------------------
isExtraField = cellfun(@(x)~ismember(fieldnames(x.TimeSeries),{'ID','Name','Keywords', ...
                            'Length','Data'}),loadedData,'UniformOutput',0);
for i = 1:2
    if any(isExtraField{i})
        theExtraFields = find(isExtraField{i});
        theFieldnames = fieldnames(loadedData{i}.TimeSeries);
        for j = 1:length(theExtraFields)
            loadedData{i}.TimeSeries = rmfield(loadedData{i}.TimeSeries,...
                                                theFieldnames{theExtraFields(j)});
            fprintf(1,'Extra field ''%s'' in %s\n',theFieldnames{theExtraFields(j)},HCTSAs{i});
        end
    end
end
theFieldnames = fieldnames(loadedData{1}.TimeSeries);

needReIndex = false; % whether you need to reindex the result (combined datasets from different indexes)

if merge_features
    % Time-series data are identical; features are disjoint

    %===============================================================================
    % Check that all time series are identical:
    %===============================================================================
    numTimeSeries = arrayfun(@(x)length(loadedData{x}.TimeSeries),1:2);
    if ~(numTimeSeries(1)==numTimeSeries(2))
        error(sprintf(['hctsa datasets contain different numbers of\n' ...
                'time series; TimeSeries IDs are not comparable.']))
    end
    numTimeSeries = numTimeSeries(1); % both the same

    % Check that all TimeSeries names match:
    namesMatch = arrayfun(@(x) strcmp(loadedData{1}.TimeSeries(x).Name,...
                                      loadedData{2}.TimeSeries(x).Name),...
                                      1:numTimeSeries);
    if ~all(namesMatch)
        keyboard
        error('The names of time series in the two files do not match');
    end

    % Ok so same number of time series in both, and all names match:
    TimeSeries = loadedData{1}.TimeSeries; % identical; keep all

    %===============================================================================
    % Construct a union of operations
    %===============================================================================
    [sameOperations, ~] = intersect({loadedData{1}.Operations.Name},{loadedData{2}.Operations.Name});
    if ~isempty(sameOperations)
        error('Some operations overlap between the two files :-/');
    end

    % All unique, so can concatenate:
    Operations = cell2struct([squeeze(struct2cell(loadedData{1}.Operations)), ...
                              squeeze(struct2cell(loadedData{2}.Operations))], ...
                                     fieldnames(loadedData{1}.Operations));

    % ------------------------------------------------------------------------------
    % Construct a union of MasterOperations
    % ------------------------------------------------------------------------------
    [sameMOperations, ~] = intersect({loadedData{1}.MasterOperations.Label},{loadedData{2}.MasterOperations.Label});
    if ~isempty(sameMOperations)
        error('Some master operations overlap between the two files :-/');
    end
    MasterOperations = cell2struct([squeeze(struct2cell(loadedData{1}.MasterOperations)), ...
                              squeeze(struct2cell(loadedData{2}.MasterOperations))], ...
                                     fieldnames(loadedData{1}.MasterOperations));

else
    %===============================================================================
    %===============================================================================
    % Time-series data are distinct; features overlap
    %===============================================================================
    %===============================================================================

    %-------------------------------------------------------------------------------
    % Construct a union of time series
    %-------------------------------------------------------------------------------
    % As a basic concatenation, then remove any duplicates

    % Fields should match the default fields, so can concatenate:
    TimeSeries = cell2struct([squeeze(struct2cell(loadedData{1}.TimeSeries)), ...
                              squeeze(struct2cell(loadedData{2}.TimeSeries))], ...
                                    theFieldnames);

    %-------------------------------------------------------------------------------
    % Check for time series duplicates
    %-------------------------------------------------------------------------------
    didTrim = false; % whether you remove time series (that appear in both hctsa data files)

    if compare_tsids % TimeSeries IDs are comparable between the two files (i.e., retrieved from the same mySQL database)
        [uniquetsids, ix_ts] = unique(vertcat(TimeSeries.ID)); % will be sorted
        TimeSeries = TimeSeries(ix_ts);
        if ~fromDatabase
            fprintf(1,'Be careful, we are assuming that time series IDs were assigned from a *single* TS_init\n')
        end
        % Check for duplicate indices:
        if length(uniquetsids) < length(TimeSeries)
            fprintf(1,'We''re assuming that TimeSeries IDs are equivalent between the two input files\n');
            fprintf(1,'We need to trim duplicate time series (with the same IDs)\n');
            fprintf(1,['(NB: This will NOT be appropriate if combinining time series from' ...
                    ' different databases, or produced using separate TS_init commands)\n']);
            fprintf(1,'Trimming %u duplicate time series to a total of %u\n', ...
                            length(TimeSeries)-length(uniquetsids),length(uniquetsids));
            didTrim = true;
        else
            fprintf(1,'All time series were distinct, we now have a total of %u.\n',length(TimeSeries));
        end
    else
        % Check that time series names are unique, and trim if not:
        [uniqueTimeSeriesNames, ix_ts] = unique({TimeSeries.Name},'stable');
        TimeSeries = TimeSeries(ix_ts);
        numUniqueTimeSeries = length(uniqueTimeSeriesNames);
        if numUniqueTimeSeries < length(TimeSeries)
            warning('%u duplicate time series names present in combined dataset -- removed',...
                                    length(TimeSeries) - numUniqueTimeSeries);
            didTrim = true; % will trim the data matrix and other such matrices with ix_ts later
        end

        % Now see if there are duplicate IDs (meaning that we need to reindex):
        uniquetsids = unique(vertcat(TimeSeries.ID));
        if length(uniquetsids) < length(TimeSeries)
            needReIndex = true;
            % This is done at the end (after saving all data)
        end
    end

    % ------------------------------------------------------------------------------
    % Construct an intersection of operations
    % ------------------------------------------------------------------------------
    % Check that the same number of operations if not from a database:
    if ~fromDatabase
        numOperations = arrayfun(@(x)length(loadedData{x}.Operations),1:2);
        if ~(numOperations(1)==numOperations(2))
            error(sprintf(['TS_init used to generate hctsa datasets with different numbers of\n' ...
                    'operations; Operation IDs are not comparable.']))
        end
        numOperations = numOperations(1); % both the same

        % Check that all operation names match:
        namesMatch = arrayfun(@(x) strcmp(loadedData{1}.Operations(x).Name,...
                                          loadedData{2}.Operations(x).Name),...
                                          1:numOperations);
        if ~all(namesMatch)
            error('TS_init used to generate hctsa datasets, and the names of operations do not match');
        end

        % Ok so same number of operations in both, and all names match:
        keepopi_1 = 1:numOperations; % range to keep for both is the same
        keepopi_2 = 1:numOperations; % range to keep for both is the same
        Operations = loadedData{1}.Operations; % keep all

    else
        % --Both datasets are from a database (assume the same database, or
        % the same operation IDs in both databases)

        % Take intersection of operation ids, and take information from first input
        [~,keepopi_1,keepopi_2] = intersect(vertcat(loadedData{1}.Operations.ID),...
                                            vertcat(loadedData{2}.Operations.ID));

        % Data from first file goes in (should be identical to keepopi_2 in the second file)
        Operations = loadedData{1}.Operations(keepopi_1);

        fprintf(1,'Keeping the %u overlapping operations.\n',length(Operations));
    end

    % ------------------------------------------------------------------------------
    % Construct an intersection of MasterOperations
    % ------------------------------------------------------------------------------
    % Take intersection, like operations -- those that are in both
    [~,keepmopi_1] = intersect(vertcat(loadedData{1}.MasterOperations.ID),...
                                        vertcat(loadedData{2}.MasterOperations.ID));
    MasterOperations = loadedData{1}.MasterOperations(keepmopi_1);

end

% ------------------------------------------------------------------------------
% 3. Data:
% ------------------------------------------------------------------------------
if isfield(loadedData{1},'TS_DataMat') && isfield(loadedData{2},'TS_DataMat')
    gotData = true;
else
    gotData = false;
end

if gotData
    % Both hctsa files contain data matrices
    fprintf(1,'Combining data matrices...');
    if merge_features
        % time series are identical, operations are distinct: simple merge:
        TS_DataMat = [loadedData{1}.TS_DataMat,loadedData{2}.TS_DataMat];
    else
        % Time series are distinct; operations made to match as intersection of two datasets
        TS_DataMat = [loadedData{1}.TS_DataMat(:,keepopi_1); loadedData{2}.TS_DataMat(:,keepopi_2)];
        if didTrim
            % Trimmed time series above, need to make sure the data matches now:
            TS_DataMat = TS_DataMat(ix_ts,:);
        end
    end
    fprintf(1,' Done.\n');
end

% ------------------------------------------------------------------------------
% 4. Quality
% ------------------------------------------------------------------------------
if isfield(loadedData{1},'TS_Quality') && isfield(loadedData{2},'TS_Quality')
    gotQuality = true;
else
    gotQuality = false;
end
if gotQuality
    % Both contain quality matrices
    fprintf(1,'Combining quality label matrices...');
    if merge_features
        TS_Quality = [loadedData{1}.TS_Quality,loadedData{2}.TS_Quality];
    else
        TS_Quality = [loadedData{1}.TS_Quality(:,keepopi_1); loadedData{2}.TS_Quality(:,keepopi_2)];
        if didTrim
            % Trimmed time series above, need to make sure the data matches now:
            TS_Quality = TS_Quality(ix_ts,:);
        end
    end
    fprintf(1,' Done.\n');
end

% ------------------------------------------------------------------------------
% 5. Calculation times
% ------------------------------------------------------------------------------
if isfield(loadedData{1},'TS_CalcTime') && isfield(loadedData{2},'TS_CalcTime')
    gotCalcTimes = true;
else
    gotCalcTimes = false;
end
if gotCalcTimes
    % Both contain Calculation time matrices
    fprintf(1,'Combining calculation time matrices...');
    if merge_features
        TS_CalcTime = [loadedData{1}.TS_CalcTime,loadedData{2}.TS_CalcTime];
    else
        TS_CalcTime = [loadedData{1}.TS_CalcTime(:,keepopi_1); loadedData{2}.TS_CalcTime(:,keepopi_2)];
        if didTrim
            % Trimmed time series above, need to make sure the data matches now:
            TS_CalcTime = TS_CalcTime(ix_ts,:);
        end
    end
    fprintf(1,' Done.\n');
end

% ------------------------------------------------------------------------------
% Save the results
% ------------------------------------------------------------------------------
% First check that the output file doesn't already exist:
fprintf(1,'A %u x %u matrix\n',size(TS_DataMat,1),size(TS_DataMat,2));
hereSheIs = which(fullfile(pwd,outputFileName));
if ~isempty(hereSheIs) % already exists
    outputFileName = [outputFileName(1:end-4),'_combined.mat'];
end
fprintf(1,'----------Saving to %s----------\n',outputFileName);

%--- Now actually save it:
save(outputFileName,'TimeSeries','Operations','MasterOperations','fromDatabase','gitInfo','-v7.3');
if gotData, save(outputFileName,'TS_DataMat','-append'); end % add data matrix
if gotQuality, save(outputFileName,'TS_Quality','-append'); end % add quality labels
if gotCalcTimes, save(outputFileName,'TS_CalcTime','-append'); end % add calculation times

%--- Tell the user what just happened:
fprintf(1,['Saved new Matlab file containing combined versions of %s' ...
                    ' and %s to %s\n'],HCTSAs{1},HCTSAs{2},outputFileName);
fprintf(1,'%s contains %u time series and %u operations.\n',outputFileName, ...
                                length(TimeSeries),length(Operations));

%-------------------------------------------------------------------------------
% ReIndex??
%-------------------------------------------------------------------------------
if needReIndex
    fprintf(1,'There are duplicate IDs in the time series -- we need to reindex\n');
    TS_ReIndex(outputFileName,'ts',true);
end

end
