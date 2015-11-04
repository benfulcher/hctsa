function TS_combine(HCTSA_loc_1,HCTSA_loc_2,compare_tsids,outputFileName)
% TS_combine   Combine two hctsa datasets.
%
% Takes a union of time series, and an intersection of operations from two hctsa
% datasets and writes the new combined dataset to a new .mat file
% Any data matrices are combined, and the structure arrays for TimeSeries and
% Operations are updated to reflect the concatenation.
%
% Note that in the case of duplicates, HCTSA_loc_1 takes precedence over
% HCTSA_loc_2.
%
% When using TS_init to generate datasets, be aware that the *same* set of
% operations and master operations must be used in both cases.
%
%---INPUTS:
% HCTSA_loc_1: the first hctsa dataset (a .mat filename)
% HCTSA_loc_2: the second hctsa dataset (a .mat filename)
% compare_tsids: (binary) whether to consider IDs in each file as the same.
%                If the two datasets to be joined are from different databases,
%                then a union of all time series results, regardless of the
%                uniqueness of their IDs (0, default).
%                However, if set to true (1, useful for different parts of a
%                dataset stored in the same mySQL database), IDs are matched so
%                that duplicate time series don't occur in the combined matrix.
% outputFileName: output to a custom .mat file ('HCTSA_loc.mat' by default).
%
%---OUTPUTS:
% Writes a new, combined .mat file (to the outputFileName)
%
%---USAGE:
% Combine two datasets into a new one:
% TS_combine('HCTSA_loc_1.mat','HCTSA_loc_2.mat');
%

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
    error('Must provide paths for two HCTSA_loc files')
end

if nargin < 3
    % If compare_tsids = 1, we assume that both are from the same database and thus
    % filter out any intersection between ts_ids in the two datasets
    compare_tsids = 0;
end
if compare_tsids == 1
    fprintf(1,['Assuming both %s and %s came from the same database so that ' ...
            'time series IDs are comparable.\nAny intersection of IDs will be filtered out.'],...
            HCTSA_loc_1,HCTSA_loc_2);
elseif compare_tsids == 0
    fprintf(1,['Assuming that %s and %s came different databases so' ...
            ' duplicate ts_ids can occur in the resulting matrix.\n'], ...
                                        HCTSA_loc_1,HCTSA_loc_2);
end

if nargin < 4
    outputFileName = 'HCTSA_loc.mat';
end
if ~strcmp(outputFileName(end-3:end),'.mat')
    error('Specify a .mat filename to output');
end


% ------------------------------------------------------------------------------
% Combine the local filenames
% ------------------------------------------------------------------------------
HCTSA_locs = {HCTSA_loc_1, HCTSA_loc_2};

% ------------------------------------------------------------------------------
% Check paths point to valid files
% ------------------------------------------------------------------------------
for i = 1:2
    if ~exist(HCTSA_locs{i},'file')
        error('Could not find %s',HCTSA_locs{i});
    end
end

% ------------------------------------------------------------------------------
% Load the two local files
% ------------------------------------------------------------------------------
fprintf(1,'Loading data...');
loadedData = cell(2,1);
for i = 1:2
    loadedData{i} = load(HCTSA_locs{i});
end
fprintf(1,' Loaded.\n');

% Give some information
for i = 1:2
    fprintf(1,'%u: The file, %s, contains information for %u time series and %u operations.\n', ...
                i,HCTSA_locs{i},length(loadedData{i}.TimeSeries),length(loadedData{i}.Operations));
end

%-------------------------------------------------------------------------------
% Check the fromDatabase flags
%-------------------------------------------------------------------------------
if loadedData{1}.fromDatabase ~= loadedData{2}.fromDatabase
    error('Weird that fromDatabase flags are inconsistent across the two HCTSA_loc files.');
else
    fromDatabase = loadedData{1}.fromDatabase;
end

% ------------------------------------------------------------------------------
% Construct a union of time series
% ------------------------------------------------------------------------------
% As a basic concatenation, then remove any duplicates

% First want to remove any additional fields
isExtraField = cellfun(@(x)~ismember(fieldnames(x.TimeSeries),{'ID','Name','Keywords', ...
                            'Length','Data'}),loadedData,'UniformOutput',0);
for i = 1:2
    if any(isExtraField{i})
        theextrafields = find(isExtraField{i});
        thefieldnames = fieldnames(loadedData{i}.TimeSeries);
        for j = 1:length(theextrafields)
            loadedData{i}.TimeSeries = rmfield(loadedData{i}.TimeSeries,...
                                                thefieldnames{theextrafields(j)});
            fprintf(1,'Extra field ''%s'' in %s\n',thefieldnames{theextrafields(j)},HCTSA_locs{i});
        end
    end
end
theFieldnames = fieldnames(loadedData{1}.TimeSeries);

% Now that fields should match the default fields, concatenate:
TimeSeries = cell2struct([squeeze(struct2cell(loadedData{1}.TimeSeries)), ...
                          squeeze(struct2cell(loadedData{2}.TimeSeries))], ...
                                theFieldnames);

%-------------------------------------------------------------------------------
% Check for time series duplicates
%-------------------------------------------------------------------------------
didTrim = 0;
[uniquetsids, ix_ts] = unique(vertcat(TimeSeries.ID));
if compare_tsids % ts_ids are comprable between the two files (i.e., retrieved from the same mySQL database)
    if ~fromDatabase
        warning('Be careful, we are assuming that time series IDs were assigned from a single TS_init (later split)')
    end
    if length(uniquetsids) < length(TimeSeries)
        fprintf(1,'We''re assuming that TimeSeries IDs are equivalent between the two input files\n');
        fprintf(1,'We need to trim duplicate time series (with the same IDs)\n');
        fprintf(1,['(NB: This will not be appropriate if combinining time series from' ...
                ' different databases, or produced using separate TS_init commands)\n']);
        fprintf(1,'Trimming %u duplicate time series to a total of %u\n', ...
                        length(TimeSeries)-length(uniquetsids),length(uniquetsids));
        TimeSeries = TimeSeries(ix_ts);
        didTrim = 1;
    else
        % Check for duplicate indices
        fprintf(1,'All time series were distinct, we now have a total of %u.\n',length(TimeSeries));
    end
else
    if length(uniquetsids) < length(TimeSeries)
        warning(['There are duplicate IDs in the time series -- you should run',
                 ' TS_ReIndex on the resulting .mat file'])
    end
end

% ------------------------------------------------------------------------------
% Construct an intersection of operations
% ------------------------------------------------------------------------------
% Check that the same number of operations if not from a database:
if ~fromDatabase
    numOperations = arrayfun(@(x)length(loadedData{x}.Operations),1:2);
    if ~(numOperations(1)==numOperations(2))
        error(['TS_init used to generate hctsa datasets with different numbers of\n' ...
                'operations; Operation IDs are not comparable.'])
    end
    numOperations = numOperations(1); % both the same

    % Check that all operation names match:
    namesMatch = arrayfun(@(x) strcmp(loadedData{1}.Name{x},loadedData{2}.Name{x}),...
                                    1:numOperations(1));
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
[allmopids,keepmopi_1] = intersect(vertcat(loadedData{1}.MasterOperations.ID),...
                                    vertcat(loadedData{2}.MasterOperations.ID));
MasterOperations = loadedData{1}.MasterOperations(keepmopi_1);

% ------------------------------------------------------------------------------
% 3. Data:
% ------------------------------------------------------------------------------
gotData = 0;
if isfield(loadedData{1},'TS_DataMat') && isfield(loadedData{2},'TS_DataMat')
    % Both contain data matrices
    gotData = 1;
    fprintf(1,'Combining data matrices...');
    TS_DataMat = [loadedData{1}.TS_DataMat(:,keepopi_1); loadedData{2}.TS_DataMat(:,keepopi_2)];
    if didTrim, TS_DataMat = TS_DataMat(ix_ts,:); end
    fprintf(1,' Done.\n');
end

% ------------------------------------------------------------------------------
% 4. Quality
% ------------------------------------------------------------------------------
gotQuality = 0;
if isfield(loadedData{1},'TS_Quality') && isfield(loadedData{2},'TS_Quality')
    gotQuality = 1;
    % Both contain Quality matrices
    fprintf(1,'Combining quality label matrices...');
    TS_Quality = [loadedData{1}.TS_Quality(:,keepopi_1); loadedData{2}.TS_Quality(:,keepopi_2)];
    if didTrim, TS_Quality = TS_Quality(ix_ts,:); end
    fprintf(1,' Done.\n');
end

% ------------------------------------------------------------------------------
% 5. Calculation times
% ------------------------------------------------------------------------------
gotCalcTimes = 0;
if isfield(loadedData{1},'TS_CalcTime') && isfield(loadedData{2},'TS_CalcTime')
    gotCalcTimes = 1;
    % Both contain Calculation time matrices
    fprintf(1,'Combining calculation time matrices...');
    TS_CalcTime = [loadedData{1}.TS_CalcTime(:,keepopi_1); loadedData{2}.TS_CalcTime(:,keepopi_2)];
    if didTrim, TS_CalcTime = TS_CalcTime(ix_ts,:); end
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
fprintf(1,['----------Saving to %s----------\n'],outputFileName);
if ~strcmp(outputFileName,'HCTSA_loc.mat')
    fprintf(1,['You''ll need to specify the custom filename ''%s'',' ...
                ' or rename to HCTSA_N.mat to run' ...
                ' normal analysis routines like TS_normalize...\n'],...
                outputFileName,outputFileName);
end

%--- Now actually save it:
save(outputFileName,'TimeSeries','Operations','MasterOperations','fromDatabase','-v7.3');
if gotData, save(outputFileName,'TS_DataMat','-append'); end % add data matrix
if gotQuality, save(outputFileName,'TS_Quality','-append'); end % add quality labels
if gotCalcTimes, save(outputFileName,'TS_CalcTime','-append'); end % add calculation times

%--- Tell the user what just happened:
fprintf(1,['Saved new Matlab file containing combined versions of %s' ...
                    ' and %s to %s\n'],HCTSA_locs{1},HCTSA_locs{2},outputFileName);
fprintf(1,'%s contains %u time series and %u operations.\n',outputFileName, ...
                                length(TimeSeries),length(Operations));

end
