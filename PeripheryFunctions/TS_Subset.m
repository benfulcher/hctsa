function [TS_DataMat,TimeSeries,Operations] = TS_Subset(whatData,ts_ids_keep,op_ids_keep,doSave,outputFileName)
% TS_Subset  Save a given subset of an hctsa dataset, based on time series and operation IDs
%
%---INPUTS:
% whatData, the source of the hctsa dataset: a .mat filename or shorthand
%           ('raw'/'norm'/..., default 'HCTSA_N.mat', cf. TS_LoadData), or a
%           struct of a pre-loaded dataset (e.g. load('HCTSA.mat')).
% ts_ids_keep, the IDs of time series to include in the subset (empty, [], to include all time series)
% op_ids_keep, the IDs of operations to include in the subset (empty, [], to include all operations)
% doSave, (binary), saves the result to file (default: true for a file source,
%           false for a struct source; a struct source with doSave=true needs
%           an explicit outputFileName).
% outputFileName, the filename to save the hctsa subset (if doSave==1).
%
%---OUTPUTS:
% TS_DataMat, TimeSeries, Operations: the subset of the hctsa dataset.
% (if doSave==true): the full subset of the hctsa dataset saved to file.
%
%---EXAMPLE USAGE:
% Import data from 'HCTSA_N.mat', then create a new dataset containing time
% series with IDs in the range 1--100, and all operations, saving the result
% back to 'HCTSA_subset.mat':
% >> TS_Subset('norm',1:100,[],1,'HCTSA_subset.mat');

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
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
    whatData = 'HCTSA_N.mat';
end
if nargin < 2
    ts_ids_keep = []; % all
end
if nargin < 3
    op_ids_keep = []; % all
end
if nargin < 4
    % Default to saving when given a file to work from; when given a pre-loaded
    % struct, default to just returning the subset in memory:
    doSave = ~isstruct(whatData);
end

if isempty(ts_ids_keep) && isempty(op_ids_keep)
    error('Nothing to subset!');
end

%-------------------------------------------------------------------------------
% Get the source dataset into a struct S (one read):
%-------------------------------------------------------------------------------
% NB: we take the matrices raw (not via TS_LoadData) and write them through
% unchanged -- a subset is a faithful slice of the source, so whatever coding
% the source uses for special/error cells (NaN from TS_Compute, or a legacy
% literal 0) is carried across as-is. The returned TS_DataMat is NaN-ified at
% the end for consistency with TS_LoadData (matching this function's historical
% return behaviour).
if isstruct(whatData)
    % Pre-loaded dataset (e.g. whatData = load('HCTSA.mat')):
    S = whatData;
    whatDataFile = '';
    if doSave && nargin < 5
        error(['TS_Subset: provide an outputFileName when subsetting a ' ...
               'pre-loaded struct with doSave=true.']);
    end
else
    switch whatData
    case {'raw','loc'}
        whatDataFile = 'HCTSA.mat';
    case {'norm','cl'}
        whatDataFile = 'HCTSA_N.mat';
    otherwise
        whatDataFile = whatData;
    end
    if ~strcmp(whatDataFile(end-3:end),'.mat')
        error('Specify a .mat filename');
    end
    if ~exist(whatDataFile,'file')
        error('%s not found',whatDataFile);
    end
    if nargin < 5
        outputFileName = regexprep(whatDataFile,'\.mat$','_subset.mat');
    end
    fprintf(1,'Loading data from %s...',whatDataFile);
    S = load(whatDataFile);
    fprintf(1,' Done.\n');
end

if doSave && ~strcmp(outputFileName(end-3:end),'.mat')
    error('Specify a .mat filename as output');
end

% Handle legacy structure-array metadata format:
if isstruct(S.TimeSeries)
    warning('Metadata stored in old structure-array format; run TS_ConvertToTable on your hctsa data file to update?');
    S.TimeSeries = struct2table(S.TimeSeries);
end
if isstruct(S.Operations)
    S.Operations = struct2table(S.Operations);
end

if ~isfield(S,'TS_DataMat')
    error('The source hctsa dataset does not contain a TS_DataMat');
end

numTimeSeries = height(S.TimeSeries);
numOperations = height(S.Operations);

%-------------------------------------------------------------------------------
% Do the subsetting:
%-------------------------------------------------------------------------------
i_keep = struct;
%--Match to TimeSeries IDs:
i_keep.TimeSeries = MatchMe(S.TimeSeries.ID,ts_ids_keep);
%--Match to Operation IDs:
i_keep.Operations = MatchMe(S.Operations.ID,op_ids_keep);

TimeSeries = S.TimeSeries(i_keep.TimeSeries,:);
Operations = S.Operations(i_keep.Operations,:);

% Subset every stored (time series)x(operations) matrix:
for theField = {'TS_DataMat','TS_Quality','TS_CalcTime'}
    if isfield(S,theField{1})
        S.(theField{1}) = S.(theField{1})(i_keep.TimeSeries,i_keep.Operations);
    end
end
TS_DataMat = S.TS_DataMat; % raw: special/error cells carried through unchanged

fprintf('The hctsa dataset now contains %u -> %u time series and %u -> %u operations.\n',...
            numTimeSeries,height(TimeSeries),numOperations,height(Operations))

% Remove group information because this will no longer be valid for sure
if ~isempty(ts_ids_keep) && ismember('Group',TimeSeries.Properties.VariableNames)
    TimeSeries(:,'Group') = [];
    fprintf('Warning: group information removed -- regenerate for subset data using TS_LabelGroups\n')
end

if doSave
    %---------------------------------------------------------------------------
    % Write a single fresh .mat file: the subsetted tables/matrices, plus every
    % other variable in the source carried forward unchanged (gitInfo,
    % MasterOperations, normalisationInfo, ...). No copyfile of the full source,
    % no save('-append') (which on -v7.3 leaves dead space for overwritten
    % variables), and only one write of exactly the subset.
    %---------------------------------------------------------------------------
    S.TimeSeries = TimeSeries;
    S.Operations = Operations;

    % Clustering info no longer valid for a subset -- reset it (if present):
    if isfield(S,'ts_clust')
        S.ts_clust = struct('distanceMetric','none','Dij',[],...
                    'ord',1:height(TimeSeries),'linkageMethod','none');
    end
    if isfield(S,'op_clust')
        S.op_clust = struct('distanceMetric','none','Dij',[],...
                    'ord',1:height(Operations),'linkageMethod','none');
    end

    save(outputFileName,'-struct','S','-v7.3');

    fprintf(1,'Data saved to %s!\n',outputFileName);
end

if nargout == 0
    % Don't display all of this info to screen if it's been saved and not stored
    clear('TS_DataMat','TimeSeries','Operations');
elseif isfield(S,'TS_Quality')
    % NaN-ify special/error cells in the returned matrix, matching TS_LoadData:
    TS_DataMat(~isfinite(TS_DataMat)) = NaN;
    TS_DataMat(S.TS_Quality > 0) = NaN;
end

%-------------------------------------------------------------------------------
function ind = MatchMe(idsAll,idsMatch)
    % Find matches:
    if isempty(idsMatch)
        ind = true(size(idsAll));
    else
        ind = ismember(idsAll,idsMatch);
    end
end

end
