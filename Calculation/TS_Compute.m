function TS_Compute(parallelHow,ts_id_range,op_id_range,computeWhat,customFile,howVocal,beVerbose,autosaveEvery)
% TS_Compute    Computes missing elements of TS_DataMat
%
%---EXAMPLE USAGE:
% TS_Compute;
%
%---INPUTS:
% parallelHow: how to use Matlab's Parallel Computing Toolbox:
%               'none'       (default): no parallelization (serial for-loops).
%               'operations': parallelize over the ~1000 master operations
%                             *within* each time series (looping time series
%                             serially). Best when there are few, long time
%                             series.
%               'series':     parallelize over time series (each worker computes
%                             whole feature vectors; operations run serially per
%                             worker). Best for "many short time series"
%                             workloads -- near-linear speedup.
%               For backwards compatibility, logical false maps to 'none' and
%               true maps to 'operations'.
% ts_id_range: a custom range of time series IDs to compute (default: [] -- compute all)
% op_id_range: a custom range of operation IDs to compute (default: [] -- compute all)
% computeWhat: whether to compute just missing values ('missing', default), or
% 				ALSO retry results that previously threw an error ('error'), or
% 				ALSO retry any result that previously did not return a good value ('bad')
% customFile: reads in and writes to a custom output file
% howVocal:   {'fast','minimal','full'}: how much output to show about the calculation of operations.
% beVerbose:  [default: false] whether to show text output (fprintf/warning messages)
% 				generated inside individual Operations, rather than suppressing it.
% autosaveEvery: how often (seconds) to write partial progress to disk during a
%               long run, so a killed job keeps its progress and can be resumed
%               with the same call. Default Inf -- OFF: results are saved once,
%               at the end (as TS_Compute has always behaved). Set e.g. 1800 for
%               long cluster runs. Rewriting the whole matrices is not free on
%               big datasets, so the effective interval is throttled up
%               automatically to keep the save overhead small.
%
%---OUTPUTS:
% Writes output to customFile (HCTSA.mat by default)
%
%---NOTES:
% When 'series' starts a parallel pool itself, workers on a *local* pool inherit
% the client's path and toolbox licences; a remote-cluster pool must be set up
% (path + licences) by the caller before calling TS_Compute.

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

% --------------------------------------------------------------------------
%% Check inputs and set defaults
% --------------------------------------------------------------------------

% How to use the Parallel Computing Toolbox?
if nargin < 1 || isempty(parallelHow)
	parallelHow = 'none';
end
if islogical(parallelHow) || isnumeric(parallelHow)
    % Backwards compatibility with the old logical doParallel argument:
    if parallelHow
        parallelHow = 'operations';
    else
        parallelHow = 'none';
    end
end
parallelHow = validatestring(lower(parallelHow),{'none','operations','series'});
doParallelOps    = strcmp(parallelHow,'operations');
doParallelSeries = strcmp(parallelHow,'series');

% Custom range of ts_ids to compute
if nargin < 2
    ts_id_range = []; % compute all ts_ids in the file by default
end

% Custom range of op_ids to compute
if nargin < 3
    op_id_range = []; % compute all op_ids in the file by default
end

if nargin < 4 || isempty(computeWhat)
	computeWhat = 'missing';
end
if ~ismember(computeWhat,{'missing','error','bad'})
	error('Unknown setting ''%s''',computeWhat);
end

% Custom HCTSA.mat file:
if nargin < 5 || isempty(customFile)
    customFile = 'raw';
end

% Be vocal?
if nargin < 6
    howVocal = 'full'; % Write back full information on all calculations to screen
    % prints every piece of code evaluated (nice for error checking)
end

% Show text output generated inside individual Operations?
if nargin < 7 || isempty(beVerbose)
    beVerbose = false;
end

% How often (seconds) to autosave partial progress. Default Inf = off (a single
% save at the end, as TS_Compute has always done):
if nargin < 8 || isempty(autosaveEvery)
    autosaveEvery = Inf;
end
lastSaveSeconds = 0; % duration of the most recent save, to throttle autosaves

% --------------------------------------------------------------------------
%% Load information from local files
% --------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,customFile] = TS_LoadData(customFile);
fileVarsStruct = whos('-file',customFile);
fileVars = {fileVarsStruct.name};
if ~all(ismember({'TimeSeries','Operations','MasterOperations','TS_DataMat'},fileVars))
	error('\nCannot compute on %s: there are missing variables',customFile);
end
MasterOperations = TS_GetFromData(customFile,'MasterOperations');
% Precompile each master operation's code string into a function handle once,
% here (outside the per-time-series loop below), instead of having eval/evalc
% re-parse the same text from scratch for every time series in
% TS_ComputeMasterLoop:
MasterOperations.Fn = cell(height(MasterOperations),1);
for i = 1:height(MasterOperations)
    MasterOperations.Fn{i} = str2func(['@(x,x_z) ',MasterOperations.Code{i}]);
end
if ismember('TS_CalcTime',fileVars)
	TS_CalcTime = TS_GetFromData(customFile,'TS_CalcTime');
end
if ismember('TS_Quality',fileVars)
	TS_Quality = TS_GetFromData(customFile,'TS_Quality');
end

%-------------------------------------------------------------------------------
% For 'series' parallelism, make sure a pool is available (fall back to serial):
%-------------------------------------------------------------------------------
if doParallelSeries
    if ~TS_InitiateParallel(true)
        warning(['Could not start a parallel pool -- computing time series ' ...
                 'serially instead of in parallel.']);
        doParallelSeries = false;
    end
end

%-------------------------------------------------------------------------------
% Tell the user about settings
%-------------------------------------------------------------------------------
if doParallelSeries
    fprintf(1,['Computation will be performed across multiple cores using ' ...
            'Matlab''s Parallel Computing Toolbox (parallelized over time series).\n']);
elseif doParallelOps
    fprintf(1,['Computation will be performed across multiple cores using ' ...
            'Matlab''s Parallel Computing Toolbox (parallelized over operations).\n']);
else % use single-threaded for loops
    fprintf(1,'Computations will be performed serially without parallelization.\n');
end

% ------------------------------------------------------------------------------
% Get indices if computing a subset
% ------------------------------------------------------------------------------
allIDs = TimeSeries.ID;
if isempty(ts_id_range)
    ts_id_range = allIDs;
    tsIndex = 1:height(TimeSeries);
else
    ts_id_range = intersect(ts_id_range,allIDs);
    tsIndex = find(ismember(allIDs,ts_id_range));
end
allIDs = Operations.ID;
if isempty(op_id_range)
    op_id_range = allIDs;
    opCompute = ones(1,height(Operations));
else
    op_id_range = intersect(op_id_range,allIDs);
    opCompute = ismember(allIDs,op_id_range)';
end
opCompute = logical(opCompute);

% Definitions
numTimeSeries = length(ts_id_range); % Number of time series
numOps = length(op_id_range); % Number of operations

% Check that some computable range exists
if numTimeSeries==0 || numOps==0
    fprintf(1,'%u time series and %u operations match the ids provided. Exiting.\n',...
						numTimeSeries,numOps);
    return
end

if (nargin < 6) && (numOps < 30)
    howVocal = 'fast';
end

%-------------------------------------------------------------------------------
fprintf(1,'[%s]: Extracting %u features from each of %u time series.\n',...
                    datestr(now),numOps,numTimeSeries);


% The times vector stores the time taken for each time series to have its
% operations calculated (for determining time remaining)
times = zeros(numTimeSeries,1);

% Initialize TS_CalcTime and TS_Quality if they don't yet exist
if ~exist('TS_CalcTime','var')
    TS_CalcTime = zeros(size(TS_DataMat));
end
if ~exist('TS_Quality','var')
    TS_Quality = zeros(size(TS_DataMat));
end

% --------------------------------------------------------------------------
%% Computation
% --------------------------------------------------------------------------
% (beVerbose's per-feature messages, propagated down from TS_CalculateFeatureVector,
% don't mix well with an animated progress bar, so skip it in that case)
showProgressBar = strcmp(howVocal,'fast') && ~beVerbose;
if showProgressBar
    % Just show a progress bar over time series (not across features for a given time series)
    % makes sense for, say, catch22, where computations are super fast
    BF_ProgressBar('new')
end

numCalc_all = zeros(numTimeSeries,1);
lastSaveTimer = tic; % time since progress was last written to disk

if doParallelSeries
    % ======================================================================
    % Parallel over time series (operations serial on each worker)
    % ======================================================================
    pool = gcp;
    blockSize = min(numTimeSeries, max(32, 8*pool.NumWorkers));

    % Broadcast the operation tables to the workers once (not per iteration):
    opsConstant  = parallel.pool.Constant(Operations);
    mopsConstant = parallel.pool.Constant(MasterOperations);

    for blockStart = 1:blockSize:numTimeSeries
        blockPos  = blockStart:min(blockStart+blockSize-1,numTimeSeries); % positions in 1:numTimeSeries
        blockRows = tsIndex(blockPos); % rows into the TS_* matrices
        nBlock    = numel(blockPos);

        % Which operations does each series in this block need?
        blockToCalc = false(nBlock,size(TS_Quality,2));
        for k = 1:nBlock
            ql = TS_Quality(blockRows(k),:);
            switch computeWhat
            case 'missing'
                blockToCalc(k,:) = opCompute & isnan(ql);
            case 'error'
                blockToCalc(k,:) = opCompute & (isnan(ql) | ql == 1);
            case 'bad'
                blockToCalc(k,:) = opCompute & (isnan(ql) | ql > 0);
            end
        end

        blockData = TimeSeries.Data(blockRows); % raw time-series vectors only

        fv = cell(nBlock,1); ct = cell(nBlock,1); cq = cell(nBlock,1);
        elapsed = zeros(nBlock,1);
        failMsg = repmat({''},nBlock,1);

        parfor k = 1:nBlock
            toCalc_k = blockToCalc(k,:);
            if ~any(toCalc_k)
                continue
            end
            seriesTimer = tic;
            try
                [fv{k},ct{k},cq{k}] = TS_CalculateFeatureVector(blockData{k},false,...
                        opsConstant.Value(toCalc_k,:),mopsConstant.Value,true,'fast',beVerbose);
            catch ME
                failMsg{k} = ME.message;
            end
            elapsed(k) = toc(seriesTimer);
        end

        % Store the block's results back into the full matrices:
        for k = 1:nBlock
            numCalc_all(blockPos(k)) = sum(blockToCalc(k,:));
            times(blockPos(k)) = elapsed(k);
            if ~isempty(failMsg{k})
                warning('Calculation for time series %u / %u failed:\n%s',...
                            blockPos(k),numTimeSeries,failMsg{k})
                continue
            end
            toCalc_k = blockToCalc(k,:);
            if ~any(toCalc_k)
                continue
            end
            TS_DataMat(blockRows(k),toCalc_k)  = fv{k};
            TS_CalcTime(blockRows(k),toCalc_k) = ct{k};
            TS_Quality(blockRows(k),toCalc_k)  = cq{k};
        end

        % Progress to the user:
        numDone = blockPos(end);
        if showProgressBar
            BF_ProgressBar(numDone/numTimeSeries,[],[],sprintf(' %.0f%% complete',100*numDone/numTimeSeries))
        elseif ~strcmp(howVocal,'minimal')
            meanTime = mean(times(1:numDone));
            fprintf(1,'- - - %u / %u time series done (~%s remaining) - - -\n',...
                    numDone,numTimeSeries,BF_TheTime((numTimeSeries-numDone)*meanTime/pool.NumWorkers,1));
        end

        % Autosave progress between blocks (if enabled). The effective interval
        % is throttled to >= 15x the last save's duration so rewriting large
        % matrices can't dominate the run:
        if numDone < numTimeSeries && any(numCalc_all > 0) && ...
                toc(lastSaveTimer) > max(autosaveEvery,15*lastSaveSeconds)
            lastSaveSeconds = saveResults(customFile,TS_DataMat,TS_CalcTime,TS_Quality,true);
            lastSaveTimer = tic;
        end
    end

else
    % ======================================================================
    % Serial over time series (operations optionally parallel)
    % ======================================================================
    for i = 1:numTimeSeries

        tsInd = tsIndex(i);
        bigTimer = tic;

        % ----
        % Which operations need calculating for this time series?:
        % ----
        qualityLabels = TS_Quality(tsInd,:); % The calculation states of any existing results for the current time series, a line of TS_Quality
              % NaN indicates a value never before calculated, 1 indicates fatal error before (try again now)

        % Determine which operations are awaiting calculation for this time series:
        switch computeWhat
        case 'missing'
            % try to compute missing values (i.e, never previously computed for this time series)
            toCalc = (opCompute & isnan(qualityLabels));
        case 'error'
            % compute missing or previously threw an error
            toCalc = (opCompute & (isnan(qualityLabels) | qualityLabels == 1));
        case 'bad'
            % compute missing, or anything that wasn't previously a good value
            toCalc = (opCompute & (isnan(qualityLabels) | qualityLabels > 0)); % Operations awaiting calculation
        end
        numCalc = sum(toCalc); % Number of operations to evaluate
        numCalc_all(i) = numCalc; % keep a record of how many were calculated at each iteration

        % -----
        % Check that all operations have a master ID attached:
        % -----
        if sum(toCalc) < numCalc
            % Error in the database structure; some operations are missing MasterID assignment
            error('?? Database structure error: some operations have not been assigned a valid master operation');
        end

        switch howVocal
        case 'minimal'
            fprintf(1,'- - - - - - - Time series %u / %u: %s - - - - -\n',i,numTimeSeries,TimeSeries.Name{tsInd});
        case 'full'
            fprintf(1,'\n\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf(1,'; ; ; : : : : ; ; ;    %s     ; ; ; : : : ; ; ;\n',datestr(now));
            fprintf(1,'- - - - - - - - - - - Time series %u / %u - - - - - - - - - - -\n',i,numTimeSeries);
        end

        if numCalc > 0 % some to calculate
            try
                [featureVector,calcTimes,calcQuality] = TS_CalculateFeatureVector(TimeSeries(tsInd,:),...
                                doParallelOps,Operations(toCalc,:),MasterOperations,true,howVocal,beVerbose);
            catch ME
                % Skip to the next time series; the entries for this time series in TS_DataMat etc. will remain NaNs
                warning('Calculation for time series %u / %u failed:\n%s',i,numTimeSeries,ME.message)
                continue
            end

            % ------------------------------------------------------------------------------
            %% Store the calculated information back to local matrices
            % ------------------------------------------------------------------------------
            TS_DataMat(tsInd,toCalc) = featureVector; % store outputs in TS_DataMat
            TS_CalcTime(tsInd,toCalc) = calcTimes; % store calculation times in TS_CalcTime
            TS_Quality(tsInd,toCalc) = calcQuality; % store quality labels in TS_Quality
            % NB: the calculation time assigned for individual operations is the total calculation
            % time taken to evaluate the master code.
        else
            fprintf(1,'Nothing calculated! All %u operations already complete!!\n',numOps);
        end

        % The time taken to calculate (or not, if numCalc = 0) all operations for this time series:
        times(i) = toc(bigTimer);

        % Update progress to the user
        switch howVocal
        case 'fast'
            if showProgressBar
                % (a time-based ETA is a poor estimate here: Operations are computed
                % in file order, and the cheap ones are clustered first -- so a
                % linear extrapolation from elapsed time is systematically optimistic
                % until the slow operations near the end. Percentage complete is
                % exact, so isn't subject to that bias.)
                BF_ProgressBar(i/numTimeSeries,[],[],sprintf(' %.0f%% complete',100*i/numTimeSeries))
            end
        case 'minimal'

        case 'full'
            if i < numTimeSeries
                fprintf(1,'- - - - - - - -  %u time series remaining - - - - - - - -\n',numTimeSeries-i);
                fprintf(1,'- - - - - - - -  %s remaining - - - - - - - - -\n', ...
                                                    BF_TheTime(((numTimeSeries-i)*mean(times(1:i))),1));
            else % The final time series
                fprintf(1,'- - - - - - - - - - All %u time series calculated! - - - - - - - - - -\n', ...
                                                            numTimeSeries);
            end
            fprintf(1,'********************************************************************\n');
        end

        % Autosave progress (if enabled; throttled to >= 15x the last save's
        % duration so rewriting large matrices can't dominate the run):
        if i < numTimeSeries && any(numCalc_all > 0) && ...
                toc(lastSaveTimer) > max(autosaveEvery,15*lastSaveSeconds)
            lastSaveSeconds = saveResults(customFile,TS_DataMat,TS_CalcTime,TS_Quality,true);
            lastSaveTimer = tic;
        end

    end
end


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%% Finished calculating!!
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
if strcmp(howVocal,'fast')
    if showProgressBar
        BF_ProgressBar('close')
    end
    fprintf(1,'Calculations complete in a total of %s.\n',BF_TheTime(sum(times),1));
else
    fprintf(1,'!! !! !! !! !! !! Calculation completed !! !! !! !! !!\n');
    fprintf(1,'[%s]: Calculations complete in a total of %s.\n',...
                        datestr(now),BF_TheTime(sum(times),1));
end

% Save back to local files (if results were computed):
if any(numCalc_all > 0)
	fprintf(1,'Saving all results to %s...',customFile);
	saveResults(customFile,TS_DataMat,TS_CalcTime,TS_Quality,false);
end

end

%===============================================================================
function saveSeconds = saveResults(customFile,TS_DataMat,TS_CalcTime,TS_Quality,isAutosave)
    % Special/error cells (fatal errors, NaN/Inf/complex outputs) are recorded
    % in TS_Quality; store them in TS_DataMat as NaN so a raw load of the file
    % can't mistake them for computed values. This is the same convention
    % TS_LoadData applies on read, applied here at the single point of truth so
    % that a dataset computed in one call and the same dataset computed in
    % several calls / blocks (batched, resumed, runslice-style) write an
    % identical TS_DataMat. (Older files may still hold the legacy 0 coding;
    % TS_LoadData maps both to NaN.)
    TS_DataMat(~isfinite(TS_DataMat)) = NaN;
    TS_DataMat(TS_Quality > 0) = NaN;

    saveTimer = tic;
    save(customFile,'TS_DataMat','TS_CalcTime','TS_Quality','-append')
    saveSeconds = toc(saveTimer);
    if isAutosave
        fprintf(1,'[%s] Autosaved progress to %s (%s).\n',datestr(now),customFile,BF_TheTime(saveSeconds));
    else
        fprintf(1,' Saved in %s.\n',BF_TheTime(saveSeconds));
    end
end
