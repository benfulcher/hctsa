function TS_compute(doParallel,ts_id_range,op_id_range,computeWhat,customFile,beVocal)
% TS_compute    Computes missing elements of TS_DataMat
%
%---EXAMPLE USAGE:
% TS_compute;
%
%---INPUTS:
% doParallel:  if 1, attempts to use the Parallel Computing Toolbox to run
%               computations in parallel over multiple cores.
% ts_id_range: a custom range of time series IDs to compute (default: [] -- compute all)
% op_id_range: a custom range of operation IDs to compute (default: [] -- compute all)
% computeWhat: whether to compute just missing values ('missing', default), or
% 				ALSO retry results that previously threw an error ('error'), or
% 				ALSO retry any result that previously did not return a good value ('bad')
% customFile: reads in and writes to a custom output file
% beVocal:     if 1, gives additional user feedback about the calculation of
%               each individual operation.
%
%---OUTPUTS:
% Writes output to customFile (HCTSA.mat by default)

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

% --------------------------------------------------------------------------
%% Check inputs and set defaults
% --------------------------------------------------------------------------

% Use Matlab's Parallel Computing toolbox?
if nargin < 1
	doParallel = 0;
end

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
    beVocal = 1; % Write back lots of information to screen
    % prints every piece of code evaluated (nice for error checking)
end

% --------------------------------------------------------------------------
%% Load information from local files
% --------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,customFile] = TS_LoadData(customFile);
fileVarsStruct = whos('-file',customFile);
fileVars = {fileVarsStruct.name};
if ~all(ismember({'TimeSeries','Operations','MasterOperations','TS_DataMat'},fileVars))
	error('\nCannot compute on %s: missing variables.',customFile);
end
load(customFile,'MasterOperations');
if ismember('TS_CalcTime',fileVars)
	load(customFile,'TS_CalcTime');
end
if ismember('TS_Quality',fileVars)
	load(customFile,'TS_Quality');
end

% ------------------------------------------------------------------------------
% Get indices if computing a subset
% ------------------------------------------------------------------------------
allIDs = [TimeSeries.ID];
if isempty(ts_id_range)
    ts_id_range = allIDs;
    tsIndex = 1:length(TimeSeries);
else
    ts_id_range = intersect(ts_id_range,allIDs);
    tsIndex = find(ismember(allIDs,ts_id_range));
    % tsIndex = arrayfun(@(x)find(allIDs==x,1),ts_id_range);
end
allIDs = [Operations.ID];
if isempty(op_id_range)
    op_id_range = allIDs;
    opCompute = ones(1,length(Operations));
else
    op_id_range = intersect(op_id_range,allIDs);
    opCompute = ismember(allIDs,op_id_range);
end

% Definitions
numTimeSeries = length(ts_id_range); % Number of time series
numOps = length(op_id_range); % Number of operations

% Check that some computable range exists
if numTimeSeries==0 || numOps==0
    fprintf(1,'%u time series and %u operations match the ids provided. Exiting.\n',...
						numTimeSeries,numOps);
    return
end

fprintf(1,['Calculation has begun on %s using %u datasets ' ...
                            'and %u operations\n'],datestr(now),numTimeSeries,numOps);


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
numCalc_all = zeros(numTimeSeries,1);

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
    if length([Operations(toCalc).MasterID]) < numCalc
        % Error in the database structure; some operations are missing MasterID assignment
        error('Database structure error: some operations have not been assigned a valid master operation');
    end

	fprintf(1,'\n\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
	fprintf(1,'; ; ; : : : : ; ; ;    %s     ; ; ; : : : ; ; ;\n',datestr(now));
	fprintf(1,'- - - - - - - - - - - Time series %u / %u - - - - - - - - - - -\n',i,numTimeSeries);

    if numCalc > 0 % some to calculate
		try
	        [featureVector,calcTimes,calcQuality] = TS_CalculateFeatureVector(TimeSeries(tsInd),doParallel,Operations(toCalc),MasterOperations,1,beVocal);
		catch
			% skip to the next time series; the entries for this time series in TS_DataMat etc. will remain NaNs
			warning('Calculation for time series %u / %u failed',i,numTimeSeries)
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
    times(i) = toc(bigTimer); clear bigTimer

    if i < numTimeSeries
        fprintf(1,'- - - - - - - -  %u time series remaining - - - - - - - -\n',numTimeSeries-i);
    	fprintf(1,'- - - - - - - -  %s remaining - - - - - - - - -\n', ...
                                        	BF_thetime(((numTimeSeries-i)*mean(times(1:i))),1));
    else % The final time series
        fprintf(1,'- - - - - - - - - - All %u time series calculated! - - - - - - - - - -\n', ...
                                                    numTimeSeries);
    end
    fprintf(1,'********************************************************************\n');

end


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%% Finished calculating!!
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
fprintf(1,['!! !! !! !! !! !! Calculation completed at %s !! !! ' ...
                                                '!! !! !!\n'],datestr(now));
fprintf(1,'Calculations complete in a total of %s.\n',BF_thetime(sum(times),1));

% Save back to local files (if results were computed):
if any(numCalc_all > 0)
	fprintf(1,'Saving all results to %s...',customFile);
	saveTimer = tic;
	save(customFile,'TS_DataMat','TS_CalcTime','TS_Quality','-append')
	fprintf(1,' Saved in %s.\n',BF_thetime(toc(saveTimer)));
	clear saveTimer
end

fprintf(1,'Calculation complete!\n');

end
