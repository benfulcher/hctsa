function [featureVector,calcTimes,calcQuality] = TS_CalculateFeatureVector(tsStruct,doParallel,Operations,MasterOperations,codeSpecial,beVocal)
% TS_CalculateFeatureVector	Compute a feature vector from an input time series
%
%---INPUTS:
% tsStruct, either (i) a vector of time-series data to compute
% 				or (ii) an hctsa-style structure of time-series data
% doParallel, (binary) whether to compute the features using parallel processing
% Operations, an hctsa-style structure array of Operations
% MasterOperations, an hctsa-style structure array of MasterOperations
% codeSpecial, whether to code special values with quality labels (for mySQL database)
% 				this makes sure the featureVector is all real numbers, and any
%				special values (like NaNs, Infs, etc.) are coded with corresponding
% 				labels in calcQuality.
% 				codeSpecial = 0: special values are kept in the feature vector,
% 								and any errors are coded as NaN.
% 				codeSpecial = 1: featureVector is all real numbers, and is set to
% 							     zero where any special-valued outputs occur.
% beVocal, whether to give user feedback on the computation.
%
%---OUTPUTS:
% featureVector, the feature vector obtained by running MasterOperations and
% 					retrieving all features defined in the Operations structure
% 					array on the time series data given in tsStruct.
% calcTimes, corresponding calculation times for each feature in featureVector
% calcQuality, quality labels of each calculation (e.g., coding for NaNs, Infs, etc.)
%
%---USAGE:
% Quickly compute a feature vector of a 500-long random time-series using the
% default operation library using parallel processing:
% >> features = TS_CalculateFeatureVector(randn(500,1));

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

%-------------------------------------------------------------------------------
% Check Inputs
%-------------------------------------------------------------------------------
if isnumeric(tsStruct)
	% Provided data only without metadata:
	tsData = tsStruct;
	tsStruct = struct('Name','Input Timeseries', ...
						'Data',tsData, ...
						'ID',1, ...
						'Length',length(tsData));
end

if nargin < 2
	doParallel = 1;
end

if nargin < 3 || isempty(Operations) || ischar(Operations)
	% Use a default library:
	if nargin >=3 && ischar(Operations)
		theINPfile = Operations;
	else
		theINPfile = 'INP_ops.txt';
	end
	Operations = SQL_add('ops', theINPfile, 0, 0)';
end
if isnumeric(Operations)
	error('Provide an input file or a structure array of Operations');
end

if nargin < 4 || isempty(MasterOperations)
	% Use the default library:
	MasterOperations = SQL_add('mops', 'INP_mops.txt', 0, 0)';
end

% Need to link operations to masters if not already supplied:
if nargin < 4
	[Operations, MasterOperations] = TS_LinkOperationsWithMasters(Operations,MasterOperations);
end

% Whether to code up special-valued outputs
if nargin < 5
	codeSpecial = 0;
end

% Whether to give information out to screen
if nargin < 6
	beVocal = 1;
end

%-------------------------------------------------------------------------------
% Check Statistics toolbox is available (needed throughout hctsa, including for
% zscoring)
%-------------------------------------------------------------------------------
BF_CheckToolbox('statistics_toolbox');

% ------------------------------------------------------------------------------
%% Open parallel processing worker pool
% ------------------------------------------------------------------------------
if doParallel
    % Check that a parallel worker pool is open (if not attempt to initiate it):
	doParallel = TS_InitiateParallel(0);
end
if doParallel
	fprintf(1,['Computation will be performed across multiple cores' ...
			' using Matlab''s Parallel Computing Toolbox.\n']);
else % use single-threaded for loops
	fprintf(1,'Computations will be performed serially without parallelization.\n');
end

% --------------------------------------------------------------------------
%% Basic checking on the data
% --------------------------------------------------------------------------
% (Univariate and [N x 1])
x = tsStruct.Data;
if size(x,2) ~= 1
	if size(x,1) == 1
		fprintf(1,['***** The time series %s is a row vector. Not sure how it snuck through the cracks, but I ' ...
								'need a column vector...\n'],tsStruct.Name);
		fprintf(1,'I''ll transpose it for you for now....\n');
		x = x';
	else
		fprintf(1,'******************************************************************************************\n');
		error('ERROR WITH ''%s'' -- is it multivariate or something weird? Skipping!\n',tsStruct.Name);
	end
end

% --------------------------------------------------------------------------
%% Display information
% --------------------------------------------------------------------------
numCalc = length(Operations); % Number of features to calculate
if numCalc == 0
	error('Nothing to calculate :/')
end

fprintf(1,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
fprintf(1,'Preparing to calculate %s\nts_id = %u, N = %u samples\nComputing %u operations.\n', ...
				tsStruct.Name,tsStruct.ID,tsStruct.Length,numCalc);
fprintf(1,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n');

% Initialize variables:
featureVector = zeros(numCalc,1); % Output of each operation
calcQuality = zeros(numCalc,1); % Quality of output from each operation
calcTimes = ones(numCalc,1)*NaN; % Calculation time for each operation

% --------------------------------------------------------------------------
%% Pre-Processing
% --------------------------------------------------------------------------
% y is a z-scored transformation of the time series:
y = zscore(x);

% So we now have the raw time series x and the z-scored time series y.
% Operations take these as inputs.

fullTimer = tic;

% --------------------------------------------------------------------------
%% Evaluate all master operation functions (maybe in parallel)
% --------------------------------------------------------------------------
% Because of parallelization, we have to evaluate all the master functions *first*
% Check through the metrics to determine which master functions are relevant for this run

% Put the output from each Master operation in an element of MasterOutput
MasterOutput = cell(length(MasterOperations),1); % Ouput structures
MasterCalcTime = zeros(length(MasterOperations),1); % Calculation times for each master operation

Master_IDs_calc = unique([Operations.MasterID]); % Master_IDs that need to be calculated
Master_ind_calc = arrayfun(@(x)find([MasterOperations.ID]==x,1),Master_IDs_calc); % Indicies of MasterOperations that need to be calculated
numMopsToCalc = length(Master_IDs_calc); % Number of master operations to calculate

% Index sliced variables to minimize the communication overhead in the parallel processing
par_MasterOpCodeCalc = {MasterOperations(Master_ind_calc).Code}; % Cell array of strings of Code to evaluate
par_mop_ids = [MasterOperations(Master_ind_calc).ID]; % mop_id for each master operation

fprintf(1,'Evaluating %u master operations...\n',length(Master_IDs_calc));

% Store in temporary variables for parfor loop then map back later
MasterOutput_tmp = cell(numMopsToCalc,1);
MasterCalcTime_tmp = zeros(numMopsToCalc,1);

% ----
% Evaluate all the master operations
% ----
TimeSeries_i_ID = tsStruct.ID; % Make a PARFOR-friendly version of the ID
masterTimer = tic;
if doParallel
	parfor jj = 1:numMopsToCalc % PARFOR Loop
		[MasterOutput_tmp{jj}, MasterCalcTime_tmp(jj)] = ...
					TS_compute_masterloop(x,y,par_MasterOpCodeCalc{jj}, ...
								par_mop_ids(jj),numMopsToCalc,beVocal,TimeSeries_i_ID,jj);
	end
else
	for jj = 1:numMopsToCalc % Normal FOR Loop
		[MasterOutput_tmp{jj}, MasterCalcTime_tmp(jj)] = ...
					TS_compute_masterloop(x,y,par_MasterOpCodeCalc{jj}, ...
								par_mop_ids(jj),numMopsToCalc,beVocal,TimeSeries_i_ID,jj);
	end
end

% Map from temporary versions to the full versions:
MasterOutput(Master_ind_calc) = MasterOutput_tmp;
MasterCalcTime(Master_ind_calc) = MasterCalcTime_tmp;

fprintf(1,'%u master operations evaluated in %s ///\n\n',...
					numMopsToCalc,BF_thetime(toc(masterTimer)));
clear masterTimer

% --------------------------------------------------------------------------
%% Assign all the results to the corresponding operations
% --------------------------------------------------------------------------
% Set sliced version of matching indicies across the range toCalc
% Indices of MasterOperations corresponding to each Operation (i.e., each index of toCalc)
par_OperationMasterInd = arrayfun(@(x)find([MasterOperations.ID]==x,1),[Operations.MasterID]);
par_MasterOperationsLabel = {MasterOperations.Label}; % Master labels
par_OperationCodeString = {Operations.CodeString}; % Code string for each operation to calculate (i.e., in toCalc)

if doParallel
	parfor jj = 1:numCalc
		[featureVector(jj), calcQuality(jj), calcTimes(jj)] = TS_compute_oploop(MasterOutput{par_OperationMasterInd(jj)}, ...
									   MasterCalcTime(par_OperationMasterInd(jj)), ...
									   par_MasterOperationsLabel{par_OperationMasterInd(jj)}, ...
									   par_OperationCodeString{jj});
	end
else
	for jj = 1:numCalc
		try
			[featureVector(jj), calcQuality(jj), calcTimes(jj)] = TS_compute_oploop(MasterOutput{par_OperationMasterInd(jj)}, ...
									   MasterCalcTime(par_OperationMasterInd(jj)), ...
									   par_MasterOperationsLabel{par_OperationMasterInd(jj)}, ...
									   par_OperationCodeString{jj});
		catch
			fprintf(1,'---Error with %s\n',par_OperationCodeString{jj});
			if (MasterOperations(par_OperationMasterInd(jj)).ID == 0)
				error(['The operations database is corrupt: there is no link ' ...
						'from ''%s'' to a master code'], par_OperationCodeString{jj});
			else
				fprintf(1,'Error retrieving element %s from %s.\n', ...
					par_OperationCodeString{jj}, par_MasterOperationsLabel{par_OperationMasterInd(jj)});
			end
		end
	end
end

% --------------------------------------------------------------------------
%% Code special values:
% --------------------------------------------------------------------------
% (*) Errorless calculation: calcQuality = 0, featureVector = <real number>

% (*) Output = NaN: calcQuality = 2, set featureVector = 0
RR = isnan(featureVector); % NaN
if any(RR)
	calcQuality(RR) = 2;
	if codeSpecial
		featureVector(RR) = 0;
	end
end

% (*) Output = Inf: calcQuality = 3, set featureVector = 0
RR = (isinf(featureVector) & featureVector > 0); % Inf
if any(RR)
	calcQuality(RR) = 3;
	if codeSpecial
		featureVector(RR) = 0;
	end
end

% (*) Output = -Inf: calcQuality = 4, set featureVector = 0
RR = (isinf(featureVector) & featureVector < 0);
if any(RR)
	calcQuality(RR) = 4;
	if codeSpecial
		featureVector(RR) = 0;
	end
end

% (*) Output is a complex number: calcQuality = 5, set featureVector = 0
RR = (imag(featureVector)~=0);
if any(RR)
	calcQuality(RR) = 5;
	if codeSpecial
		featureVector(RR) = 0;
	end
end

% (*) Fatal error: calcQuality = 1, featureVector = 0
% 		(this is done already in the code above)
if ~codeSpecial
	% If you want special values in your feature vector (~codeSpecial), then
	% we need to put errors back into the feature vector as NaNs:
	RR = (calcQuality==1);
	featureVector(RR) = NaN;
end

% --------------------------------------------------------------------------
%% Calculation complete: display information about this time series calculation
% --------------------------------------------------------------------------

% Calculate statistics for writing to file/screen
% The number of calculated operations that returned real outputs without errors, numGood:
numGood = sum(calcQuality == 0);
% The number of errors encountered, numErrors:
numErrors = sum(calcQuality == 1);
% The number of other special outputs, numSpecial:
numSpecial = sum(calcQuality > 1);

fprintf(1,'********************************************************************\n');
fprintf(1,'; ; ; : : : : ; ; ; ;   %s    ; ; ; ; : : : ; ; ;\n',datestr(now));
fprintf(1,'Calculation complete for %s (ts_id = %u, N = %u)\n', ...
					tsStruct.Name,tsStruct.ID,tsStruct.Length);
fprintf(1,'%u real-valued outputs, %u errors, %u special-valued outputs stored.\n',...
					numGood,numErrors,numSpecial);
fprintf(1,'All %u calculations for this time series took %s.\n',numCalc,BF_thetime(toc(fullTimer),1));
fprintf(1,'********************************************************************\n');

end
