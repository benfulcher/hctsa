function [featureVector,calcTimes,calcQuality] = TS_CalculateFeatureVector(tsStruct,doParallel,Operations,MasterOperations,codeSpecial,howVocal)
% TS_CalculateFeatureVector	Compute a feature vector from an input time series
%
%---INPUTS:
% tsStruct, either (i) a vector of time-series data to compute
% 				or (ii) an hctsa-style table of time-series data
% doParallel, (binary) whether to compute the features using parallel processing
% Operations, an hctsa-style table of Operations
% MasterOperations, an hctsa-style table of MasterOperations
% codeSpecial, whether to code special values with quality labels (for mySQL database)
% 				this makes sure the featureVector is all real numbers, and any
%				special values (like NaNs, Infs, etc.) are coded with corresponding
% 				labels in calcQuality.
% 				codeSpecial = false [default]: special values are kept in the
% 							feature vector, and any errors are coded as NaN.
% 				codeSpecial = true: featureVector is all real numbers, and is set to
% 							     zero where any special-valued outputs occur.
% howVocal, {'minimal' or 'full'}: how to give user feedback on the computation.
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
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs
%-------------------------------------------------------------------------------
if isnumeric(tsStruct)
	% Provided data only without metadata:
	tsData = tsStruct;
	tsStruct = struct('Name','Input Timeseries', ...
						'Data',tsData, ...
						'ID',1, ...
						'Length',length(tsData));
elseif istable(tsStruct)
	tsStruct = table2struct(tsStruct);
end

if nargin < 2
    fprintf(1,'Computing features in parallel by default\n');
	doParallel = true;
end

if nargin < 3 || isempty(Operations)
    fprintf(1,'Importing the default set of time-series features\n');
    theINPfile = 'INP_ops.txt';
    Operations = SQL_Add('ops',theINPfile,false,false);
elseif ischar(Operations)
    theINPfile = Operations;
    Operations = SQL_Add('ops',theINPfile,false,false);
end
if isnumeric(Operations)
	error('Provide an input file or a structure array of Operations');
end

if nargin < 4 || isempty(MasterOperations)
	% Use the default library:
	MasterOperations = SQL_Add('mops','INP_mops.txt',false,false);
	% Need to link operations to masters if not already supplied:
	[Operations,MasterOperations] = TS_LinkOperationsWithMasters(Operations,MasterOperations);
end

% Whether to code up special-valued outputs
if nargin < 5
	codeSpecial = false;
end

% Whether to give information out to screen
if nargin < 6
	howVocal = 'minimal';
end

%-------------------------------------------------------------------------------
% Check Statistics toolbox is available (needed throughout hctsa, including for
% z-scoring)
%-------------------------------------------------------------------------------
BF_CheckToolbox('statistics_toolbox');

% ------------------------------------------------------------------------------
%% Open parallel processing worker pool
% ------------------------------------------------------------------------------
if doParallel
    % Check that a parallel worker pool is open (if not, attempt to initiate it):
	doParallel = TS_InitiateParallel(false);
end
% Tell the user about it:
if strcmp(howVocal,'full')
    if doParallel
        fprintf(1,['Computation will be performed across multiple cores' ...
                ' using Matlab''s Parallel Computing Toolbox.\n']);
    else % use single-threaded for loops
        fprintf(1,'Computations will be performed serially without parallelization.\n');
    end
end

% --------------------------------------------------------------------------
%% Basic checking on the data
% --------------------------------------------------------------------------
% (x a N x 1 column vector)
x = tsStruct.Data;
if size(x,2) ~= 1
	if size(x,1) == 1
		fprintf(1,['***** The time series %s is a row vector. Not sure how it snuck',...
								' through the processing cracks, but I need' ...
								' a column vector...\n'],tsStruct.Name);
		fprintf(1,'I''ll transpose it for you for now....\n');
		x = x';
	else
		fprintf(1,'******************************************************************************************\n');
		error('ERROR WITH ''%s'' -- is it multivariate or something weird? Skipping!\n',tsStruct.Name);
	end
end
% Basics checks on data type:
if ~isa(x,'numeric')
    error('The time series provided must be numerical data')
elseif isa(x,'integer')
    error('Methods in hctsa are generally not well-suited to integer-valued time-series data. Convert to double.')
elseif isa(x,'single')
    warning('Your data is provided in single precision. Converting to double for compatibility with hctsa methods.')
    x = double(x);
end
% (x contains no special values)
if ~all(isfinite(x))
	error('ERROR WITH ''%s'' -- contains non-finite values',tsStruct.Name);
end
if all(x==x(1))
	warning('Data are a constant -- there is no information to derive from the time series; many features will fail')
end

% --------------------------------------------------------------------------
%% Display information
% --------------------------------------------------------------------------
numCalc = height(Operations); % Number of features to calculate
if numCalc == 0
	error('Nothing to calculate :-/')
end

if strcmp(howVocal,'full')
    fprintf(1,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
    fprintf(1,'Preparing to calculate %s\nts_id = %u, N = %u samples\nComputing %u features.\n', ...
                    tsStruct.Name,tsStruct.ID,tsStruct.Length,numCalc);
    fprintf(1,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n');
end

% Initialize variables:
featureVector = zeros(numCalc,1); % Output of each operation
calcQuality = zeros(numCalc,1); % Quality of output from each operation
calcTimes = nan(numCalc,1); % Calculation time for each operation

% --------------------------------------------------------------------------
%% Pre-Processing
% --------------------------------------------------------------------------
% x_z is a z-scored transformation of the time series:
x_z = zscore(x);

% So we now have the raw time series x and the z-scored time series x_z.
% Operations take these as inputs.

fullTimer = tic;

% --------------------------------------------------------------------------
%% Evaluate all master operation functions (maybe in parallel)
% --------------------------------------------------------------------------
% Because of parallelization, we have to evaluate all the master functions *first*
% Check through the metrics to determine which master functions are relevant for this run

% Put the output from each Master operation in an element of MasterOutput
MasterOutput = cell(height(MasterOperations),1); % Ouput structures
MasterCalcTime = zeros(height(MasterOperations),1); % Calculation times for each master operation

Master_IDs_calc = unique(Operations.MasterID); % Master_IDs that need to be calculated
Master_ind_calc = arrayfun(@(x)find(MasterOperations.ID==x,1),Master_IDs_calc); % Indicies of MasterOperations that need to be calculated
numMopsToCalc = length(Master_IDs_calc); % Number of master operations to calculate

% Index sliced variables to minimize the communication overhead in the parallel processing
par_MasterOpCodeCalc = MasterOperations.Code(Master_ind_calc); % String array of strings of Code to evaluate
par_mop_ids = MasterOperations.ID(Master_ind_calc); % mop_id for each master operation

% Store in temporary variables for parfor loop then map back later
MasterOutput_tmp = cell(numMopsToCalc,1);
MasterCalcTime_tmp = zeros(numMopsToCalc,1);

% ----
% Evaluate all the master operations
% ----
if strcmp(howVocal,'full')
    fprintf(1,'Evaluating %u master operations...\n',length(Master_IDs_calc));
end
TimeSeries_i_ID = tsStruct.ID; % Make a PARFOR-friendly version of the ID
masterTimer = tic;
if doParallel
    % PARFOR Loop (parallel)
	parfor jj = 1:numMopsToCalc
		[MasterOutput_tmp{jj},MasterCalcTime_tmp(jj)] = ...
			TS_ComputeMasterLoop(x,x_z,par_MasterOpCodeCalc{jj}, ...
				par_mop_ids(jj),numMopsToCalc,howVocal,TimeSeries_i_ID,jj);
	end
else
    % Normal FOR Loop (serial)
    if strcmp(howVocal,'minimal')
        BF_ProgressBar('new');
    end
	for jj = 1:numMopsToCalc
		[MasterOutput_tmp{jj},MasterCalcTime_tmp(jj)] = ...
			TS_ComputeMasterLoop(x,x_z,par_MasterOpCodeCalc{jj}, ...
				par_mop_ids(jj),numMopsToCalc,howVocal,TimeSeries_i_ID,jj);

        if strcmp(howVocal,'minimal')
            BF_ProgressBar(jj/numMopsToCalc);
        end
	end
    if strcmp(howVocal,'minimal')
        BF_ProgressBar('close');
    end
end

% Map from temporary versions to the full versions:
MasterOutput(Master_ind_calc) = MasterOutput_tmp;
MasterCalcTime(Master_ind_calc) = MasterCalcTime_tmp;

if strcmp(howVocal,'full')
    fprintf(1,'%u time-series analysis functions evaluated in %s ///\n\n',...
                        numMopsToCalc,BF_TheTime(toc(masterTimer)));
end
clear('masterTimer')

% --------------------------------------------------------------------------
%% Assign all the results to the corresponding operations
% --------------------------------------------------------------------------
% Set sliced version of matching indicies across the range toCalc
% Indices of MasterOperations corresponding to each Operation (i.e., each index of toCalc)
MasterOp_ind = arrayfun(@(x)find(MasterOperations.ID==x,1),Operations.MasterID);

for jj = 1:numCalc
	try
		[featureVector(jj),calcQuality(jj),calcTimes(jj)] = ...
			TS_ComputeOpLoop(MasterOutput{MasterOp_ind(jj)}, ...
				MasterCalcTime(MasterOp_ind(jj)), ...
				MasterOperations.Label{MasterOp_ind(jj)}, ... % Master label
				Operations.CodeString{jj}); % Code string for each operation to calculate (i.e., in toCalc)
	catch
		fprintf(1,'---Error with %s\n',Operations.CodeString{jj});
		if (MasterOperations.ID(MasterOp_ind(jj)) == 0)
			error(['The operations database is corrupt: there is no link ' ...
					'from ''%s'' to a master code'],Operations.CodeString{jj});
		else
			fprintf(1,'Error retrieving element %s from %s.\n', ...
				Operations.CodeString{jj},MasterOperations.Label{MasterOp_ind(jj)});
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

switch howVocal
case 'full'
    fprintf(1,'********************************************************************\n');
    fprintf(1,'; ; ; : : : : ; ; ; ;   %s    ; ; ; ; : : : ; ; ;\n',datestr(now));
    fprintf(1,'Calculation complete for %s (ts_id = %u, N = %u)\n', ...
                        tsStruct.Name,tsStruct.ID,tsStruct.Length);
    fprintf(1,'%u real-valued outputs, %u errors, %u special-valued outputs stored.\n',...
                        numGood,numErrors,numSpecial);
    fprintf(1,'All %u calculations for this time series took %s.\n',numCalc,BF_TheTime(toc(fullTimer),1));
    fprintf(1,'********************************************************************\n');
case 'minimal'
    fprintf(1,'Calculation complete for %s (ts_id = %u, N = %u) in %s\n', ...
                        tsStruct.Name,tsStruct.ID,tsStruct.Length,BF_TheTime(toc(fullTimer),1));
    fprintf(1,'%u real-valued outputs, %u errors, %u special-valued outputs stored.\n\n',...
                        numGood,numErrors,numSpecial);
case 'fast'
    % (none)
end

end
