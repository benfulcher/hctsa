% --------------------------------------------------------------------------
% TSQ_brawn
% --------------------------------------------------------------------------
% 
% This function fills in the missing elements of TS_DataMat, from HCTSA_loc.mat
% (retrieved from the database using TSQ_prepared).
% The function systematically calculates these missing elements (in parallel
% over operations for each time series if specified).
% 
%---INPUTS:
% doLog:      if 1 (0 by default) writes to a log file.
% doParallel: if 1, attempts to use the Parallel Computing Toolbox to run
%             computations in parallel over multiple cores.
% bevocal:    if 1, gives additional user feedback.
% 
%---OUTPUTS:
% Writes output to file, over HCTSA_loc.mat
%
% ------------------------------------------------------------------------------
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

function TSQ_brawn(doLog,doParallel,bevocal)
    
% --------------------------------------------------------------------------
%% Check inputs and set defaults
% --------------------------------------------------------------------------

% Log to file
if nargin < 1
    % By default, do not log to file, write to screen (if bevocal)
	doLog = 0;
end
if doLog
	fn = sprintf('HCTSA_brawn_%s.log',datestr(now,30));
	fid = fopen(fn,'w','n');
	fprintf(1,'Calculation details will be logged to %s\n',fn);
else
    % Write output to screen rather than .log file
    fid = 1;
end

% Use Matlab's Parallel Computing toolbox?
if nargin < 2
	doParallel = 0;
end
if doParallel
	fprintf(fid,'Computation will be performed across multiple cores using Matlab''s Parallel Computing Toolbox\n')
else % use single-threaded for loops
	fprintf(fid,'Computations will be performed serially without parallelization\n')
end

% Be vocal?
if nargin < 3
    bevocal = 1; % Write back lots of information to screen
    % prints every piece of code evaluated (nice for error checking)
end


% --------------------------------------------------------------------------
%% Load information from local files
% --------------------------------------------------------------------------
fprintf(fid,'Reading in data from HCTSA_loc.mat...');
load('HCTSA_loc.mat')
fprintf(fid,' Done.\n');

% Definitions
NumTimeSeries = length(TimeSeries); % Number of time series
NumOps = length(Operations); % Number of operations

fprintf(fid,['Calculation has begun on %s using %u datasets ' ...
                            'and %u operations\n'],datestr(now),NumTimeSeries,NumOps);

% ------------------------------------------------------------------------------
%% Open parallel processing worker pool
% ------------------------------------------------------------------------------
if doParallel
    % first check that the user can use the Parallel Computing Toolbox:
    heyo = which('matlabpool');
    if isempty(heyo)
        fprintf(1,['Parallel Computing Toolbox not found -- ' ...
                        'cannot perform computations across multiple cores\n'])
        doParallel = 0;
    else
        if (matlabpool('size') == 0)
        	matlabpool open;
            fprintf(fid,['Matlab parallel processing pool opened with %u ' ...
                                    'and ready to go'],matlabpool('size'))
        else
            fprintf(fid,['Matlab parallel processing pool already open. ' ...
                                        'Size: %u\n'],matlabpool('size'))
        end
    end
end


% Times stores the time taken for each time series to have its operations
% calculated (for determining time remaining)
times = zeros(NumTimeSeries,1); 
lastSavedTime = 0; % Last saved time

% --------------------------------------------------------------------------
%% Computation
% --------------------------------------------------------------------------
for i = 1:NumTimeSeries
	bigtimer = tic;

	qq = TS_Quality(i,:); % The calculation states of any existing results for the current time series, a line of TS_Quality
					   	  % NaN indicates a value never before calculated, 1 indicates fatal error before (try again now)
    tcal = (isnan(qq) | qq == 1); % Operations awaiting calculation
    ncal = sum(tcal); % Number of operations to evaluate
    
    % Check that all operations have a master ID attached:
    if length([Operations(tcal).MasterID]) < ncal
        % Error in the database structure; some 
    end
    
    if ncal > 0 % some to calculate
        ffi = zeros(ncal,1); % Output of each operation
		qqi = zeros(ncal,1); % Quality of output from each operation
		cti = ones(ncal,1)*NaN; % Calculation time for each operation
        
        % Load the time-series data as x
        x = TimeSeries(i).Data;
        
        % --------------------------------------------------------------------------
		%% Basic checking on x
        % --------------------------------------------------------------------------
		% Univariate and [N x 1]
		if size(x,2) ~= 1
			if size(x,1) == 1
                fprintf(fid,['***** The time series %s is a row vector. Not sure how it snuck through the cracks, but I ' ...
                                            'need a column vector...\n'],TimeSeries(i).FileName);
				fprintf(fid,'I''ll transpose it for you for now....\n');
				x = x';
			else
				fprintf(fid,'******************************************************************************************\n')
                fprintf(fid,['MASSIVE ERROR WITH THIS TIME SERIES!!!: %s -- is it multivariate or something weird???.' ...
                                                                    ' Skipping it!!\n'], TimeSeries(i).FileName);
				fprintf(fid,'******************************************************************************************\n');
				continue % skip to the next time series; the entries for this time series in TS_DataMat etc. will remain NaNs
            end
		end

        % --------------------------------------------------------------------------
        %% Pre-Processing
        % --------------------------------------------------------------------------
		% y is a z-scored transformation of the time series
        % z-score without using a Statistics Toolbox license (i.e., the 'zscore' function):
		y = BF_zscore(x);

		% So we now have the raw time series x and the z-scored time series y.
		% Operations take these as inputs.

        % --------------------------------------------------------------------------
		%% Display information
        % --------------------------------------------------------------------------
		WhichTimeSeries = which(TimeSeries(i).FileName);
		fprintf(fid,'\n\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n')
	    fprintf(fid,'; ; ; : : : : ; ; ; ;    %s     ; ; ; ; : : : ; ; ;\n',datestr(now))
	    fprintf(fid,'- - - - - - - -  Loaded time series %u / %u - - - - - - - -\n',i,NumTimeSeries)
		fprintf(fid,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n')
		fprintf(fid,'%s\n',WhichTimeSeries)
		fprintf(fid,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n')
		fprintf(fid,'Preparing to calculate %s\nts_id = %u, N = %u samples\nComputing %u / %u operations.\n', ...
                            		TimeSeries(i).FileName,TimeSeries(i).ID,TimeSeries(i).Length,ncal,NumOps)
	    fprintf(fid,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')

        % --------------------------------------------------------------------------
		%% Evaluate all master operation functions (maybe in parallel)
        % --------------------------------------------------------------------------
		% Because of parallelization, we have to evaluate all the master functions *first*
		% Check through the metrics to determine which master functions are relevant for this run

		% Put the output from each Master operation in an element of MasterOutput
		MasterOutput = cell(length(MasterOperations),1); % Ouput structures
		MasterCalcTime = zeros(length(MasterOperations),1); % Calculation times for each master operation
		
		Master_IDs_calc = unique([Operations(tcal).MasterID]); % Master_IDs that need to be calculated
        Master_ind_calc = arrayfun(@(x)find([MasterOperations.ID]==x,1),Master_IDs_calc); % Indicies of MasterOperations that need to be calculated
		nMtocalc = length(Master_IDs_calc); % Number of master operations to calculate
        
        % Index sliced variables to minimize the communication overhead in the parallel processing
        par_MasterOperationCodeCalculate = {MasterOperations(Master_ind_calc).Code}; % Cell array of strings of Code to evaluate
        % par_OperationMasterID = [Operations(tcal).MasterID]; % Master_IDs corresponding to each Operation
        
		fprintf(fid,'Evaluating %u master operations...\n',length(Master_IDs_calc));
		
	    % Store in temporary variables for parfor loop then map back later
        MasterOutput_tmp = cell(nMtocalc,1);
        MasterCalcTime_tmp = zeros(nMtocalc,1);
		
		% Evaluate all the master operations
        TimeSeries_i_ID = TimeSeries(i).ID; % Make a PARFOR-friendly version of the ID
        master_timer = tic;
		if doParallel
            parfor jj = 1:nMtocalc % PARFOR Loop
                [MasterOutput_tmp{jj}, MasterCalcTime_tmp(jj)] = ...
                                TSQ_brawn_masterloop(x,y,par_MasterOperationCodeCalculate{jj}, ...
                                                        fid,bevocal,TimeSeries_i_ID);
            end
        else
            for jj = 1:nMtocalc % Normal FOR Loop
                [MasterOutput_tmp{jj}, MasterCalcTime_tmp(jj)] = ...
                                TSQ_brawn_masterloop(x,y,par_MasterOperationCodeCalculate{jj}, ...
                                                        fid,bevocal,TimeSeries_i_ID);
            end
		end
		
        % Map back from temporary versions to the full versions
        MasterOutput(Master_ind_calc) = MasterOutput_tmp;
        MasterCalcTime(Master_ind_calc) = MasterCalcTime_tmp;
		
		fprintf(fid,'%u master operations evaluated in %s ///\n\n',nMtocalc,BF_thetime(toc(master_timer)));
        clear master_timer
        
        
        % Set sliced version of matching indicies across the range tcal
        % Indices of MasterOperations corresponding to each Operation (i.e., each index of tcal)
        par_OperationMasterInd = arrayfun(@(x)find([MasterOperations.ID]==x,1),[Operations(tcal).MasterID]);
        par_MasterOperationsLabel = {MasterOperations.Label}; % Master labels
        par_OperationCodeString = {Operations(tcal).CodeString}; % Code string for each operation to calculate (i.e., in tcal)
        
        % --------------------------------------------------------------------------
		%% Assign all the results to the corresponding operations
        % --------------------------------------------------------------------------
		if doParallel
	        parfor jj = 1:ncal
                [ffi(jj), qqi(jj), cti(jj)] = TSQ_brawn_oploop(MasterOutput{par_OperationMasterInd(jj)}, ...
                                                   MasterCalcTime(par_OperationMasterInd(jj)), ...
                                                   par_MasterOperationsLabel{par_OperationMasterInd(jj)}, ...
                                                   par_OperationCodeString{jj},fid);
            end
		else
            for jj = 1:ncal
                try
                    [ffi(jj), qqi(jj), cti(jj)] = TSQ_brawn_oploop(MasterOutput{par_OperationMasterInd(jj)}, ...
                                                       MasterCalcTime(par_OperationMasterInd(jj)), ...
                                                       par_MasterOperationsLabel{par_OperationMasterInd(jj)}, ...
                                                       par_OperationCodeString{jj},fid);
                catch
                    fprintf(1,'---Error with %s\n',par_OperationCodeString{jj});
                    keyboard
                    if (MasterOperations(par_OperationMasterInd(jj)).ID == 0)
                        error(['The operations database is corrupt: there is no link ' ...
                                        'from ''%s'' to a master code'], ...
                                        par_OperationCodeString{jj});
                    else
                        fprintf(1,'Error retrieving element %s from %s\n',par_OperationCodeString{jj}, ...
                                        par_MasterOperationsLabel{par_OperationMasterInd(jj)})
                        RA_keyboard
                    end
                end
            end
        end
        
        % --------------------------------------------------------------------------
		%% Code special values:
        % --------------------------------------------------------------------------
		% (*) Errorless calculation: q = 0, Output = <real number>
		% (*) Fatal error: q = 1, Output = 0; (this is done already in the code above)

		% (*) Output = NaN: q = 2, Output = 0
		RR = isnan(ffi); % NaN
		if any(RR)
			qqi(RR) = 2; ffi(RR) = 0;
		end

		% (*) Output = Inf: q = 3, Output = 0
		RR = (isinf(ffi) & ffi > 0); % Inf
		if any(RR)
			qqi(RR) = 3; ffi(RR) = 0;
		end
		
        % (*) Output = -Inf: q = 4, Output = 0
		RR = (isinf(ffi) & ffi < 0);
		if any(RR)
			qqi(RR) = 4; ffi(RR) = 0;
		end
        
		% (*) Output is a complex number: q = 5, Output = 0
		RR = (imag(ffi)~=0);
		if any(RR)
			qqi(RR) = 5; ffi(RR) = 0;
		end

        % ------------------------------------------------------------------------------
		%% Store the calculated information back to local matrices
        % ------------------------------------------------------------------------------
        TS_DataMat(i,tcal) = ffi; % store outputs in TS_DataMat
		TS_CalcTime(i,tcal) = cti; % store calculation times in TS_CalcTime
		TS_Quality(i,tcal) = qqi; % store quality labels in TS_Quality
		% NB: the calculation time assigned for individual operations is the total calculation
		% time taken to evaluate the master code.

        % Calculate statistics for writing to file/screen
        % The number of calculated operations that returned real outputs without errors, numGood:
		numGood = sum(qqi == 0);
        % The number of fatal errors encountered, numErrors:
		numErrors = sum(qqi == 1);
        % The number of other special outputs, numSpecial:
		numSpecial = sum(qqi > 1);
    end
    
    % The time taken to calculate (or not, if ncal = 0) all operations for this time series:
    times(i) = toc(bigtimer); clear bigtimer


    % --------------------------------------------------------------------------
    %% Calculation complete: display information about this time series calculation
    % --------------------------------------------------------------------------
	fprintf(fid,'********************************************************************\n')
    fprintf(fid,'; ; ; : : : : ; ; ; ;   %s    ; ; ; ; : : : ; ; ;\n',datestr(now))
    fprintf(fid,'oOoOo Calculation complete for %s oOoOo\nts_id = %u, N = %u\n', ...
                            TimeSeries(i).FileName,TimeSeries(i).ID,TimeSeries(i).Length);
    if ncal > 0 % Some amount of calculation was performed
	    fprintf(fid,'%u real-valued outputs, %u errors, %u special-valued outputs stored. [%u / %u]\n',...
	     					numGood,numErrors,numSpecial,ncal,NumOps);
		fprintf(fid,'All calculations for this time series took %s.\n',BF_thetime(times(i),1));
    else
    	fprintf(fid,'Nothing calculated! All %u operations already complete!!  0O0O0O0O0O0\n',NumOps);
    end
    if i < NumTimeSeries
        fprintf(fid,'- - - - - - - -  %u time series remaining - - - - - - - -\n',NumTimeSeries-i);
    	fprintf(fid,'- - - - - - - -  %s remaining - - - - - - - - -\n', ...
                                        	BF_thetime(((NumTimeSeries-i)*mean(times(1:i))),1));
    else % The final time series
        fprintf(fid,'- - - - - - - - - - All %u time series calculated! - - - - - - - - - -\n', ...
                                                    NumTimeSeries);
    end
    fprintf(fid,'********************************************************************\n');


 	%     % --------------------------------------------------------------------------
 	% %% Save to HCTSA_loc every 10 minutes in case of disaster
 	%     % --------------------------------------------------------------------------
 	%     if (sum(times(1:i)) - lastSavedTime) > 60*10
 	% 		% It's been more than 10 mins since last save
 	%         fprintf(fid,'Not finished calculations yet, but saving progress so far to file...')
 	%
 	%         save('HCTSA_loc.mat','TS_DataMat','TS_CalcTime','TS_Quality','TimeSeries', ...
 	%                                     'Operations','MasterOperations','-v7.3')
 	%
 	%         fprintf(fid,' All saved.\n')
 	%         lastSavedTime = sum(times(1:i)); % The last saved time
 	%     end
end


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%% Finished calculating!!
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
fprintf(fid,['!! !! !! !! !! !! !!  Calculation completed at %s !! !! !! ' ...
                                                '!! !! !!\n'],datestr(now))
fprintf(fid,'Calculations took a total of %s.\n',BF_thetime(sum(times),1))

% Save the local files for subsequent upload to the mySQL database
fprintf(1,'Saving all results to HCTSA_loc.mat...')
save('HCTSA_loc.mat','TS_DataMat','TS_CalcTime','TS_Quality','TimeSeries', ...
                                'Operations','MasterOperations','-v7.3')
fprintf(fid,' All saved.\n')

% Close the .log file:
if doLog
	fclose(fid);
end

fprintf(1,['Calculation complete: you can now run TSQ_agglomerate to upload' ...
                                ' the results to a mySQL database\n'])

end