function TSQ_brawn(tolog,toparallel,dbname,agglomerate)
% This is a routine to fill in the blank parts of TS_loc, as prepared by TSQ_prepared
% It reads in the output files produced by TS_prepared and then goes through and systematically
% calculates the missing bits (in parallel over the metrics for each time series)

%% INPUT (function argument):
% (()) tolog: if 1 (0 by default) writes to a log file
%% INPUTS (files):
% The output of TSQ_prepared: TS_loc.mat, TS_guides.mat
%% OUTPUTS (to file):

% 12/11/2000 ~ Ben Fulcher ~ Added tolog to functionality
% 13/11/2009 ~ Ben Fulcher ~ Added noparallel to functionality, mainly for testing purposes
% 							, but could also be for single-core machines
% 30/11/2009 ~ Ben Fulcher ~ Completely rehauled for mySQL usage
% 11/1/2010 ~ Ben Fulcher ~ Quickly added dbc argument -- TSQ_agglomerate writes back to this after calculation
% 12/4/2011 ~ Ben Fulcher ~ added agglomerate argument -- can specify 0 to
%                               stop autoagglomeration

%%% FOREPLAY

%% Check valid inputs / set defaults
% Log to file
if nargin < 1
	tolog = 0; % by default, do not log to file
end
if tolog
	fn = sprintf('TS_brawn_%s.log',datestr(now,30));
	fid = fopen(fn,'w','n');
	fprintf(1,'Calculation details will be logged to %s\n',fn);
else
    fid = 1; % write output to screen rather than .log file
end

% Parallel processing
if nargin < 2
	toparallel = 0;
end
if toparallel
	fprintf(fid,'Computation will be performed across multiple cores using Matlab''s Parallel Computing Toolbox\n')
else % use single-threaded for loops
	fprintf(fid,'Computations will be performed serially without parallelization\n')
end

% mySQL database
if nargin < 3
	dbname = ''; % use default database by default
end

% Write back results when finished?
if nargin < 4
    agglomerate = 1; % yes, agglomerate by default
end

%% Read in information from local files
fprintf(fid,'Reading in data and guides from file...');
load TS_loc.mat TS_loc
fprintf(fid,' TS_loc, '); % data
load TS_loc_ct.mat TS_loc_ct
fprintf(fid,' TS_loc_ct, '); % calculation times
load TS_loc_q.mat TS_loc_q
fprintf(fid,' TS_loc_q,'); % quality codes
load TS_loc_guides.mat ts_ids_keep tsf tskw tsl m_ids_keep mcode mlab mkw mlink Mmid Mmlab Mmcode
fprintf(fid,'TS_loc_guides. All loaded.\n'); % guides

%% Definitions
nts = length(ts_ids_keep); % number of time series
nm = length(m_ids_keep); % number of metrics
nMm = length(Mmlab); % number of master functions

%%% Let's get going
fprintf(fid,'Calculation has begun on %s using %u datasets and %u operations\n',datestr(now),nts,nm);

%% Open parallel processing worker pool
if toparallel
	if matlabpool('size')==0
		matlabpool open;
		fprintf(fid,'MATLAB parallel processing pool opened with %u and ready to go',matlabpool('size'))
	else
		fprintf(fid,'MATLAB parallel processing pool already open. Size: %u\n',matlabpool('size'))
	end
end


times = zeros(nts,1); % stores the time taken for each time series to have its metrics calculated (for determining time remaining)
lst = 0; % last saved time
bevocal = 0; % print every piece of code evaluated (nice for error checking)

for i = 1:nts
	bigtimer = tic;

	qq = TS_loc_q(i,:); % the calculation states of any existing results for the current time series, a line of TS_loc_q
					   	% NaN indicates a value never before calculated, 1 indicates fatal error before -- try again this time
    tcal = (isnan(qq) | qq == 1); % time series awaiting calculation
    ncal = sum(tcal); % number of time series to calculate
    
    if ncal > 0 % some to calculate
        ffi = zeros(ncal,1); % output of metrics
		qqi = zeros(ncal,1); % label of quality of outputs
		cti = ones(ncal,1)*NaN; % calculation times
        x = dlmread(tsf{i}); % this is the raw time series from the file
        
		% write to log
		whichtsf = which(tsf{i});
		if tolog
			fprintf(fid,'\n\n%s\n','=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=');
		    fprintf(fid,'%s\n',['; ; ; : : : : ; ; ; ;    ' datestr(now) '     ; ; ; ; : : : ; ; ;']);
			fprintf(fid,'%s\n',['Loaded ' whichtsf '  (' num2str(tsl(i)) ' samples)']);
		    fprintf(fid,'%s\n',['    - - - - - - - -  ' num2str(i) '  /  ' num2str(nts) '   - - - - - - - -']);
			fprintf(fid,'%s\n',['    - - - - - - ' benrighttime(((nts-i)*mean(times(1:i)))) ' remaining - - - - - - -']);
		    fprintf(fid,'%s\n\n','=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=');
		end

		%% Basic checking on x
		% univariate and Nx1
		if size(x,2)~=1
			if size(x,1)==1
				fprintf(1,'***** The time series %s is a row vector. Not sure how it snuck through the gates, but I need a column vector...\n',tsf{i});
				fprintf(1,'I''ll transpose it for you for now....\n');
				if tolog, fprintf(fid,'*** %s a 1xN vector... Transposed.',tsf{i}); end
				x = x';
			else
				fprintf(1,'******************************************************************************************\n')
				fprintf(1,'MASSIVE ERROR WITH THIS TIME SERIES!!!: %s -- is it multivariate or something weird???. Skipping it!!\n', tsf{i});
				fprintf(1,'******************************************************************************************\n');
				if tolog, fprintf(fid,'*** MASSIVE ERROR: A MULTIVARIATE(?) TIME SERIES... SKIPPING...\n'); end
				continue % skip to the next time series; the entries for this time series in TS_loc etc. will remain NaNs
            end
		end

        %% Pre-Processing
		% Keep a z-scored version of the time series as y
		y = (x-mean(x))/std(x); % z-scoring without using a Statistics Toolbox license (the zscore function)

		% So we now have the raw time series x and the z-scored time series y. Operations take these as inputs.
        % index sliced variables to minimize the communication overhead in the parallel processing
        partsfi  = tsf{i}; parmcode = mcode(tcal); parmlab = mlab(tcal); parmlink = mlink(tcal);

		fprintf(1,'Calculating %s (ts_id = %u, N = %u) [%u / %u]\n',partsfi,ts_ids_keep(i),tsl(i),ncal,nm)
		if tolog
			fprintf(fid,'Calculating %u (/%u) operations\n',ncal,nm);
		end
		
		%% Evaluate master functions in parallel
		% Because of parallelization, we have to evaluate all the master functions *first*
		% Check through the metrics to determine which master functions are relevant for this run
		ip = find(parmlink > 0); % a range of tcal -- that point to a master function
		if ~isempty(ip) % there are some pointers
			% For each implicated master function -- put the structure in this element
			
			Moutput = cell(nMm,1); % Ouput structures
			Mcts = zeros(nMm,1); % Calculation times for the master structures
			
			Mitocalc = unique(parmlink(parmlink > 0)); % indicies of masters (Mmlab/Mmcode) that need to be called (i.e., they'll be called later)
			nMtocalc = length(Mitocalc); % number of implicated master functions

			fprintf(fid,'Evaluating %u master operations first...\n',nMtocalc);
			
		    % store in temporary variables then map back
            MoutputPAR = cell(nMtocalc,1);
            MctsPAR = zeros(nMtocalc,1);
			
			% Evaluate the master functions
			if toparallel
	            parfor j = 1:nMtocalc
					jj = Mitocalc(j); % the index of the master array to evaluate
					if bevocal, fprintf(fid,'%s\n',Mmcode{jj}); end % display for error checking
					try
						mastertimer = tic;
						% disp(Mmcode{jj}) % FOR DEBUGGING -- SHOW CODE BEFORE EVALUATING TO DETERMINE PROBLEMS
						MoutputPAR{j} = pareval(x,y,Mmcode{jj});
						% for not-applicable/'real NaN', outputs a NaN, otherwise should output a
						% structure with components to be called below by pointer metrics.
						MctsPAR(j) = toc(mastertimer);
                    catch emsg
						fprintf(fid,'**** ERROR EVALUATING MASTER FILE: %s (%s)\n',Mmlab{jj},Mmcode{jj});
                        fprintf(fid,'%s\n',emsg.message)
						% masterdat{j} didn't evaluate properly -- remains an empty cell entry.
					end
                end
			else % just a single-thread for loop
	            for j = 1:nMtocalc
					jj = Mitocalc(j); % the index of the master array to evaluate
					if bevocal, fprintf(fid,'%s\n',Mmcode{jj}); end % for error checking
					try
						mastertimer = tic;
						MoutputPAR{j} = pareval(x,y,Mmcode{jj});
						% for not-applicable/'real NaN', outputs a NaN, otherwise should output a
						% structure with components to be called below by pointer metrics.
						MctsPAR(j) = toc(mastertimer);
                    catch emsg
						fprintf(fid,'**** ERROR EVALUATING MASTER FILE: %s (%s)\n',Mmlab{jj},Mmcode{jj});
                        fprintf(fid,'%s\n',emsg.message)
                        % masterdat{j} didn't evaluate properly -- remains an empty cell entry.
					end
                end
			end
			
			% map back from PAR versions to real versions
            Moutput(Mitocalc) = MoutputPAR;
            Mcts(Mitocalc) = MctsPAR;
			
			fprintf(fid,'Master operations evaluated.\n');
		else
			% No master metrics need to be calculed.
			Moutput = {}; Mcts = {}; % This initiaition is necessary for the next parfor loop
			fprintf(fid,'No master operations.\n');
		end


		%% GO GO GO!!!
		if toparallel
	        parfor j = 1:ncal
				if parmlink(j) > 0 % pointer to a master function
					try
						% retrive from master structure:
						if ~isstruct(Moutput{parmlink(j)}) && isnan(Moutput{parmlink(j)});
							% All master function outputs are to be set to real NaNs (unsuitable)
							ffi(j) = NaN;
		                    cti(j) = Mcts(parmlink(j));
						else % (normal -- retrieve required element from master structure)
                            [~,thest] = strtok(parmcode{j},'.'); thest = thest(2:end); % remove the '.'
		                    ffi(j) = parevalM(Moutput{parmlink(j)},['themasterdat.' thest]);
							qqi(j) = 0; % good, real-valued output
		                    cti(j) = Mcts(parmlink(j));
						end
                    catch emsg
						fprintf(fid,'Error evaluating link to Master structure %s by %s\n',Mmlab{parmlink(j)},parmcode{j});
                        fprintf(fid,'%s\n',emsg.message)
                        ffi(j) = 0;
						qqi(j) = 1; % fatal error code
					end
				else % A regular, single-output operation
                    if bevocal, fprintf(fid,'%s\n',parmcode{j}); end % for error checking
		            try
						operationtimer = tic;
						ffi(j) = pareval(x,y,parmcode{j});
						qqi(j) = 0;
						cti(j) = toc(operationtimer);
		            catch
		                fprintf(fid,'Fatal error %s || %s\n',partsfi,parmcode{j});
						ffi(j) = 0; qqi(j) = 1; % fatal error code = 1
		            end
				end
	        end
		else
	        for j = 1:ncal
				if parmlink(j) > 0 % pointer to a master function
					try
						% retrive from master structure:
						if ~isstruct(Moutput{parmlink(j)}) && isnan(Moutput{parmlink(j)});
							% All master function outputs are to be set to real NaNs (unsuitable)
							ffi(j) = NaN;
		                    cti(j) = Mcts(parmlink(j));
						else % (normal -- retrieve required element from master structure)
                            [~,thest] = strtok(parmcode{j},'.'); thest = thest(2:end); % remove the '.'
		                    ffi(j) = parevalM(Moutput{parmlink(j)},['themasterdat.' thest]);
							qqi(j) = 0; % good, real-valued output
		                    cti(j) = Mcts(parmlink(j));
						end
	                catch
						fprintf(fid,'Error evaluating link to Master structure %s by %s\n',Mmlab{parmlink(j)},parmcode{j});
                        fprintf(fid,'%s\n',emsg.message)
						ffi(j) = 0; qqi(j) = 1; % Fatal error code: 1
					end
				else % A regular, single-output operation
                    if bevocal, fprintf(fid,'%s\n',parmcode{j}); end % for error checking
		            try
						operationtimer = tic;
						ffi(j) = pareval(x,y,parmcode{j});
						qqi(j) = 0;
						cti(j) = toc(operationtimer);
		            catch
		                fprintf(fid,'Fatal error %s || %s\n',partsfi,parmcode{j});
						ffi(j) = 0; qqi(j) = 1; % fatal error code = 1
		            end
				end
	        end
		end
        
		%% Error coding -- qqi

		% (*) Errorless calculation: q = 0, output = <real number>
		% (*) Fatal error: q = 1, output = 0; (this is done already in the code above)

		% (*) output = NaN: q = 2, output = 0
		RR = isnan(ffi); % NaN
		if any(RR)
			qqi(RR) = 2; ffi(RR) = 0;
		end

		% (*) output = Inf: q = 3, output = 0
		RR = (isinf(ffi) & ffi > 0); % Inf
		if any(RR)
			qqi(RR) = 3; ffi(RR) = 0;
		end
		
        % (*) output = -Inf: q = 4, output = 0
		RR = (isinf(ffi) & ffi < 0);
		if any(RR)
			qqi(RR) = 4; ffi(RR) = 0;
		end
        
		% (*) output is a complex number: q = 5, output = 0
		RR = (imag(ffi)~=0);
		if any(RR)
			qqi(RR) = 5; ffi(RR) = 0;
		end

		
		% Store the calculated information back to local matrices
        TS_loc(i,tcal) = ffi; % store outputs in TS_loc
		TS_loc_ct(i,tcal) = cti; % store calculation times in TS_loc_ct
		TS_loc_q(i,tcal) = qqi; % store quality labels in TS_loc_q
		% Note that the calculation time for pointers is the total calculation time for the full master structure.		

        % Calculate statistics for writing to file/screen
		ngood = sum(qqi==0); % the number of calculated metrics that returned real outputs without errors
		nerror = sum(qqi==1); % number of fatal errors encountered
		nother = sum(qqi>1); % number of other special outputs
    end
    
    times(i) = toc(bigtimer); % the time taken to calculate (or not, if ncal = 0) operations for this time series

    % Print information about this time series calculation
	fprintf(fid,'****************************************************************\n')
    fprintf(fid,'; ; ; : : : : ; ; ; ;    %s     ; ; ; ; : : : ; ; ;\n',datestr(now))
    fprintf(fid,'oOoOoOoOoOo  %s (%u)  oOoOoOoOoOoOo\n',tsf{i},tsl(i));
    if ncal > 0 % Some amount of calculation was performed
	    fprintf(fid,'%u real-valued outputs, %u errors, %u other outputs stored. [%u / %u]\n',...
	     					ngood,nerror,nother,ncal,nm);
		fprintf(fid,'Calculation for this data set took %s.\n',benrighttime(times(i)));
    else
    	fprintf(fid,'Nothing calculated! All %u operations already complete!!  0O0O0O0O0O0',nm);
    end
    fprintf(fid,'    - - - - - - - -  %u / %u   - - - - - - - -\n',i,nts);
	fprintf(fid,'    - - - - - - - -  %s remaining - - - - - - - - -\n',benrighttime(((nts-i)*mean(times(1:i)))));
    fprintf(fid,'****************************************************************\n');


 	% save TS_loc to file every 10 minutes
    if sum(times(1:i))-lst > 60*10; % it's been more than 10 mins since last save
        fprintf(fid,'Not finished calculations yet, but saving progress so far to file...\n')
		save('TS_loc.mat','TS_loc','-v7.3');        fprintf(fid,'TS_loc.mat')
		save('TS_loc_ct.mat','TS_loc_ct','-v7.3');  fprintf(fid,', TS_loc_ct.mat')
		save('TS_loc_q.mat','TS_loc_q','-v7.3');    fprintf(fid,', TS_loc_q.mat. All saved.\n')
        lst = sum(times(1:i)); % the last saved time (1)
    end
        
end

%%% Finished calculating!! Aftermath:
fprintf(fid,'!! !! !! !! !! !! !!  Calculation completed at %s !! !! !! !! !! !!\n',datestr(now))
fprintf(fid,'Calculations took a total of %s.\n',benrighttime(sum(times)))

% Aave the local files for subsequent integration into the storage system
fprintf(1,'Saving all results to local .mat files:')
save('TS_loc.mat','TS_loc','-v7.3');        fprintf(fid,' TS_loc.mat')
save('TS_loc_ct.mat','TS_loc_ct','-v7.3');  fprintf(fid, ', TS_loc_ct.mat')
save('TS_loc_q.mat','TS_loc_q','-v7.3');    fprintf(fid,', TS_loc_q.mat. All saved.\n')

if tolog, fclose(fid); end % close the .log file


if agglomerate
    fprintf(1,'Calculation done: Calling TSQ_agglomerate to store back results\n')
    log = 0; % don't bother with a log file
    onlyempty = 0; % only write back to NULL entries in the database
    TSQ_agglomerate(log,onlyempty,dbname)
else
    fprintf(1,'Calculation complete: you can now run TSQ_agglomerate to store results to a mySQL database\n')
end

end