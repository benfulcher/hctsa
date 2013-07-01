function TSQ_brawn(tolog,toparallel,bevocal)
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

% Use Matlab's Parallel Computing toolbox?
if nargin < 2
	toparallel = 0;
end
if toparallel
	fprintf(fid,'Computation will be performed across multiple cores using Matlab''s Parallel Computing Toolbox\n')
else % use single-threaded for loops
	fprintf(fid,'Computations will be performed serially without parallelization\n')
end

% be vocal?
if nargin < 3
    bevocal = 1; % write back lots of information to screen
    % prints every piece of code evaluated (nice for error checking)
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
    % first check that the user can use the Parallel Computing Toolbox:
    heyo = which('matlabpool');
    if isempty(heyo)
        fprintf(1,'Parallel Computing Toolbox not found -- cannot perform computations across multiple cores\n')
        toparallel = 0;
    else
        if matlabpool('size')==0
        	matlabpool open;
        	fprintf(fid,'MATLAB parallel processing pool opened with %u and ready to go',matlabpool('size'))
        else
        	fprintf(fid,'MATLAB parallel processing pool already open. Size: %u\n',matlabpool('size'))
        end
    end
end

times = zeros(nts,1); % stores the time taken for each time series to have its metrics calculated (for determining time remaining)
lst = 0; % last saved time

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
        
		%% Basic checking on x
		% univariate and Nx1
		if size(x,2) ~= 1
			if size(x,1) == 1
				fprintf(fid,'***** The time series %s is a row vector. Not sure how it snuck through the gates, but I need a column vector...\n',tsf{i});
				fprintf(fid,'I''ll transpose it for you for now....\n');
				x = x';
			else
				fprintf(fid,'******************************************************************************************\n')
				fprintf(fid,'MASSIVE ERROR WITH THIS TIME SERIES!!!: %s -- is it multivariate or something weird???. Skipping it!!\n', tsf{i});
				fprintf(fid,'******************************************************************************************\n');
				continue % skip to the next time series; the entries for this time series in TS_loc etc. will remain NaNs
            end
		end

        %% Pre-Processing
		% y is a z-scored transformation of the time series
		y = (x-mean(x))/std(x); % z-scoring without using a Statistics Toolbox license (i.e., the zscore function)

		% So we now have the raw time series x and the z-scored time series y. Operations take these as inputs.
        % index sliced variables to minimize the communication overhead in the parallel processing
        partsfi  = tsf{i}; parmcode = mcode(tcal); parmlab = mlab(tcal); parmlink = mlink(tcal);

		% Display information
		whichtsf = which(tsf{i});
		fprintf(fid,'\n\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n')
	    fprintf(fid,'; ; ; : : : : ; ; ; ;    %s     ; ; ; ; : : : ; ; ;\n',datestr(now))
	    fprintf(fid,'- - - - - - - -  Loaded time series %u / %u - - - - - - - -\n',i,nts)
		fprintf(fid,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n')
		fprintf(fid,'%s\n',whichtsf)
		fprintf(fid,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n')
		fprintf(fid,'Preparing to calculate %s (ts_id = %u, N = %u samples) [computing %u / %u operations]\n',partsfi,ts_ids_keep(i),tsl(i),ncal,nm)
	    fprintf(fid,'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')

		%% Evaluate master functions in parallel
		% Because of parallelization, we have to evaluate all the master functions *first*
		% Check through the metrics to determine which master functions are relevant for this run

		% For each implicated master operation -- put the structure in this element
		Moutput = cell(nMm,1); % Ouput structures
		Mcts = zeros(nMm,1); % Calculation times for the master structures
		
		Mitocalc = unique(parmlink(parmlink > 0)); % indicies of masters (Mmlab/Mmcode) that need to be called
		nMtocalc = length(Mitocalc); % number of master operations to calculate

		fprintf(fid,'Evaluating %u master operations...\n',nMtocalc);
		
	    % Store in temporary variables for parfor loop then map back later
        Moutput_tmp = cell(nMtocalc,1);
        Mcts_tmp = zeros(nMtocalc,1);
		
		% Evaluate all the master operations
		if toparallel
            parfor j = 1:nMtocalc
                [Moutput_tmp{j}, Mcts_tmp(j)] = TSQ_brawn_masterloop(x,y,Mmcode{Mitocalc(j)},fid,bevocal);
            end
        else
            for j = 1:nMtocalc
                [Moutput_tmp{j}, Mcts_tmp(j)] = TSQ_brawn_masterloop(x,y,Mmcode{Mitocalc(j)},fid,bevocal);
            end
		end
		
		% map back from temporary versions to the real versions
        Moutput(Mitocalc) = Moutput_tmp;
        Mcts(Mitocalc) = Mcts_tmp;
		
		fprintf(fid,'%u master operations evaluated ///\n\n',nMtocalc);

		%% Now assign all the results to the right operations
		if toparallel
	        parfor j = 1:ncal
                [ffi(j), qqi(j), cti(j)] = TSQ_brawn_oploop(x, y, parmlink(j), Moutput,...
                                                    Mcts,Mmlab,parmcode{j},fid);
            end
            % TSQ_brawn_oploop(x, y, moplink, Moutput, Mcts, Mmlab, parmcodej, partsfi, fid, bevocal)
		else
            for j = 1:ncal
                try
                    [ffi(j), qqi(j), cti(j)] = TSQ_brawn_oploop(x, y, parmlink(j), Moutput,...
                                                Mcts,Mmlab,parmcode{j},fid);
                catch
                    fprintf(1,'Error retrieving element %s from %s\n',parmcode{j},Mmlab{parmlink(j)})
                    keyboard
                end
                % [ffi(j),qqi(j),cti(j)] = TSQ_brawn_oploop(x, y, parmlink(j), Moutput{parmlink(j)},...
                %                                 Mcts(parmlink(j)),Mmlab{parmlink(j)},fid,bevocal);
                % When all operations have masters, we can make this nicer, more like ^
            end
        end
        %             for j = 1:ncal
        % if parmlink(j) > 0 % pointer to a master function
        %     try
        %         % retrive from master structure:
        %         if ~isstruct(Moutput{parmlink(j)}) && isnan(Moutput{parmlink(j)});
        %             % All master function outputs are to be set to real NaNs (unsuitable)
        %             ffi(j) = NaN;
        %                             cti(j) = Mcts(parmlink(j));
        %         else % (normal -- retrieve required element from master structure)
        %                             [~,thest] = strtok(parmcode{j},'.'); thest = thest(2:end); % remove the '.'
        %                             ffi(j) = parevalM(Moutput{parmlink(j)},['themasterdat.' thest]);
        %             qqi(j) = 0; % good, real-valued output
        %                             cti(j) = Mcts(parmlink(j));
        %         end
        %                     catch emsg
        %         fprintf(fid,'Error evaluating link to Master structure %s by %s\n',Mmlab{parmlink(j)},parmcode{j});
        %                         fprintf(fid,'%s\n',emsg.message)
        %         ffi(j) = 0; qqi(j) = 1; % Fatal error code: 1
        %     end
        % else % A regular, single-output operation
        %                     if bevocal, fprintf(fid,'%s\n',parmcode{j}); end % for error checking
        %                     try
        %         operationtimer = tic;
        %         ffi(j) = pareval(x,y,parmcode{j});
        %         cti(j) = toc(operationtimer);
        %         qqi(j) = 0;
        %                     catch
        %                         fprintf(fid,'Fatal error %s || %s\n',partsfi,parmcode{j});
        %         ffi(j) = 0; qqi(j) = 1; % fatal error code = 1
        %                     end
        % end
        %             end
        
        
		%% Code special values:

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
		ngood = sum(qqi == 0); % the number of calculated metrics that returned real outputs without errors
		nerror = sum(qqi == 1); % number of fatal errors encountered
		nother = sum(qqi > 1); % number of other special outputs
    end
    
    times(i) = toc(bigtimer); % the time taken to calculate (or not, if ncal = 0) operations for this time series


    % Calculation complete: print information about this time series calculation
	fprintf(fid,'********************************************************************\n')
    fprintf(fid,'; ; ; : : : : ; ; ; ;    %s     ; ; ; ; : : : ; ; ;\n',datestr(now))
    fprintf(fid,'oOoOoOo Calculation complete for %s (N = %u samples)  oOoOoOoOo\n',tsf{i},tsl(i));
    if ncal > 0 % Some amount of calculation was performed
	    fprintf(fid,'%u real-valued outputs, %u errors, %u other outputs stored. [%u / %u]\n',...
	     					ngood,nerror,nother,ncal,nm);
		fprintf(fid,'Calculations for this time series took %s.\n',benrighttime(times(i),1));
    else
    	fprintf(fid,'Nothing calculated! All %u operations already complete!!  0O0O0O0O0O0\n',nm);
    end
    if i < nts
        fprintf(fid,'- - - - - - - -  %u time series remaining - - - - - - - -\n',nts-i);
    	fprintf(fid,'- - - - - - - -  %s remaining - - - - - - - - -\n',benrighttime(((nts-i)*mean(times(1:i))),1));
    else % the final time series
        fprintf(fid,'- - - - - - - -  All time series calculated! - - - - - - - -\n',nts-i);
    end
    fprintf(fid,'********************************************************************\n');

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
fprintf(fid,'Calculations took a total of %s.\n',benrighttime(sum(times),1))

% Aave the local files for subsequent integration into the storage system
fprintf(1,'Saving all results to local .mat files:')
save('TS_loc.mat','TS_loc','-v7.3');        fprintf(fid,' TS_loc.mat')
save('TS_loc_ct.mat','TS_loc_ct','-v7.3');  fprintf(fid, ', TS_loc_ct.mat')
save('TS_loc_q.mat','TS_loc_q','-v7.3');    fprintf(fid,', TS_loc_q.mat. All saved.\n')

if tolog, fclose(fid); end % close the .log file

% if agglomerate
%     fprintf(1,'\n\nCalculation done: Calling TSQ_agglomerate to store back results\n')
%     writewhat = 'null'; % only write back to NULL entries in the database
%     TSQ_agglomerate(writewhat,dbname,tolog)
% else
fprintf(1,'Calculation complete: you can now run TSQ_agglomerate to store results to a mySQL database\n')
% end

end