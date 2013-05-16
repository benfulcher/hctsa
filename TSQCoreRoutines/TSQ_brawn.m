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
% (1) tolog: log to file
if nargin<1
	tolog = 0; % by default, do not log to file
end
if tolog
	disp('<TS_brawn>: Logging progess to file');
	fn = ['TS_brawn_' datestr(now,30) '.log'];
	fid = fopen(fn,'w','n');
	disp(['Log file created: ' fn]);
else
	disp('<TS_brawn>: By your command: NOT LOGGING PROGRESS TO FILE?!?!');
end

% (2) Parallel Processing
if nargin<2
	toparallel = 1;
end
if toparallel
	disp('Using parallel processing');
	if tolog
		fprintf(fid, '%s\n', 'Using parallel processing');
	end
else % use single-threaded for loops
	disp('Using single-threaded for loops without any parallel processing on your advice');
	disp('JUST TO REPEAT: WE''RE NOT USING PARALLELIZATION!!!!!');
	if tolog
		fprintf(fid, '%s\n', 'Using single-threaded for loops without any parallel processing');
	end
end

if nargin<3
	dbname = '';
end

if nargin<4
    agglomerate = 1; 
end

%% Update path
% assume path is up to date
% TSQ_updatepath % provides path information for all data, metrics, functions, etc.
% addpath('H:\MATLAB\Time_Series_Resources\TSQ_routines\Core');

%% Read in information from local files
disp('TS_brawn >> Reading in data and guides from file...');
% (As prepared by TS_prepared2.m)

% load the relevent subsegment of the storage: TS_loc
load TS_loc.mat TS_loc
disp('TS_loc loaded');

% load the relevent subsegment of the complete STORE_cts matrix (calculation times)
load TS_loc_ct.mat TS_loc_ct
disp('TS_loc_ct loaded');

% load the quality information of time series: TS_loc_q
load TS_loc_q.mat TS_loc_q
disp('TS_loc_q loaded');

% load the relevant information about the time series
load TS_loc_guides.mat ts_ids_keep tsf tskw tsl m_ids_keep mcode mlab mkw mpoint mlink Mmid Mmlab Mmcode  % tsf tsl tskw tsprep tsmap nts
disp('TS_loc_guides loaded');

%% Definitions
nts = length(ts_ids_keep); % number of time series
nm = length(m_ids_keep); % number of metrics
nMm = length(Mmlab); % number of master functions

%%% THE BRAWN %%%
disp(['Calculation has begun using ' num2str(nts) ' datasets and ' num2str(nm) ' metrics']);
disp(datestr(now));
if tolog
    fprintf(fid,'%s\n',['*** TS_brawn run on ' datestr(now) ' using ' num2str(nts) ' time series and ' num2str(nm) ' metrics']);
end

%% Open parallel processing worker pool
if toparallel
	if matlabpool('size')==0
		matlabpool open;
		disp(['MATLAB parallel processing pool opened with ' num2str(matlabpool('size')) ' and ready to go']);
	else
		disp(['MATLAB parallel processing pool already open. Size: ' num2str(matlabpool('size'))])
	end
end


times = zeros(nts,1); % stores the time taken for each time series to have its metrics calculated (for determining time remaining)
lst1 = 0; % last saved time 1
lst2 = 0; % last saved time 2
for i = 1:nts
	bigtimer = tic;

	qq = TS_loc_q(i,:); % the calculation states of any existing results for the current time series, a line of TS_loc_q
					   	 % NaN indicates a value never before calculated, 1 indicates fatal error before -- try again this time
    tcal = find(isnan(qq) | qq == 1); % indicies of time series awaiting calculation

    if ~isempty(tcal); % some to calculate
        ffi = zeros(length(tcal),1); % output of metrics
		qqi = zeros(length(tcal),1); % label of quality of outputs
		cti = ones(length(tcal),1)*NaN; % calculation times
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
			if size(x,1)==1; % we just need to transpose
				disp(['***** THIS TIME SERIES!!!: ' tsf{i} ' -- IS A 1xN VECTOR. I''D PREFER A Nx1 VECTOR, TO BE HONEST...']);
				disp('I''ll transpose it for you just this once.... but you should try running TS_prep_dimfid.m on your data files');
				if tolog, fprintf(fid,'%s\n',['*** ERROR: a 1xN vector... Transposed.']); end
				x=x';
			else
				disp('******************************************************************************************');
				disp(['MASSIVE ERROR WITH THIS TIME SERIES!!!: ' tsf{i} ' -- IT IS NOT an Nx1 VECTOR. SKIPPING!!']);
				disp('******************************************************************************************');
				if tolog, fprintf(fid,'%s\n',['*** MASSIVE ERROR: A MULTIVARIATE(?) TIME SERIES... SKIPPING...']); end
				continue % the entries in TS_loc and STORE_data will remain NaNs
            end
		end

        %% Pre-Processing
		% AS IT STANDS, PRE-PROCESSING IS NOW ALL JUST Z-SCORING
        % y = zscore(x); % remove the mean and normalize the variance using
                         % the Statistics Toolbox
		y = (x-mean(x))/std(x); % zscoring without the statistics toolbox
		if tolog, fprintf(fid,'%s\n','STANDARD PREPROCESSING: Z-SCORING'); end

		% So we now have the raw time series x and the zscored time series y. The metrics take these as inputs.
        % index sliced variables to minimize the communication overhead in the parallel processing
        partsfi  = tsf{i}; parmcode = mcode(tcal); parmlab = mlab(tcal); parmlink = mlink(tcal);

		disp(['Calculating ' partsfi ' (ts_id = ' num2str(ts_ids_keep(i)) ', N = ' num2str(tsl(i)) ') [' num2str(length(tcal)) ' / ' num2str(nm) ']'])
		if tolog
			fprintf(fid,'%s\n',['Calculating ' num2str(length(tcal)) ' (/' num2str(nm) ') operations']);
		end
		
		
		%% Evaluate master functions in parallel
		% Because of parallelization, we have to evaluate all the master functions *first*
		% Check through the metrics to determine which master functions are relevant for this run
		ip = find(parmlink > 0); % a range of tcal -- that point to a master function
		if ~isempty(ip) % there are some pointers
			
			% masterpool = cell(length(ip),2);
			% for j=1:length(ip) % loops over all pointer metrics
			% 	tmpi = strfind(parmcode{ip(j)},'.');
			% 	masterpool{j,1} = parmcode{ip(j)}(1:tmpi-1); % removes the structure field label; the master label
			% 	masterpool{j,2} = char(mMf(strmatch(masterpool{j,1},mMl,'exact'))); % the code to call
			%             end
			%             
			%             [meismasterl a b] = unique(masterpool(:,1)); % these are the master labels to point to later
			%             meismasterf = masterpool(a,2); % these are the master functions that need to be evaluated
            
			%             if any(arrayfun(@(x)strcmp(x,''),masterpool(:,2)))
			%                 char(masterpool(arrayfun(@(x)strcmp(x,''),masterpool(:,2)),1))
			%                 disp('INP_mets master files have not been formatted correctly. Exiting');
			% keyboard; return
			%             end
			%             if length(a)~=length(unique(masterpool(:,2)))
			%                 % this is an error -- there is something wrong with how the
			%                 % INP_mets file has been prepared
			%                 [meismasterl a1 b1] = unique(masterpool(:,1)); % these are the master labels to point to later
			%                 [meismasterf a2 b2] = unique(masterpool(:,2)); % these are the master functions that need to be evaluated
			%                 char(masterpool(setxor(a1,a2),:)) % there are the problem ones; ideally a1 and and a2 are identical
			%                 disp('INP_mets master files have not been formatted correctly. Exiting');
			% keyboard; return
			%             end
            
			

			% For each implicated master function -- put the structure in this element
			
			% Ouput structures stored in Moutput -- unused master functions left with empty elements that should not be referred to
			Moutput = cell(nMm,1);
			% Calculation times for the master structures stored in Mcts
			Mcts = zeros(nMm,1);
			
			% Mitocalc are the indicies of masters (Mmlab/Mmcode) that need to be called (i.e., they'll be called later)
			Mitocalc = unique(parmlink(parmlink>0));
			nMtocalc = length(Mitocalc); % number of implicated master functions

			disp(['Master metrics are being evaluated first: all ' num2str(nMtocalc) ' of them...']);
			if tolog, fprintf(fid,'%s\n',['Evaluating ' num2str(nMtocalc) ' master metrics']); end % to log
			
		    % store in temporary variables then map back
            MoutputPAR = cell(nMtocalc,1);
            MctsPAR = zeros(nMtocalc,1);
			
			
			% Evaluate the master functions
			if toparallel
	            parfor j = 1:nMtocalc
					jj = Mitocalc(j); % the index of the master array to evaluate
					try
						mastertimer = tic;
						% disp(Mmcode{jj}) % FOR DEBUGGING -- SHOW CODE BEFORE EVALUATING TO DETERMINE PROBLEMS
						MoutputPAR{j} = pareval(x,y,Mmcode{jj});
						% for not-applicable/'real NaN', outputs a NaN, otherwise should output a
						% structure with components to be called below by pointer metrics.
						MctsPAR(j) = toc(mastertimer);
                    catch emsg
						disp(['  ** ** ERROR EVALUATING MASTER FILE: ' Mmlab{jj}]);
                        disp(emsg.message)
						% masterdat{j} didn't evaluate properly -- remains an empty cell entry.
					end
                end
			else % just a single-thread for loop
	            for j = 1:nMtocalc
					jj = Mitocalc(j); % the index of the master array to evaluate
					disp(Mmcode{jj}); % for error checking
					try
						mastertimer = tic;
						MoutputPAR{j} = pareval(x,y,Mmcode{jj});
						% for not-applicable/'real NaN', outputs a NaN, otherwise should output a
						% structure with components to be called below by pointer metrics.
						MctsPAR(j) = toc(mastertimer);
                    catch emsg
						disp(['  ** ** ERROR EVALUATING MASTER FILE: ' Mmlab{jj}]);
						disp(emsg.message)
                        % masterdat{j} didn't evaluate properly -- remains an empty cell entry.
					end
                end
			end
			
			% map back from PAR versions to real versions
			% keyboard
            Moutput(Mitocalc) = MoutputPAR;
            Mcts(Mitocalc) = MctsPAR;
			
			
			disp('Master metrics evaluated.');
		else
			% no master metrics need to be calculed.
			Moutput={}; Mcts={}; % This initiaition is necessary for the next parfor loop
		end


		%% GO GO GO!!!
		if toparallel
	        parfor j = 1:length(tcal)
				if parmlink(j)>0 % pointer to a master function
					% tmpi=strfind(parmcode{j},'.');
					% theml=parmcode{j}(1:tmpi-1); % this is the master label
					% thest=parmcode{j}(tmpi+1:end); % this is the structure element to retrieve
					% linky=strmatch(theml,meismasterl,'exact'); % find the master structure to link to
					try
						% retrive from master structure:
						if ~isstruct(Moutput{parmlink(j)}) && isnan(Moutput{parmlink(j)});
							% All master function outputs are to be set to real NaNs (unsuitable)
							ffi(j) = NaN;
		                    cti(j) = Mcts(parmlink(j));
						else % (normal -- retrieve required element from master structure)
							thedot = strfind(parmcode{j},'.');
							thest = parmcode{j}(thedot+1:end);
		                    ffi(j) = parevalM(Moutput{parmlink(j)},['themasterdat.' thest]);
							qqi(j) = 0; % good output
		                    cti(j) = Mcts(parmlink(j));
							% if isinf(ffi(j)), ffi(j) = NaN; end % real Inf stored as if a real NaN...
						end
                    catch emsg
						disp(['Error evaluating link to Master structure ' Mmlab{parmlink(j)} ' by ' parmcode{j}]);
% 						disp(emsg.message)
                        ffi(j) = 0;
						qqi(j) = 1; % fatal error code
					end
				else % A regular metric
% 					disp(parmcode{j}); % for error checking
		            try
						metrictimer = tic;
						ffi(j) = pareval(x,y,parmcode{j});
						qqi(j) = 0;
						cti(j) = toc(metrictimer);
						% if isinf(ffi(j)), ffi(j)=NaN; end % real Inf stored as if real NaN...
		            catch
		                disp(['Fatal error ' partsfi ' || ' parmcode{j} '; stored as NaN']);
						ffi(j) = 0;
						qqi(j) = 1; % fatal error code
		            end
				end
	        end
		else
	        for j = 1:length(tcal)
				if parmlink(j)>0 % pointer to a master function
					% tmpi=strfind(parmcode{j},'.');
					% theml=parmcode{j}(1:tmpi-1); % this is the master label
					% thest=parmcode{j}(tmpi+1:end); % this is the structure element to retrieve
					% linky=strmatch(theml,meismasterl,'exact'); % find the master structure to link to
					try
						% retrive from master structure:
						if ~isstruct(Moutput{parmlink(j)}) && isnan(Moutput{parmlink(j)});
							% All master function outputs are to be set to real NaNs (unsuitable)
							ffi(j) = NaN;
		                    cti(j) = Mcts(parmlink(j));
						else % (normal -- retrieve required element from master structure)
							thedot = strfind(parmcode{j},'.');
							thest = parmcode{j}(thedot+1:end);
		                    ffi(j) = parevalM(Moutput{parmlink(j)},['themasterdat.' thest]);
							qqi(j) = 0; % good output
		                    cti(j) = Mcts(parmlink(j));
							% if isinf(ffi(j)), ffi(j) = NaN; end % real Inf stored as if a real NaN...
						end
	                catch
						disp(['Error evaluating link to Master structure ' Mmlab{parmlink(j)} ' by ' parmcode{j}]);
						ffi(j) = 0;
						qqi(j) = 1; % fatal error code
					end
				else % A regular metric
					disp(parmcode{j}); % for error checking
		            try
						metrictimer = tic;
						ffi(j) = pareval(x,y,parmcode{j});
						qqi(j) = 0;
						cti(j) = toc(metrictimer);
						% if isinf(ffi(j)), ffi(j)=NaN; end % real Inf stored as if real NaN...
		            catch
		                disp(['Fatal error ' partsfi ' || ' parmcode{j}]);
						ffi(j) = 0;
						qqi(j) = 1; % fatal error code
		            end
				end
	        end
		end
        
		%% Treat error coding -- qqi
		% (*) Errorless calculation: q=0, output=real number (default)
		% (*) Fatal error: q=1, output = 0; (this is done already in the code above)
		% (*) Real NaN: q=2, output = 0;
		RR = isnan(ffi); % NaN
		if ~isempty(RR)
			qqi(RR) = 2;
			ffi(RR) = 0;
		end
		% (*) Real Inf: q=3, output = 0;
		RR = (isinf(ffi) & ffi>0); % Inf
		if ~isempty(RR)
			qqi(RR) = 3;
			ffi(RR) = 0;
		end
		% (*) Real -Inf: q=4, output=0;
		RR = (isinf(ffi) & ffi<0);
		if ~isempty(RR)
			qqi(RR) = 4;
			ffi(RR) = 0;
		end
		% (*) Complex Number: q=5
		RR = imag(ffi)~=0;
		if ~isempty(RR)
			qqi(RR) = 5;
			ffi(RR) = 0;
		end

		
		% Store this time series information back to global structures
		
        TS_loc(i,tcal) = ffi; % store new ones in TS_loc
		TS_loc_ct(i,tcal) = cti; % store calculation times in TS_loc_ct
		TS_loc_q(i,tcal) = qqi; % store quality labels in TS_loc_q
		% Note that the calculation time for pointers is the total calculation time for the full master structure.		


		% Set time for this data set, for run-timing purposes
		times(i) = toc(bigtimer);
		
		% LOG
		ngood = sum(qqi==0); % the number of calculated metrics that returned real outputs without errors
		nerror = sum(qqi==1);
		nother = sum(qqi>1);

	    disp('****************************************************************')
	    disp(['; ; ; : : : : ; ; ; ;    ' datestr(now) '     ; ; ; ; : : : ; ; ;'])
	    disp(['; ; ; : : : : ; ; ; ; ; ;   ts_id: ' num2str(ts_ids_keep(i)) '     ; ; ; ; : ; ; ; ; ; ; ;'])
	    disp(['oOoOoOoOoOo  ' partsfi ' (' num2str(tsl(i)) ')  oOoOoOoOoOoOo'])
	    disp([num2str(ngood) ' good outputs... ' num2str(nerror) ' errors;  ' ...
	     		num2str(nother) ' other ouputs stored. [ ' num2str(length(tcal)) ' / ' num2str(nm) ' ]']);
		disp(['Calculation for this data set took ' benrighttime(times(i))])
	    disp(['    - - - - - - - - - - -  ' num2str(i) '  /  ' num2str(nts) '   - - - - - - - - - - -'])
		disp(['    - - - - - - - -  ' benrighttime(((nts-i)*mean(times(1:i)))) ' remaining - - - - - - - - -'])
	    fprintf('%s\n\n\n','****************************************************************')
		
		% File log
		if tolog
			% Summary
			fprintf(fid,'%s\n','****************************************************************');
		    fprintf(fid,'%s\n',['; ; ; : : : : ; ; ; ;    ' datestr(now) '     ; ; ; ; : : : ; ; ;']);
		    fprintf(fid,'%s\n',['oOoOoOoOoOo  ' partsfi ' (' num2str(tsl(i)) ')  oOoOoOoOoOoOo']);
		    fprintf(fid,'%s\n',[num2str(ngood) ' real outputs...' num2str(nerror) ' errors;  ' ...
		     					num2str(nother) ' other outputs stored. [ ' num2str(length(tcal)) ' / ' num2str(nm) ' ]']);
			fprintf(fid,'%s\n',['Calculation for this data set took ' benrighttime(times(i))]);
		    fprintf(fid,'%s\n',['    - - - - - - - -  ' num2str(i) '  /  ' num2str(nts) '   - - - - - - - -']);
			fprintf(fid,'%s\n',['    - - - - - - - -  ' benrighttime(((nts-i)*mean(times(1:i)))) ' remaining - - - - - - - - -']);
		    fprintf(fid,'%s\n\n','****************************************************************');
			
			% NaNs
			% rz = find(isnan(ffi));
			% if ~isempty(rz)
			% 	fprintf(fid,'%s\n',['(**)' num2str(length(rz)) ' NaNs / ' num2str(length(tcal)) ' calculated:']);
			% 	for k=1:length(rz)
			% 		fprintf(fid,'%s\n',parmlab{rz(k)});
			% 	end
			% else
			% 	fprintf(fid,'%s\n',['(**) no NaNs ( 0 NaNs /' num2str(length(tcal)) ' calculated)']);
			% end
			% % Infs
			% rz = find(isinf(ffi));
			% if ~isempty(rz)
			% 	fprintf(fid,'%s\n',['(**)' num2str(length(rz)) ' Infs / ' num2str(length(tcal)) ' calculated:']);
			% 	for k=1:length(rz)
			% 		fprintf(fid,'%s\n',parmlab{rz(k)});
			% 	end
			% else
			% 	fprintf(fid,'%s\n',['(**) no Infs calculated ( 0 /' num2str(length(tcal)) ')']);
			% end
			fprintf(fid,'%s\n\n\n\n','****************************************************************');
		end
		
    else % nothing to calculate (empty tcal)
		% Timer/updater;
	    times(i) = toc(bigtimer);
	

		disp('****************************************************************')
	    disp(['; ; ; : : : : ; ; ; ;    ' datestr(now) '     ; ; ; ; : : : ; ; ;'])
	    disp(['oOoOoOoOoOo  ' tsf{i} ' (' num2str(tsl(i)) ')  oOoOoOoOoOoOo'])
		disp(['0O0O0O0O0O0 nothing calculated! All ' num2str(nm) ' already complete!!  0O0O0O0O0O0'])
	    disp(['    - - - - - - - -  ' num2str(i) '  /  ' num2str(nts) '   - - - - - - - -'])
		disp(['    - - - - - - - -  ' benrighttime(((nts-i)*mean(times(1:i)))) ' remaining - - - - - - - - -'])
	    disp('****************************************************************')
	
		if tolog
			fprintf(fid,'%s\n','****************************************************************');
		    fprintf(fid,'%s\n',['; ; ; : : : : ; ; ; ;    ' datestr(now) '     ; ; ; ; : : : ; ; ;']);
		    fprintf(fid,'%s\n',['oOoOoOoOoOo  ' tsf{i} ' (' num2str(tsl(i)) ')  oOoOoOoOoOoOo']);
			fprintf(fid,'%s\n',['0O0O0O0O0O0 nothing calculated! All ' num2str(nm) ' already complete!!  0O0O0O0O0O0']);
		    fprintf(fid,'%s\n',['    - - - - - - - -  ' num2str(i) '  /  ' num2str(nts) '   - - - - - - - -']);
			fprintf(fid,'%s\n',['    - - - - - - - -  ' benrighttime(((nts-i)*mean(times(1:i)))) ' remaining - - - - - - - - -']);
		    fprintf(fid,'%s\n\n\n\n','****************************************************************');
		end 
    end



 	% save TS_loc to file every 10 minutes
    if sum(times(1:i))-lst1>60*10;
		save('TS_loc.mat','TS_loc','-v7.3')
		save('TS_loc_ct.mat','TS_loc_ct','-v7.3')
		save('TS_loc_q.mat','TS_loc_q','-v7.3')
        lst1=sum(times(1:i)); % the last saved time (1)
    end
    
    if sum(times(1:i))-lst2>6*3600; % every six hours save backups -- in case files become corrupt -- can retrieve
		disp('Saving Backup...');
		timenow = datestr(now,30);

        save(['TS_loc_' timenow '.mat'],'TS_loc','-v7.3');
		save(['TS_loc_ct_' timenow '.mat'],'TS_loc_ct','-v7.3');
		save(['TS_loc_q_' timenow '.mat'],'TS_loc_q','-v7.3');
		disp(['******* saved TS_loc_' timenow '.mat (' benrighttime(lst2) ') **********'])
        disp(['******* saved TS_loc_ct_' timenow '.mat (' benrighttime(lst2) ') **********'])
        disp(['******* saved TS_loc_q_' timenow '.mat (' benrighttime(lst2) ') **********'])

		if tolog
			fprintf(fid,'%s\n',['******* saved TS_loc_' timenow '.mat (' benrighttime(lst2) ') **********']);
	        fprintf(fid,'%s\n',['******* saved TS_loc_ct_' timenow '.mat (' benrighttime(lst2) ') **********']);
	        fprintf(fid,'%s\n',['******* saved TS_loc_q_' timenow '.mat (' benrighttime(lst2) ') **********']);
		end

        lst2 = sum(times(1:i)); % the last saved time (2)
    end
    
end
% if toparallel
% 	matlabpool close
% end

%%% Aftermath
disp(['! ! ! !! !! !! !! ! ! ! !  DONE AT  ' datestr(now) '     ! ! ! ! !! !! !! ! ! !'])
disp(['The calculation took ' benrighttime(sum(times))])

if tolog
	fprintf(fid,'%s\n',['Calculations took ' benrighttime(sum(times)) ' in total']);
	fprintf(fid,'%s\n',['You''ve been a great audience. Goodnight.']);
	fclose(fid);
end

% save the local files for subsequent integration into the storage system
save('TS_loc.mat','TS_loc','-v7.3');        disp('saved TS_loc')
save('TS_loc_ct.mat','TS_loc_ct','-v7.3');  disp('saved TS_loc_ct')
save('TS_loc_q.mat','TS_loc_q','-v7.3');    disp('saved TS_loc_q')

if agglomerate
    disp('Calculation done: Calling TSQ_agglomerate to store back results')
    % Save outputs back to database dbc
    TSQ_agglomerate(0,dbname) % 0 -- don't both with a log file
else
    disp(['Calculation done: you''ve told me not to but you should call '...
            'TSQ_agglomerate to store back results'])
end

% Update percentage calculated statistics in TimeSeries, TimeSeriesKeywords, Operations, and OperationKeywords tables
% reply = input(['Update percentage calculated, mean calculation time statistics? [''y'' for yes, ''t'' for just time series]'], 's');
% if strcmp(reply,'y')
% 	SQL_fillfromresults(ts_ids_keep,m_ids_keep,[1,1,1],[1,1,1,1],dbname) % Updates percentagecalc, percentagegood, meancalctime
% elseif strcmp(reply,'t')
% 	SQL_fillfromresults(ts_ids_keep,[],[1,1,1],[1,1,0,0],dbname) % Updates percentagecalc, percentagegood, meancalctime
% end

end