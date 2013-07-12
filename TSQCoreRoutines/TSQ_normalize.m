function TSQ_normalize(normopt,trimopt,subs,trainset)
% Updated version of TS_normalize: Ben Fulcher August 2009
% Reads in from TS_loc and guides; writes to TS_loc_N and guides
% These will be different to TS_loc because this is *already trimmed*
% as well as normalized. TS_cluster now has the same size as TS_loc_N, but
% just reordered.
% There are options for how to normalize (normopt);
% and also how to trim bad entries (trimopt).
% [opt] subs: only normalize and trim a subset of TS_loc. This can be used,
% for example, to analyze just a subset of the full space, which can
% subsequently be clustered and further subsetted using TS_cluster2...
% For example, can choose a subset using SUB_autolabel2 to get only sound
% time series. Then can try and distinguish 'animalsounds' from amongst the
% other sound data using TS_prune4
% subs in the form {[rowrange],[columnrange]} (to keep, from TS_loc).
% Ben Fulcher August 2009
% Ben Fulcher 9/12/2009 -- Updated for new mySQL setup
% Ben Fulcher 12/9/2011 -- added trainset input -- indices of training
%                           indicies (will only derive normalization
%                           transformation on these indices)

% Reads in local TS_loc data, and saves the `normalized' version as TS_loc_N
% The normalization is all about a rescaling to the [0,1] interval for visualization and clustering
% Normalization rules are imported from the metric file, as saved in TS_guide_met.mat

% If set one of the trimopts{1} to 1, will have no bad values in your
% matrix.

% Ben Fulcher 26/1/2011 -- Added incorporation of TS_loc_q to determination
% of bad entries, rather than simply values missing in the matrix.
% Ben Fulcher 28/1/2011 -- Removed filtering by high covariance -- can do
% this elsewhere (e.g., TSQ_preprocess), added TS_loc_q_N

%% Argument Checking, Preliminaries
if nargin < 1 || isempty(normopt)
    fprintf(1,'Using the default, scaled quantile-based sigmoidal transform -- ''scaledSQzscore''\n')
%     normopt = 'scaledsigmoid';
    normopt = 'scaledSQzscore';
end

if nargin < 2 || isempty(trimopt)
    trimopt = [0.90, 1]; % (default): remove less than 90%-good time series, & then less than 
                        % 100%-good metrics.
end
if iscell(trimopt) % still using the cell input of previous TSQ_normalize
    trimopt = trimopt{1};
end
fprintf(1,['Removing time series with more than %.2f%% special-valued outputs, ' ...
            'and operations with more than %.2f%%-special-valued outputs\n'],trimopt(1)*100,trimopt(2)*100);

if nargin < 3
    subs = {}; % don't subset
end

if nargin < 4
    trainset = []; % normalize on the full set
end

%% Read in information from local files
% (As prepared by TS_prepare.m)
fprintf(1,'Reading in local files...');
load TS_loc.mat TS_loc % the local data matrix retrieved from the database
F = TS_loc; clear TS_loc; fprintf(1,' TS_loc');
load TS_loc_q.mat TS_loc_q; fprintf(1,', TS_loc_q'); % quality labels for F
load TS_loc_guides.mat; fprintf(1,', TS_loc_guides.\n');
% m_ids_keep, ts_ids_keep, tsf, tskw, tsl, mcode, mlab, mkw, mpoint, mlink, Mmid, Mmlab, Mmcode

%% (0) SUBSET USING GIVEN INDICIES
if ~isempty(subs)
    kr0 = subs{1}; % rows to keep (0)
    if isempty(kr0),
        kr0 = 1:size(F,1);
    else
        fprintf(1,'Filtered down time series by given subset; from %u to %u.\n',...
                    size(F,1),length(kr0))
        F = F(kr0,:);
        TS_loc_q = TS_loc_q(kr0,:);
    end
    
    kc0 = subs{2}; % columns to keep (0)
    if isempty(kc0),
        kc0 = 1:size(F,2);
    else
        fprintf(1,'Filtered down operations by given subset; from %u to %u.\n',...
            size(F,2),length(kc0))
        F = F(:,kc0);
        TS_loc_q = TS_loc_q(:,kc0);
    end
else
    kr0 = 1:size(F,1);
    kc0 = 1:size(F,2);
end

%% (1) TRIM THE BASTARD

% (i) NaNs in TS_loc mean values uncalculated in the matrix.
F(~isfinite(F)) = NaN; % convert all nonfinite values to NaNs for consistency
% need to incorporate knowledge of bad entries in TS_loc_q and filter these
% out:
notgood = (TS_loc_q > 0); % QualityCode > 0 means some special value (NaN, Inf, error, ...)
F(TS_loc_q > 0) = NaN;
fprintf(1,'There were %u special values in TS_loc_q\n',sum(notgood(:)))
% now all bad values are NaNs, we can get on with filtering them out


% (i) Filter based on proportion of bad entries. If either threshold is 1,
% the resulting matrix is guaranteed to be free from bad values entirely.
[badr, badc] = find(isnan(F));
thresh_r = trimopt(1); thresh_c = trimopt(2);
if thresh_r > 0 % if 1, then even the worst are included
    [badr, ~, rj] = unique(badr); % neat code, but really slow to do this
%     unique... Loop instead
    % (ii) remove rows with more than a proportion thresh_r bad values
    badrp = zeros(length(badr),1); % stores the number of bad entries
    for i = 1:length(badr)
        badrp(i) = length(find(rj==i));
    end
    badrp = badrp/size(F,2);
    xkr1 = badr(badrp >= 1 - thresh_r); % don't keep rows (1) if fewer good values than thresh_r
    kr1 = setxor(1:size(F,1),xkr1);

    if ~isempty(kr1)
        if ~isempty(xkr1)
            fprintf(1,['Removed time series with fewer than %4.2f%% good values:'...
                            ' from %u to %u\n'],thresh_r*100,size(F,1),length(kr1))
            % display filtered times series to screen:
            fprintf(1,'Lost %u time series: %s\n',size(F,1)-length(kr1),bencat(tsf(xkr1),','))
        else
            fprintf(1,'All %u time series had greater than %4.2f%% good values. Keeping them all.\n', ...
                            size(F,1),thresh_r*100)
        end
        F = F(kr1,:); % ********************* KR1 ***********************
        TS_loc_q = TS_loc_q(kr1,:);
    else
        error('No time series had more than %4.2f%% good values.',thresh_r*100)
    end
else
    fprintf(1,'No filtering of time series based on proportion of bad values.\n')
    kr1 = (1:size(F,1));
end

if thresh_c > 0
    if thresh_r > 0 && ~isempty(kr1) % did row filtering and removed some
        [~, badc] = find(isnan(F)); % have to recalculate indicies
    end
    [badc, ~, cj] = unique(badc);
    % (iii) Remove metrics that are more than thresh_c bad
    badcp = zeros(length(badc),1); % stores the number of bad entries
    for i = 1:length(badc), badcp(i) = length(find(cj==i)); end
    badcp = badcp/size(F,1);
    xkc1 = badc(badcp >= 1-thresh_c); % don't keep columns if fewer good values than thresh_c
    kc1 = setxor(1:size(F,2),xkc1); % keep columns (1)
    
    if ~isempty(kc1)
        if ~isempty(xkc1)
            fprintf(1,'Removed operations with fewer than %5.2f%% good values: from %u to %u\n',thresh_c*100,size(F,2),length(kc1))
            fprintf(1,'Lost %u operations: %s\n',size(F,2)-length(kc1),bencat(mlab(xkc1),','))
        else
            fprintf(1,['All operations had greater than %5.2f%% good values; ' ...
                    'keeping them all :-)'],thresh_c*100)
        end

        F = F(:,kc1); % *********************** kc1 *********************
        TS_loc_q = TS_loc_q(:,kc1);
        
    else
        error('No operations had fewer than %u%% good values.',thresh_c*100)
    end
else
    % fprintf(1,'No filtering of operations based on proportion of bad values\n')
    kc1 = (1:size(F,2));
end

% (ii) Remove operations that are constant across the time series dataset
if size(F,1) > 1 % otherwise just a single time series remains and all will be constant!
    crap_op = zeros(size(F,2),1);
    for j = 1:size(F,2)
        crap_op(j) = (range(F(~isnan(F(:,j)),j)) < eps);
    end
    kc2 = find(crap_op==0); % kept column (2)

    if ~isempty(kc2)
        if length(kc2) < size(F,2)
            fprintf(1,'Filtered down operations with near-constant outputs: from %u to %u',...
                             size(F,2),length(kc2))
            F = F(:,kc2); % ********************* KC2 **********************
            TS_loc_q = TS_loc_q(:,kc2);
        end
    else
        error('All %u operations produced constant outputs on the %u time series?!',length(kc2),size(F,1))
    end
else
    % just one time series remains: keep all operations
    kc2 = ones(1,size(F,2));
end

% (ii) Remove time series with constant feature vectors
crap_ts = zeros(size(F,1),1);
for j = 1:size(F,1)
    crap_ts(j) = (range(F(j,~isnan(F(j,:)))) < eps);
end
kr2 = find(crap_ts == 0); % kept column (2)

if ~isempty(kr2)
    if (length(kr2) < size(F,1))
        fprintf(1,'Filtered out time series with constant feature vectors: from %u to %u\n',...
                            size(F,1),length(kr2))
        F = F(kr2,:); % ********************* KR2 **********************
        TS_loc_q = TS_loc_q(kr2,:);
    end
else
    error('All time series have constant feature vectors?!')
end

% % (iii) Trim based on high covariances
% cov_thresh_r = trimopt{2}(1);
% cov_thresh_c = trimopt{2}(2);
% if cov_thresh_r>0
%     F_c=corrcoef(F'); % I *think* this transpose is appropriate
%     [xi xj]=find(abs(F_c)>cov_thresh_r);
%     x = [xi xj];
%     r=setxor(1:size(x,1),find(x(:,1)>=x(:,2))); % do not include diagonal or below-diagonal entries
%     if ~isempty(r) % some rows are highly correlated
%         x=x(r,:);
%         kr3=setxor(1:size(F,1),unique(x(:,1)));
%         % have a mx2 list of indicies pairing rows with common variation across the time series
%         % strategy is to go through and eliminate the first one each time:
%         disp(['Removed time series with covariance greater than ' num2str(cov_thresh_r*100) ...
%             ' : from ' num2str(size(F,1)) ' to ' num2str(length(kr3))])
%         if ~isempty(kr3)
%             F=F(kr3,:);        % ********************* kr3 *******************************
%         else
%             disp(['No time series meet the ' num2str(thresh_cov_r) ' covariance criterion. Exiting.'])
%             return
%         end
%     else
%         disp('No filtering of time series based on high covariance')
%         kr3 = (1:size(F,1));
%     end
% 
% else
%     disp('No filtering of time series based on high covariance')
%     kr3 = (1:size(F,1));
% end
% 
% if cov_thresh_c>0
%     F_c=corrcoef(Fc); % I *think* this transpose is appropriate
%     [xi xj]=find(abs(F_c)>cov_thresh_c);
%     x=[xi xj];
%     r=setxor(1:size(x,1),find(x(:,1)>=x(:,2))); % do not include diagonal or below-diagonal entries
%     if ~isempty(r) % some rows are highly correlated
%         x=x(r,:);
%         kc3=setxor(1:size(F,2),unique(x(:,1)));
%         % have a mx2 list of indicies pairing rows with common variation across the time series
%         % strategy is to go through and eliminate the first one each time:
%         disp(['Removed time series with covariance greater than ' num2str(cov_thresh_c*100) ...
%             ' : from ' num2str(size(F,2)) ' to ' num2str(length(kc3))])
%         if ~isempty(kc3)
%             F=F(:,kc3);    % ********************* KC3 ***************************
%         else
%             disp(['No metrics meet the ' num2str(thresh_cov_c) ' covariance criterion. Exiting.'])
%             return
%         end
%     else
%         disp('No filtering of metrics based on high covariance')
%         kc3=1:size(F,2);
%     end
% 
% else
%     disp('No filtering of metrics based on high covariance')
%     kc3 = (1:size(F,2));
% end



%% Update the labels and stores
% Time series
kr_tot = kr0(kr1(kr2)); ntsn = length(kr_tot);

ts_ids_keepn = ts_ids_keep(kr_tot);
tsfn = tsf(kr_tot);
tskwn = tskw(kr_tot);
tsln = tsl(kr_tot);
% tsfn=tsf(kr_tot); tsln=tsl(kr_tot); tskwn=tskw(kr_tot);
% tsprepn=tsprep(kr_tot); tsmapn=tsmap(kr_tot);

if ~isempty(trainset)
    trainset = intersect(trainset,kr_tot); % re-adjust training indices, too
end

% Metrics
kc_tot = kc0(kc1(kc2)); nmn = length(kc_tot);

m_ids_keepn = m_ids_keep(kc_tot);
mcoden = mcode(kc_tot);
mlabn = mlab(kc_tot);
mkwn = mkw(kc_tot);
% mpointn = mpoint(kc_tot);
mlinkn = mlink(kc_tot);


% Can't change masters -- in an ideal world, you would check to see if any masters are never pointed to (i.e., not in mlink)
% and recalibrate the indexing, but I'm not going to bother.


% m_ids_keep, ts_ids_keep, tsf, tskw, tsl, mcode, mlab, mkw, mpoint, mlink, Mmid, Mmlab, Mmcode
fprintf(1,'We now have %u time series and %u operations in play\n',ntsn,nmn)
fprintf(1,'%u bad entries (%4.2f%%) in the %ux%u data matrix\n',sum(isnan(F(:))),sum(isnan(F(:)))/length(F(:))*100,size(F,1),size(F,2))

%% NORMALIZE THE LITTLE SHIT

if ismember(normopt,{'nothing','none'})
    fprintf(1,'NO NORMALIZING ON MY WATCH!\n')
else
    if isempty(trainset)
        % no training subset
        fprintf(1,'Normalizing a %u x %u object. Your patience is greatly appreciated...\n',ntsn,nmn)
        F = SUB_normalizematrix(F,normopt);
    else
        % retrieve a subset
        disp(['Normalizing a ' num2str(ntsn) 'x' num2str(nmn) ' object using ' num2str(length(trainset)) ...
            ' to train the transformation! Your patience is greatly appreciated...'])
        F = SUB_normalizematrix(F,normopt,trainset);
    end
    fprintf(1,'Normalized: the data matrix contains %u special-valued elements\n',sum(isnan(F(:))))
end

%% Remove bad entries
% these can be created by feature vectors that are constant after e.g., the
% sigmoid transform -- a bit of a weird thing to do if pre-filtering by 
% percentage...
nancol = zeros(size(F,2),1); %column of all NaNs
for i = 1:size(F,2)
    nancol(i) = all(isnan(F(:,i)));
end
if all(nancol) % all columns are NaNs
    error('After normalization, all columns were bad-values... :(');
elseif any(nancol) % there are columns that are all NaNs
    kc = find(nancol==0);
    F = F(:,kc);
    TS_loc_q = TS_loc_q(:,kc);
    fprintf(1,'We just removed (%u/%u) all-NaN columns from after normalization\n',sum(nancol),length(nancol));
    m_ids_keepn = m_ids_keepn(kc);
    mcoden = mcoden(kc);
    mlabn = mlabn(kc);
    mkwn = mkwn(kc);
    % mpointn = mpointn(kc);
    mlinkn = mlinkn(kc);
end



%% Now, make sure the columns are still good
% check again for constant columns after normalization
kc = find(range(F) ~= 0); % (NaN or positive)
if ~isempty(kc) && length(kc) < size(F,2)
    fprintf(1,'Post-normalization filtering of operations with constant outputs: from %u to %u', ...
                        size(F,2),length(kc))
    m_ids_keepn = m_ids_keepn(kc);
    mcoden = mcoden(kc);
    mlabn = mlabn(kc);
    mkwn = mkwn(kc);
    % mpointn = mpointn(kc);
    mlinkn = mlinkn(kc);
    F = F(:,kc);
    TS_loc_q = TS_loc_q(:,kc);
end

fprintf(1,'%u bad entries (%4.2f%%) in the %ux%u data matrix\n',sum(isnan(F(:))),sum(isnan(F(:)))/length(F(:))*100,size(F,1),size(F,2))

%% Done -- save output to file
fprintf(1,'Saving the trimmed, normalized data to local files: ')
TS_loc_N = F;
save('TS_loc_N.mat','TS_loc_N'); fprintf(1,'TS_loc_N')
% Save qualities in normalized matrix TS_loc_q_N
TS_loc_q_N = TS_loc_q; clear TS_loc_q;
save('TS_loc_q_N.mat','TS_loc_q_N'); fprintf(1,', TS_loc_q_N');

% contains all the clustered time series and metric information, respectively.
% Note that the clustering is only done on 'good' metrics and so nmcl<=nm
save('TS_loc_guides_N.mat','m_ids_keepn','ts_ids_keepn','tsfn','tskwn','tsln','mcoden','mlabn','mkwn','mlinkn','Mmid','Mmlab','Mmcode','-v7.3');
fprintf(1,', TS_loc_guides_N. All saved!\n');

end