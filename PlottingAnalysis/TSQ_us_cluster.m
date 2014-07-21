% ------------------------------------------------------------------------------
% TSQ_us_cluster
% ------------------------------------------------------------------------------
% 
% Perform unsupervised clustering on a matrix using a given method.
% 
% Loads a data matrix and clustering options and outputs a clustering of the
% indicies of this data matrix.
% 
% Quite alot of unnecessary baggage lives in this code. Apologies for that.
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

function [ackwgs, acgi, TS_DataMat_cl] = TSQ_us_cluster(norcl,ClusterMethod,ClusterParams)

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
% 1) norcl: can be the data matrix, or a string:
%                  'norm' -- loads HCTSA_N,
%                  'cl'   -- loads HCTSA_cl,
if nargin < 1,
    norcl = ''; % not necessary
end

% 2) ClusterMethod: a string specifying the type of clustering method to use
if nargin < 2 || isempty(ClusterMethod),
    ClusterMethod = 'linkage';
end

% 3) ClusterParams: specify the parameters for the clustering method
if nargin < 3
    ClusterParams = {}; % Defaults are specified within each method
end

% ------------------------------------------------------------------------------
%% Get the data
% ------------------------------------------------------------------------------
if ischar(norcl)
    switch norcl
    case 'cl'
        fprintf(1,'Loading HCTSA_cl.mat...');
        load('HCTSA_cl.mat','TS_DataMat')
        fprintf(1,' Loaded.\n');
    case 'norm'
        fprintf(1,'Loading HCTSA_N.mat...');
        load('HCTSA_N.mat','TS_DataMat')
        fprintf(1,' Loaded.\n');
    end
else % Input is a matrix to be clustered -- call it TS_DataMat.
    TS_DataMat = norcl;
end

% ------------------------------------------------------------------------------
%% Do the unsupervised clustering
% ------------------------------------------------------------------------------
switch ClusterMethod
	case 'linkage'
% 		fprintf(1,'Using inbuilt matlab linkage clustering\n');
        % parameter is a cell:
        % {DistanceMetric, LinkageMethod, showdend, clustth, savetofile, customR}
        % Better to make this a structure in future...
        %% Check inputs
        % ** DistanceMetric
        if length(ClusterParams)>=1 && ~isempty(ClusterParams{1})
            DistanceMetric = ClusterParams{1};
        else
            DistanceMetric = 'euclidean';
        end
        
        % ** LinkageMethod
        if length(ClusterParams)>=2 && ~isempty(ClusterParams{2})
            LinkageMethod = ClusterParams{2};
        else
            LinkageMethod = 'average';
        end
        
        fprintf(1,'Using %s linkage clustering on %s distances\n',LinkageMethod,DistanceMetric);
        
        % ** showdend
        if length(ClusterParams)>=3 && ~isempty(ClusterParams{3})
            showdend = ClusterParams{3};
        else
            showdend = 0;
        end
        if showdend == 0
            fprintf(1,'Suppressing dendrogram output\n')
        end
        
        % ** clustthresh -- many ways of doing this -- current way searches
        % until number of clusters is less than some threshold (but
        % obtained using inconsistent criterion). Another way is to just
        % specify a number of clusters and use the distance criterion
        % (i.e., just snips the dendrogram off at some threshold)...
        if length(ClusterParams) >= 4 && ~isempty(ClusterParams{4})
            clustth = ClusterParams{4}; % method (string), max # clusters (integer)
            % e.g., {'cutoff',10} % will get (max) 10 clusters using cutoff method
            % e.g., {'maxnclust',10} will get 10 clusters using distance criterion
        else
            clustth = []; % don't do clustering
        end
        if ~isempty(clustth)
            clusterM = clustth{1}; % METHOD for forming clusters
            clusterN = clustth{2}; % number of clusters to form
        end
        
        % ** savetofile -- if want to save output of linkage clustering to
        % file, to be read in later by other routines.
        % if 1, saves to file; if 2, loads from file, if zero, doesn't do
        % either -- (outputs only to function outputs).
        % can also specify a filename to load from as savetofile
        if length(ClusterParams)>=5 && ~isempty(ClusterParams{5})
            savetofile = ClusterParams{5};
        else
            savetofile = 0; % don't save output of linkage clustering to file
        end
        
		% output is just one group with a given ordering
		
		%% Linkage
        if isstruct(savetofile) % specify R and links within a structure
            R = savetofile.R;
            links = savetofile.links;
            clear savetofile
        elseif ischar(savetofile) % specify a filename containing R and links
            FileName = savetofile;
            load(FileName,'R','links'); % load custom distance and linkage information
        elseif (savetofile == 2) % this means to load from file
            if strcmp(TSorOps,'ts')
                FileName = 'TS_guide_clinksr.mat';
            else
                FileName = 'TS_guide_clinksc.mat';
            end
            load(FileName,'R','links');
        else
            % Pairwise distances
            if strcmp(DistanceMetric,'abscorr') % custom distance function
                if any(isnan(TS_DataMat(:)));
                    fprintf(1,'NaNs found in the input matrix. Distance calculations will probably be SLOW...\n')
                    % Use BF_pdist to calculate distances even when NaNs are present
                    R = BF_pdist(TS_DataMat,'corr',1);
                else % all good values -- can do this using pdist which is very fast
                    R = pdist(TS_DataMat,'corr');
                end
                R = 1 - abs(1-R);
                R(R < 0) = 0; % Sometimes get numerical error putting entries slightly under 0... (dirty fix but ok)
                fprintf(1,'abscorr transformation :: R between %f (0) -- %f (1)',min(R),max(R))
            else
                if any(isnan(TS_DataMat(:))) % NaNs: need to do this the slow way:
                    fprintf(1,'NaNs found in the input matrix. Distance calculations will probably be SLOW...\n')
                    % Use BF_pdist to calculate distances even when NaNs are present
                    R = BF_pdist(TS_DataMat,DistanceMetric,1);
                else
                    R = pdist(TS_DataMat,DistanceMetric);
                end
            end
            
            % links
            links = linkage(R,LinkageMethod);

            % Display cophentic correlation (goodness of linkage)
            cpc = cophenet(links,R);
            fprintf(1,'FYI, the cophenetic correlation is %f.\n',cpc)

            % Save to file
            if (savetofile == 1)
                if strcmp(TSorOps,'ts')
                    FileName = 'TS_guide_clinksr.mat';
                else
                    FileName = 'TS_guide_clinksc.mat';
                end
                fprintf(1,'Saving the linkage information as ''%s''...',FileName)
                save(FileName,'R','links','-v7.3')
                fprintf(1,' Done.\n');
                if nargout < 1
                    return % don't bother doing the rest if all we wanted was this
                end
            end
        end

        
        % ------------------------------------------------------------------------------
        %% Cluster
        % ------------------------------------------------------------------------------
        % extracts a discrete clustering from the hierarchy obtained above
		if ~isempty(clustth) % do this clustering
            switch clusterM
            case 'cutoffN'
                % specify number of clusters by inconsistent measure
                depth = 2; % depth down hierarchy to look
                criterion = 'inconsistent';
                
                cr = (0.1:0.1:10);
                nc = length(cr);
                nclusters = zeros(nc,1);
                for i = 1:nc
                    c = cr(i);
                    T = cluster(links,'cutoff',c,'depth',depth,'criterion',criterion);
                    nclusters(i) = max(T);
                    if nclusters(i) <= clusterN % we've got it!
                        break
                    end
                end
                
                % plot the result
                figure('color','w'); box('on');
                plot(cr(1:i),nclusters(1:i),'.-k')
                xlabel('cutoff value, c')
                ylabel('# clusters')
                nclusters = nclusters(i);

                % we've reached our cluster threshold! (or the end of cr)
                % acgi contains indices for members of each cluster
                disp(['Clustering at cutoff ' num2str(c) ' with ' num2str(nclusters) ' clusters'])
                
            case 'cutoff'
                % just do cutoff clustering
                % for this method clusterN is the cutoff value, rather than
                % the actual number of clusters.
                
                depth = 2; % depth down hierarchy to look
                criterion = 'inconsistent';
                T = cluster(links,'cutoff',clusterN,'depth',depth,'criterion',criterion);
                nclusters = max(T);
                
            case 'maxnclust'
                % just do distance-based clustering
                T = cluster(links,'maxclust',clusterN);
                nclusters = max(T);
                fprintf(1,'Distance-based clustering with %u clusters\n',nclusters)
            otherwise
                error('Unknown clustering method ''%s''',clusterM);
            end
            
            acgi = cell(nclusters,1);
            ackwgs = cell(nclusters,1);
            gil = zeros(nclusters,1);
            R = squareform(R); % could be fancy to save memory, but I think we can handle it...
            for j = 1:nclusters
                ackwgs{j} = ['AGG_C' num2str(j)];
                acgi{j} = find(T==j);
                % reorder in terms to put members closest to 'cluster
                % centre' (chosen by *mean* of group's feature vectors) first
                if length(acgi{j}) > 1 % more than one member
                    % Cluster center has minimum mean distance to other points:
                    [~,ix] = sort(sum(R(acgi{j},acgi{j})),'ascend');
                    acgi{j} = acgi{j}(ix);
                end
            end
            
            % Reorder by decreasing cluster size
            ClusterSize = cellfun(@length,acgi);
            [~,ix] = sort(ClusterSize,'descend');
            ackwgs = ackwgs(ix);
            acgi = acgi(ix);
            
        else % Don't do agglomerative clustering, just return the dendrogram ordering
            figure('color','w');
            if ~showdend, set(gcf,'Visible','off'); end % suppress figure output
            if size(TS_DataMat,1) < 1000 % small enough to try optimalleaforder
                try
%                     ord = bensoptimalleaforder(links,R); % NEW!
                    ord = optimalleaforder(links,R); % NEW!
                    [~,~,ord] = dendrogram(links,0,'r',ord);
                    fprintf(1,'Used optimalleaforder!\n')
                catch
                    beep
                    fprintf(1,'optimalleaforder was not used :(\n')
                    [~,~,ord] = dendrogram(links,0);
                end
            else
                fprintf(1,'Too big for optimalleaforder\n')
                [~,~,ord] = dendrogram(links,0);
            end
            ackwgs = {[LinkageMethod '_' DistanceMetric '_linkage']};
            acgi = ord; % outputs one cluster with an ordering given by the linkage clustering
            if ~showdend, close; end
        end
        
	case 'kmeans_matlab'
        %% Check the inputs
        disp('Using matlab''s statistics toolbox kmeans clustering')
        % ClusterParams specifies {k,distancemeasure,nrep,starts}
        % ** k
        if ~isempty(ClusterParams) && ~iscell(ClusterParams), ClusterParams = {ClusterParams}; end

        if length(ClusterParams)>=1
            k = ClusterParams{1};
        else
            k = 2;
            disp('forming 2 clusters using kmeans')
        end

        % ** distance measure
        if length(ClusterParams)>=2
            distancemeasure = ClusterParams{2};
            if strcmp(distancemeasure,'Euclidean')
                distancemeasure = 'sqEuclidean';
            end
        else
            distancemeasure = 'sqEuclidean';
            disp('using euclidean distance for kmeans')
        end

        % ** nrep
        % number of replicates
        if length(ClusterParams)>=3
            nrep = ClusterParams{3};
        else
            nrep = 1;
        end
		
		% ** starts
		% how to initialize the algorithm
		if length(ClusterParams)>=4
			starts = ClusterParams{4};
		else
			starts = 'sample'; % samples from the data matrix
        end

        
        %% Specify the model
        [idx,~,~,D] = kmeans(TS_DataMat, k, 'dist',distancemeasure, 'replicates',nrep,...
        							'start',starts, 'emptyaction','singleton', 'display','off');

        ackwgs = cell(k,1); % keywords
		acgi = cell(k,1); % indicies
        for i=1:k
            ackwgs{i} = ['KMEANS_C' num2str(i)];
            acgi{i} = find(idx==i);
		    d = D(acgi{i},i); % distances of points in this cluster to this cluster's centroid
		    [~,ix] = sort(d,'ascend');
		    acgi{i} = acgi{i}(ix); % those that are closest to the cluster centroid are listed first
        end
        
    case 'kmedoids'
        %% Check the inputs
        disp('Using Ben''s implementation of kmedoids')
        % ClusterParams specifies k, distancemeasure, nrep
        % ** k
        if isfield(ClusterParams,'k')
            k = ClusterParams.k;
        else
            k = 2;
            disp('Forming 2 clusters using kmediods')
        end
        if isfield(ClusterParams,'dmth')
            DistanceMetric = ClusterParams.dmth;
        else
            DistanceMetric = 'Euclidean';
        end
        if isfield(ClusterParams,'maxIter')
            maxIter = ClusterParams.maxIter;
        else
            maxIter = 50;
        end
        if isfield(ClusterParams,'nrep')
            nrep = ClusterParams.nrep;
        else
            nrep = 20;
        end
        if isfield(ClusterParams,'file')
            whatwithfile = ClusterParams.file; % filename to retrieve, or integer to specify
        else
            whatwithfile = 0; % calculate distance matrix now, don't save to file
        end
        if isfield(ClusterParams,'errmeas')
            errmeas = ClusterParams.errmeas;
        else
            errmeas = 'sum';
        end
        
        % retrieve distances if necessary
        if ischar(whatwithfile) % load from custom filename
            FileName = whatwithfile;
            LoadTimer = tic;
            load(FileName,'R'); % load custom distance ('links' not needed)
            fprintf(1,'Loaded R from %s in %s\n',FileName,BF_thetime(toc(LoadTimer)));
            clear LoadTimer
        elseif (whatwithfile == 2) % load from TS_guide_clinks
            if strcmp(TSorOps,'ts')
                FileName = 'TS_guide_clinksr.mat';
            else
                FileName = 'TS_guide_clinksc.mat';
            end
            load(FileName,'R'); % don't need 'links'
        elseif ndims(whatwithfile)==2 && length(whatwithfile)>1 % you've supplied R in file field!
            R = whatwithfile; clear whatwithfile
        else
            % Pairwise distances not provided: calculate them now
            if strcmp(DistanceMetric,'abscorr') % special distance function
                R = pdist(TS_DataMat,'correlation');
                R = 1-abs(1-R);
                R(R<0) = 0;% sometimes get numerical error
                fprintf(1,'R between %4.2f (0) - %4.2f (1)',min(R),max(R))
            else
                R = pdist(TS_DataMat,DistanceMetric);
            end
        end

        %% Run the algorithm
        if size(R,1)==1
            R = squareform(R); % we want a full matrix for benkmedoids
        end
        
        [~,~,errs,acgi] = benkmedoids(R,k,maxIter,nrep,errmeas);
        disp(['k-Medoids successful, with ' num2str(sum(errs)) ' ''error'''])
        % output acgi is already ordered by sum of distances to cluster
        % centre
        
        ackwgs = cell(k,1); % keywords
        for i=1:k
            ackwgs{i} = ['KMEDOIDS_C' num2str(i)];
        end

    case 'gmm'
        %% Gaussian Mixture Modeling
        % use Gaussian Mixture modeling from the Statistics Toolbox
        % Ben Fulcher 6/7/2010
        
        if ~isempty(ClusterParams) && ~iscell(ClusterParams), ClusterParams = {ClusterParams}; end

        if length(ClusterParams)>=1
            k = ClusterParams{1};
        else
            k = 2;
            disp('forming 2 clusters using a mixture of Gaussians')
        end
        
        options = statset('Display','final');
        gm = gmdistribution.fit(TS_DataMat, k, 'Options', options);
        % now assign clusters
        idx = cluster(gm,TS_DataMat);
        
%         for i=1:k
%             ackwgs{i} = ['KMEANS_C' num2str(i)];
%             acgi{i} = find(idx==i);
% 		    d = D(acgi{i},i); % distances of points in this cluster to this cluster's centroid
% 		    [d_s ix] = sort(d,'ascend');
% 		    acgi{i} = acgi{i}(ix); % those that are closest to the cluster centroid are listed first
%         end
        ackwgs = cell(k,1); % keywords
		acgi = cell(k,1); % indicies
        P = posterior(gm,TS_DataMat); % posterior under each mixture component
        for i = 1:k
            ackwgs{i} = ['GMM_C' num2str(i)];
            acgi{i} = find(idx==i);
            [~,ix] = sort(P(acgi{i},i),'ascend');
            acgi{i} = acgi{i}(ix);
            % now reorder based on distance to cluster centre (probably
            % better to do by posterior...?
%             cc = gm.mu(1,:); % the mean of this gaussian mixture
%             component
%             dd = zeros(length(acgi{i}),1);
%             for j = 1:length(acgi{i})
%                 dd(k) = sqrt(sum((F(acgi{i}(j),:)-cc).^2)); % euclidean distances
%             end
            
        end
        
        
    case 'spectral'
        %% Spectral Clustering: check the inputs
        disp('Using the spider package''s spectral clustering')
        % doesn't work so well on large numbers of features
        % ClusterParams specifies {k,sigma}
        % ** k, number of clusters
        if length(ClusterParams)>=1
            k = ClusterParams{1};
        else
            k = 2; % 2 clusters
        end
        
        % ** sigma, scale of exponential
        if length(ClusterParams)>=2
            sigma = ClusterParams{2};
        else
            sigma = 0.05;
        end
        
        %% Specify and train the model
        a = spectral;
        a.k = k;
        a.sigma = sigma;
        
        d = data(TS_DataMat);
        d.Y = []; % unsupervised
        
        [r,a] = train(a,d); % do the clustering
        
        %% Package Output
        ackwgs = cell(k,1);
        acgi = cell(k,1);
        
        for i=1:k
            ackwgs{i} = ['SPECTRAL_C' num2str(i)];
            acgi{i} = find(r.X==i); % unordered within each cluster
        end
	otherwise
		disp([ClusterMethod ' -- an invalid clustering option'])
		return
end


% Also output the clustered input matrix
if nargout > 2
	if iscell(acgi)
		TS_DataMat_cl = TS_DataMat(vertcat(acgi{:}),:);
	else
		TS_DataMat_cl = TS_DataMat(acgi,:);
	end
end


end