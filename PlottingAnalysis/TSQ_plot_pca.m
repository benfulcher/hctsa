% ------------------------------------------------------------------------------
% TSQ_plot_pca
% ------------------------------------------------------------------------------
%
% Calculates then plots a lower-dimensional feature-based representation of the
% data (e.g., using PCA).
% 
%---HISTORY:
% [Previously called TSQ_dimred]
% Ben Fulcher 31/3/2010 -- new classMeth option to specify classification
%                           method -- i.e., built in linear/quadratic; or
%                           svm-based method, etc.
% Ben Fulcher 18/4/2010 -- justus specifies whether to do PCA on the full matrix,
% 						   or just the groups specified in the given subset
% Ben Fulcher 13/7/2010 -- removed justus option!! Trying to clean up an
%                           inconsistency in labeling. Added whatData input
% Ben Fulcher 29/10/2010 -- added annotateParams: can annotate time series
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

function TSQ_plot_pca(whatData,TsorOps,showDist,classMeth,annotateParams)

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
    fprintf(1,'Getting data from HCTSA_N\n'); 
end

if nargin < 2 || isempty(TsorOps)
    fprintf(1,'Using for time series\n');
    TsorOps = 'ts';
end
if ~any(ismember(TsorOps,{'ops','ts'}));
    error('Specify either operations (''ops'') or time series (''ts'').');
end

% % Specify a keyword labeling of the data, kwgs:
% if nargin < 3
%     kwgs = {};
%     fprintf(1,'No assignment of the data? Ok up to you...\n');
% end
% 
% % Specify group indices, gi:
% if nargin < 4
%     gi = [];
%     fprintf(1,'Will obtain group indices from file\n');
% end

if nargin < 3 || isempty(showDist)
    showDist = 1;
end

if nargin < 4 || isempty(classMeth)
    classMeth = 'linclass';
    fprintf(1,'No discriminant\n'); 
end

if nargin < 5 || isempty(annotateParams)
    % Annotate 10 points by default
    fprintf(1,'Annotating 10 points by default with time series segment and filenames\n');
    annotateParams = struct('n',10,'textAnnotation','fileName');
end
if ~isstruct(annotateParams)
    annotateParams = struct('n',annotateParams);
end

% ------------------------------------------------------------------------------
%% Load the data and group labeling from file
% ------------------------------------------------------------------------------
if strcmp(whatData,'cl') || strcmp(whatData,'norm') 
    % Retrive data from local files
    switch whatData
    case 'cl'
        whatDataFile = 'HCTSA_cl.mat';
    case 'norm'
        whatDataFile = 'HCTSA_N.mat';
    end
    fprintf(1,'Loading data and grouping information from %s...',whatDataFile);
    load(whatDataFile,'TS_DataMat');
    if strcmp(TsorOps,'ts')
        load(whatDataFile,'Operations')
        dimensionLabels = {Operations.Name}; clear Operations % We just need their names
        if isstruct(annotateParams) || length(annotateParams) > 1 || annotateParams > 0
            load(whatDataFile,'TimeSeries')
            DataLabels = {TimeSeries.FileName};
            data_ids = [TimeSeries.ID];
            TimeSeriesData = {TimeSeries.Data};
            if isfield(TimeSeries,'Group')
                load(whatDataFile,'GroupNames')
                dataGroups = [TimeSeries.Group];
            else
                % No groups assigned:
                dataGroups = {};
                GroupNames = {'all data'};
                % fprintf(1,'\n');
                % error('No groups assigned -- Use TSQ_LabelGroups to provide group information.')
            end
            clear('TimeSeries'); % we no longer need you
        end
    else
        load(whatDataFile,'TimeSeries')
        dimensionLabels = {TimeSeries.FileName}; clear TimeSeries
    end
    fprintf(1,' Loaded.\n');
else
    % The user provided data yourself
    if ~isfield(whatData,'DataMat')
        error('No field ''DataMat'' provided in the data input')
    elseif ~isfield(whatData,'Groups')
        error('No field ''Groups'' provided in the data input')
    end
    TS_DataMat = whatData.DataMat;
    dataGroups = whatData.Groups;
    if isfield(whatData,'DimLabels')
        dimensionLabels = whatData.DimLabels;
    else
        dimensionLabels = {};
    end
    if isfield(whatData,'GroupNames')
        GroupNames = whatData.GroupNames;
    else
        GroupNames = {};
    end
    if isfield(whatData,'DataLabels')
        DataLabels = whatData.DataLabels;
    else
        DataLabels = {};
    end
    if isfield(whatData,'TimeSeriesData')
        TimeSeriesData = whatData.TimeSeriesData;
    else
        TimeSeriesData = {};
    end
end

if strcmp(TsorOps,'ops')
    % Take the transpose of the input data matrix for operations
    TS_DataMat = TS_DataMat';
end


% Label groups
if ~isempty(dataGroups)
    GroupIndices = BF_ToGroup(dataGroups);
else
    GroupIndices = {1:size(TS_DataMat,1)};
end
numGroups = length(GroupIndices); % Number of groups

% ------------------------------------------------------------------------------
%% Do the dimensionality reduction
% ------------------------------------------------------------------------------

% Matlab's build-in PCA
% Sort it so that when choose different set of keywords the output is consistent
% There's a strange thing in princomp that can give different scores when the
% input rows are in a different order. The geometry is the same, but they're reflected
% relative to different orderings

% Can't run PCA on data containing NaNs:
if any(isnan(TS_DataMat(:)))
    error('Data matrix contains NaNs');
end
fprintf(1,'Calculating principal components of the %u x %u data matrix...', ...
                    size(TS_DataMat,1),size(TS_DataMat,2));
[pc,score,latent] = princomp(TS_DataMat);
fprintf(1,' Done.\n');
percVar = round(latent/sum(latent)*1000)/10; % Percentage of variance explained (1 d.p.)

% Work out the main contributions to each principle component
featLabel = cell(2,2); % Feature label :: proportion of total (columns are PCs)
% toLabel = cell(2,1);
% topNcont = 1; %min(length(dimensionLabels),1);  % the top contribution to the first principle component
% for i = 1:2
%     [s1, ix1] = sort(abs(pc(:,i)),'descend');
%     featLabel{i,1} = dimensionLabels(ix1(1:topNcont));
%     featLabel{i,2} = round(abs(s1(1:topNcont)/sum(abs(s1)))*100)/100;
%     for j = 1:topNcont
%         toLabel{i} = sprintf('%s%s (%f), ',toLabel{i},featLabel{i,1}{j},featLabel{i,2}(j));;
%     end
%     toLabel{i} = toLabel{i}(1:end-2);
% end

% ------------------------------------------------------------------------------
% Plot this two-dimensional representation of the data using TSQ_plot_2d
% ------------------------------------------------------------------------------

nameString = 'PC';
for i = 1:2
    DataInfo.labels{i} = sprintf('%s %u (%.2f%%)',nameString,i,percVar(i));
end

% function TSQ_plot_2d(Features,DataInfo,TrainTest,annotateParams,keepksdensities,lossmeth,extras)

DataInfo.GroupNames = GroupNames;
DataInfo.GroupIndices = GroupIndices;
DataInfo.DataLabels = DataLabels;
DataInfo.TimeSeriesData = TimeSeriesData;

TSQ_plot_2d(score(:,1:2),DataInfo,{},annotateParams,showDist,classMeth);

end