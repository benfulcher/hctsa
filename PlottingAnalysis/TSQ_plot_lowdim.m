% ------------------------------------------------------------------------------
% TSQ_plot_lowdim
% ------------------------------------------------------------------------------
%
% Calculates then plots a lower-dimensional feature-based representation of the
% data (e.g., using PCA).
% 
%---HISTORY:
% [Previously called TSQ_dimred]
% Ben Fulcher 31/3/2010 -- new classmeth option to specify classification
%                           method -- i.e., built in linear/quadratic; or
%                           svm-based method, etc.
% Ben Fulcher 18/4/2010 -- justus specifies whether to do PCA on the full matrix,
% 						   or just the groups specified in the given subset
% Ben Fulcher 13/7/2010 -- removed justus option!! Trying to clean up an
%                           inconsistency in labeling. Added TheData input
% Ben Fulcher 29/10/2010 -- added annotatep: can annotate time series
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

function TSQ_plot_lowdim(TheData,TsorOps,classmeth,showks,annotatep)

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(TheData)
    TheData = 'norm';
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

% if nargin < 5 || isempty(DimRedMethod)
%     DimRedMethod = 'pca';
%     fprintf(1,'Using PCA\n');
% end

if nargin < 3 || isempty(classmeth)
    classmeth = 'linclass';
    fprintf(1,'No discriminant\n'); 
end

if nargin < 4 || isempty(showks)
    showks = 1;
end

if nargin < 5 || isempty(annotatep)
    annotatep = struct('n',10);
end
if ~isstruct(annotatep)
    annotatep = struct('n',annotatep);
end

% ------------------------------------------------------------------------------
%% Load the data and group labeling from file
% ------------------------------------------------------------------------------
if strcmp(TheData,'cl') || strcmp(TheData,'norm') 
    % Retrive data from local files
    switch TheData
    case 'cl'
        TheDataFile = 'HCTSA_cl.mat';
    case 'norm'
        TheDataFile = 'HCTSA_N.mat';
    end
    fprintf(1,'Loading data and grouping information from %s...',TheDataFile);
    load(TheDataFile,'TS_DataMat');
    if strcmp(TsorOps,'ts')
        load(TheDataFile,'Operations')
        dimensionLabels = {Operations.Name}; clear Operations % We just need their names
        if isstruct(annotatep) || length(annotatep) > 1 || annotatep > 0
            load(TheDataFile,'TimeSeries')
            DataLabels = {TimeSeries.FileName};
            data_ids = [TimeSeries.ID];
            TimeSeriesData = {TimeSeries.Data};
            if isfield(TimeSeries,'Group')
                load(TheDataFile,'GroupNames')
                DataGroups = [TimeSeries.Group];
            else
                fprintf(1,'\n');
                error('No groups assigned -- Use TSQ_LabelGroups.')
            end
            clear('TimeSeries'); % we no longer need you
        end
    else
        load(TheDataFile,'TimeSeries')
        dimensionLabels = {TimeSeries.FileName}; clear TimeSeries
    end
    fprintf(1,' Loaded.\n');
else
    % The user provided data yourself
    if ~isfield(TheData,'DataMat')
        error('No field ''DataMat'' provided in the data input')
    elseif ~isfield(TheData,'Groups')
        error('No field ''Groups'' provided in the data input')
    end
    TS_DataMat = TheData.DataMat;
    DataGroups = TheData.Groups;
    if isfield(TheData,'DimLabels')
        dimensionLabels = TheData.DimLabels;
    else
        dimensionLabels = {};
    end
    if isfield(TheData,'GroupNames')
        GroupNames = TheData.GroupNames;
    else
        GroupNames = {};
    end
    if isfield(TheData,'DataLabels')
        DataLabels = TheData.DataLabels;
    else
        DataLabels = {};
    end
    if isfield(TheData,'TimeSeriesData')
        TimeSeriesData = TheData.TimeSeriesData;
    else
        TimeSeriesData = {};
    end
end

if strcmp(TsorOps,'ops')
    % Take the transpose of the input data matrix for operations
    TS_DataMat = TS_DataMat';
end


GroupIndices = BF_ToGroup(DataGroups);
% if isempty(gi)
%     gi = SUB_autoLabelQ(kwgs,TsorOps,TheData);
% end
% CheckEmpty = cellfun(@isempty,gi);
% if any(CheckEmpty)
%     error('No keywords found for: %s . Exiting.',kwgs{find(CheckEmpty,1)})
% end
% if (size(kwgs,2) == 2) && (size(kwgs,1) > 1); % specified subsets of each keyword
%    kwgs = kwgs(:,1);
% end
NumGroups = length(GroupIndices); % Number of groups

% ------------------------------------------------------------------------------
%% Do the dimensionality reduction
% ------------------------------------------------------------------------------

% Matlab's build-in PCA
% Sort it so that when choose different set of keywords the output is consistent
% There's a strange thing in princomp that can give different scores when the
% input rows are in a different order. The geometry is the same, but they're reflected
% relative to different orderings
fprintf(1,'Calculating principal components of the %u x %u data matrix...', ...
                    size(TS_DataMat,1),size(TS_DataMat,2));
[pc,score,latent] = princomp(TS_DataMat);
fprintf(1,' Done.\n');
percVar = round(latent/sum(latent)*1000)/10; % Percentage of variance explained (1 d.p.)

% Work out the main contributions to each principle component
featLabel = cell(2,2); % Feature label :: proportion of total [columns are PCs)
toLabel = cell(2,1);
ngcontr = min(length(dimensionLabels),2);  % the 2 greatest contributions to the first principle component
for i = 1:2
    [s1, ix1] = sort(abs(pc(:,i)),'descend');
    featLabel{i,1} = dimensionLabels(ix1(1:ngcontr));
    featLabel{i,2} = round(abs(s1(1:ngcontr)/sum(abs(s1)))*100)/100;
    for j = 1:ngcontr
        toLabel{i} = sprintf('%s%s (%f), ',toLabel{i},featLabel{i,1}{j},featLabel{i,2}(j));;
    end
    toLabel{i} = toLabel{i}(1:end-2);
end

% ------------------------------------------------------------------------------
% Plot this two-dimensional representation of the data using TSQ_plot_2d
% ------------------------------------------------------------------------------

NameString = 'PC';
for i = 1:2
    DataInfo.labels{i} = sprintf('%s%u (%f%%) : %s',NameString,i,percVar(i),toLabel{i});
end

% function TSQ_plot_2d(Features,DataInfo,TrainTest,annotatep,keepksdensities,lossmeth,extras)

DataInfo.GroupNames = GroupNames;
DataInfo.GroupIndices = GroupIndices;
DataInfo.DataLabels = DataLabels;
DataInfo.TimeSeriesData = TimeSeriesData;

TSQ_plot_2d(score(:,1:2),DataInfo,{},annotatep,showks,classmeth);

end