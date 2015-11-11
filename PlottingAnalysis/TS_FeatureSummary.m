function TS_FeatureSummary(opID, whatData, annotateParams)
% TS_FeatureSummary   How a given feature behaves across a time-series dataset
%
% Plots the distribution of outputs of an operation across the given dataset
% and allows the user to annotate time series onto the plot to visualize
% how the operation is behaving.
%
%---INPUTS:
% opID, the operation ID to plot
% whatData, the data to visualize (HCTSA.mat by default; cf. TS_LoadData)
% annotateParams, a structure of custom plotting options

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
% Check inputs
%-------------------------------------------------------------------------------
if nargin < 1
    opID = 1;
end

if nargin < 2 || isempty(whatData)
   whatData = 'raw'; % Visualize unnormalized outputs by default
end

if nargin < 3 || isempty(annotateParams) % annotation parameters
    annotateParams = struct(); % take 10 annotatations
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
theOp = ([Operations.ID]==opID);
dataVector = TS_DataMat(:,theOp); % the outputs of interest
notNaN = find(~isnan(dataVector));
dataVector = dataVector(notNaN); % remove bad values
TimeSeries = TimeSeries(notNaN); % remove bad values
theOperation = Operations(theOp);

if isempty(dataVector)
    error('No data for %s',Operations(theOp).Name);
end

% Retrieve group names also:
fileVarsStruct = whos('-file',whatDataFile);
fileVars = {fileVarsStruct.name};
if ismember('groupNames',fileVars)
    groupNames = load(whatDataFile,'groupNames');
    groupNames = groupNames.groupNames;
else
    groupNames = {};
end

%-------------------------------------------------------------------------------
% Sort out any custom plotting settings in the annotateParams structure
%-------------------------------------------------------------------------------
if ~isfield(annotateParams,'userInput')
    annotateParams.userInput = 1; % user clicks to annotate rather than randomly chosen
end
if ~isfield(annotateParams,'textAnnotation') % what text annotation to use
    annotateParams.textAnnotation = 'Name';
end
if ~isfield(annotateParams,'n')
    annotateParams.n = min(6,length(TimeSeries));
end
if ~isfield(annotateParams,'maxL')
    annotateParams.maxL = 100;
end

% if isfield(annotateParams,'whereann') % what text annotation to use
%     textann = annotateParams.whereann;
% else
%     whereann = 'onplot'; %'onplot', 'newplot'
% end

%-------------------------------------------------------------------------------
%% Plot the kernel smoothed density
%-------------------------------------------------------------------------------
fig = figure('color','w'); box('on'); hold on
% if strcmp(whereann,'newplot')
%     subplot(3,1,[1,2]); box('on'); hold on
% end

if isfield(TimeSeries,'Group')
    numGroups = length(unique([TimeSeries.Group]));
    annotateParams.groupColors = BF_getcmap('set1',numGroups,1);;

    % Repeat for each group
    fx = cell(numGroups,1);
    lineHandles = cell(numGroups+1,1);
    tsInd = cell(numGroups,1); % keeps track of indices from TimeSeries structure

    % Global distribution:
    [fr,xr,lineHandles{1}] = BF_plot_ks(dataVector,ones(1,3)*0.5,0,1,8);

    % Distribution for each group:
    for k = 1:numGroups
        [fr,xr,lineHandles{k+1}] = BF_plot_ks(dataVector([TimeSeries.Group]==k),...
                            annotateParams.groupColors{k},0,2,12);
        fx{k} = [xr',fr'];
        tsInd{k} = find([TimeSeries.Group]==k)';
    end
    xy = vertcat(fx{:});
    tsInd = vertcat(tsInd{:});
    TimeSeries = TimeSeries(tsInd);

    % Set up legend:
    legendText = cell(length(groupNames)+1,1);
    legendText{1} = 'combined';
    legendText(2:end) = groupNames;
    legend(horzcat(lineHandles{:}),legendText)

else
    % Just run a single global one
    [fr,xr] = BF_plot_ks(dataVector,'k',0,1.5,10);
    xy = [xr',fr'];
end

fig.Position = [fig.Position(1:2),649,354];

%-------------------------------------------------------------------------------
%% Annotate time series:
%-------------------------------------------------------------------------------
xlabel(theOperation.Name,'Interpreter','none');
ylabel('Probability Density')
BF_AnnotatePoints(xy,TimeSeries,annotateParams);
title(sprintf('[%u] %s (%u-sample annotations)',opID,theOperation.Name, ...
                annotateParams.maxL),'Interpreter','none');
xlabel('Outputs','Interpreter','none');

end
