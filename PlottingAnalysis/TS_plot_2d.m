function TS_plot_2d(Features,TimeSeries,featureLabels,groupNames,annotateParams,showDistr,whatClassifier)
% TS_plot_2d   Plots a dataset in a two-dimensional space.
%
% e.g., The space of two chosen features, or two principal components.
%
%---INPUTS:
% Features, an Nx2 vector of where to plot each of the N data objects in the
%           two-dimensional space
%
% TimeSeries, structure array for time series metadata
%
% featureLabels, cell of labels for each feature
%
% groupNames, cell containing a label for each class of timeseries
%
% annotateParams, a structure containing all the information about how to annotate
%           data points. Fields can include:
%               - n, the number of data points to annotate
%               - userInput, 0: randomly selected datapoints, 1: user clicks to annotate datapoints
%               - fdim, 1x2 vector with width and height of time series as fraction of plot size
%               - maxL, maximum length of annotated time series
%               - textAnnotation: 'Name', 'tsid', or 'none' to annotate this data
%               - cmap, a cell of colors, with elements for each group
%               - theMarkerSize, a custom marker size
%               - theLineWidth: line width for annotated time series
%
% showDistr, if 1 (default), plots marginal density estimates for each variable
%                   (above and to the right of the plot), otherwise set to 0.
%
% whatClassifier, can select a classifier to fit to the different classes (e.g.,
%               'linclass' for a linear classifier).

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------

% Features should be a Nx2 vector of where to plot each of the N data objects in the two-dimensional space
if nargin < 1
    error('You must provide two-dimensional feature vectors for the data.')
end

if nargin < 3 || isempty(featureLabels)
    featureLabels = {'',''};
end

if nargin < 5 || isempty(annotateParams)
    annotateParams = struct('n',0); % don't annotate
end

% By default, plot kernel density estimates above and on the side of the plot:
if nargin < 6 || isempty(showDistr)
    showDistr = 1;
end

if nargin < 7 || isempty(whatClassifier)
    whatClassifier = 'svm_linear';
end

makeFigure = 1; % default is to plot on a brand new figure('color','w')

%-------------------------------------------------------------------------------
% Preliminaries
%-------------------------------------------------------------------------------
if isfield(TimeSeries,'Group')
    groupLabels = [TimeSeries.Group]'; % Convert GroupIndices to group form
    numClasses = length(unique(groupLabels));
else
    % No group information assigned to time series
    numClasses = 1;
    groupLabels = ones(length(TimeSeries),1);
end

% ------------------------------------------------------------------------------
%% Plot
% ------------------------------------------------------------------------------
if makeFigure % can set extras.makeFigure = 0 to plot within a given setting
    f = figure('color','w'); box('on'); % white figure
    f.Position = [f.Position(1), f.Position(2), 600, 550];
end

% Set colors
if (numClasses == 1)
    groupColors = {[0,0,0]}; % Just use black...
else
    if isstruct(annotateParams) && isfield(annotateParams,'cmap')
        if ischar(annotateParams.cmap)
            groupColors = BF_getcmap(annotateParams.cmap,numClasses,1);
        else
            groupColors = annotateParams.cmap; % specify the cell itself
        end
    else
        groupColors = GiveMeColors(numClasses);
    end
end
annotateParams.groupColors = groupColors;

% ------------------------------------------------------------------------------
%% Plot distributions
% ------------------------------------------------------------------------------
if showDistr
    % Top distribution (marginal of first feature)
    subplot(4,4,1:3); hold on; box('on')
    maxx = 0; minn = 100;
    for i = 1:numClasses
        fr = BF_plot_ks(Features(groupLabels==i,1),groupColors{i},0);
        maxx = max([maxx,fr]); minn = min([minn,fr]);
    end
    axTop = gca;
    set(axTop,'XTickLabel',[]);
    set(axTop,'YTickLabel',[]);
    set(axTop,'ylim',[minn,maxx]);

    % Side distribution (marginal of second feature)
    subplot(4,4,[8,12,16]); hold on; box('on')
    maxx = 0; minn = 100;
    for i = 1:numClasses
        fr = BF_plot_ks(Features(groupLabels==i,2),groupColors{i},1);
        maxx = max([maxx,fr]); minn = min([minn,fr]);
    end
    axSide = gca;
    set(axSide,'XTickLabel',[]);
    set(axSide,'YTickLabel',[]);
    set(axSide,'xlim',[minn,maxx]);
end

% ------------------------------------------------------------------------------
%% Set up a 2D plot
% ------------------------------------------------------------------------------
if showDistr
    subplot(4,4,[5:7,9:11,13:15]); box('on');
    axMain = gca;
end
hold on;

if isfield(annotateParams,'theMarkerSize');
    theMarkerSize = annotateParams.theMarkerSize; % specify custom marker size
else
    theMarkerSize = 12; % Marker size for '.'
end

for i = 1:numClasses
    plot(Features(groupLabels==i,1),Features(groupLabels==i,2),...
                '.','color',groupColors{i},'MarkerSize',theMarkerSize)
end

% Link axes
if showDistr
    linkaxes([axMain,axTop],'x');
    linkaxes([axMain,axSide],'y');
end

% ------------------------------------------------------------------------------
% Set Legend
%-------------------------------------------------------------------------------
if numClasses > 1
    legendText = cell(numClasses,1);
    for i = 1:numClasses
        if ~isempty(groupNames)
            legendText{i} = sprintf('%s (%u)',groupNames{i},sum(groupLabels==i));
        else
            legendText{i} = sprintf('Group %u (%u)',i,sum(groupLabels==i));
        end
    end
    legend(legendText,'interpreter','none');
end

%-------------------------------------------------------------------------------
% Annotate points:
%-------------------------------------------------------------------------------
% Label axes first without classification rates so the user can see what they're doing when annotating
labelAxes(0);
if isfield(TimeSeries,'Data')
    % Only attempt to annotate if time-series data is provided

    % Produce xy points
    xy = Features;

    % Go-go-go:
    BF_AnnotatePoints(xy,TimeSeries,annotateParams);
end

% ------------------------------------------------------------------------------
%% Do classification and plot a classify boundary?
% ------------------------------------------------------------------------------
didClassify = 0;
if numClasses > 1
    % Compute the in-sample classification rate:
    classRate = zeros(3,1); % classRate1, classRate2, classRateboth
    try
        fprintf(1,'Estimating %u-class classification rates for each feature (and in combination)...\n',numClasses);
        classRate(1) = GiveMeCfn(whatClassifier,Features(:,1),groupLabels,Features(:,1),groupLabels,numClasses);
        classRate(2) = GiveMeCfn(whatClassifier,Features(:,2),groupLabels,Features(:,2),groupLabels,numClasses);
        [classRate(3),~,whatLoss] = GiveMeCfn(whatClassifier,Features,groupLabels,Features(:,1:2),groupLabels,numClasses);
        % Record that classification was performed successfully:
        didClassify = 1;
        fprintf(1,'%s in 2-d space: %.2f%%\n',whatLoss,classRate(3));
    catch emsg
        fprintf(1,'\nLinear classification rates not computed\n(%s)\n',emsg.message);
        classRate(:) = NaN;
    end

    % Also plot an SVM classification boundary:
    if numClasses < 5
        fprintf(1,'Estimating classification boundaries...\n');
        try
            % Train the model (in-sample):
            [~,Mdl] = GiveMeCfn(whatClassifier,Features,groupLabels,Features,groupLabels,numClasses);

            % Predict scores over the 150x150 grid through space
            gridInc = 150;
            [x1Grid,x2Grid] = meshgrid(linspace(min(Features(:,1)),max(Features(:,1)),gridInc),...
                                       linspace(min(Features(:,2)),max(Features(:,2)),gridInc));
            fullGrid = [x1Grid(:),x2Grid(:)];
            predLabels = predict(Mdl,fullGrid);

            % For each class plot the contour of their region of the space:
            if numClasses==2
                contour(x1Grid,x2Grid,reshape(predLabels-1.5,size(x1Grid)),[0 0],'--k','LineWidth',2);
            else
                for i = 1:numClasses
                    % isMostProbable = scores(:,i) > max(scores(:,setxor(1:size(scores,2),i)),[],2);
                    isMostProbable = (predLabels==i);
                    if ~all(isMostProbable==isMostProbable(1));
                        contour(x1Grid,x2Grid,reshape(isMostProbable,size(x1Grid)),[0.5 0.5],'-','LineWidth',2,'color',groupColors{i});
                    end
                end
            end
        catch emsg
            warning('Error fitting classification model in 2-d space')
            keyboard
        end
    end
end

%-------------------------------------------------------------------------------
% Relabel axes with classification rates and set title if classification
% performed successfully:
%-------------------------------------------------------------------------------
if didClassify
    labelAxes(1);
    title(sprintf('Combined classification rate (%s) = %.2f%%',whatClassifier, ...
                    round(classRate(3,1))),'interpreter','none');
end

%-------------------------------------------------------------------------------
function labelAxes(didClassify)
    %-------------------------------------------------------------------------------
    %% Label Axes
    % ------------------------------------------------------------------------------
    if didClassify
        labelText = cell(2,1);
        for i = 1:2
            labelText{i} = sprintf('%s (acc = %.2f %%)',featureLabels{i}, ...
                                    round(classRate(i,1)));
        end
    else
        labelText = featureLabels;
    end

    xlabel(labelText{1},'interpreter','none')
    ylabel(labelText{2},'interpreter','none')
end

end
