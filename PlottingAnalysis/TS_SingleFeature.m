function TS_SingleFeature(whatData,featID,makeViolin,makeNewFigure,whatStat,beVocal)
% TS_SingleFeature  Plot distributions for a single feature given a feature ID
%
%---INPUTS:
% whatData: the data to load in (cf. TS_LoadData)
% featID: the ID of the feature to plot
% makeViolin: makes a violin plot instead of overlapping kernel-smoothed distributions
% makeNewFigure: generates a new figure
% whatStat: can provide an already-computed stat for the feature (otherwise will
%           compute a simple linear classification based metric)

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs:
if nargin < 3
    makeViolin = 0;
end
if nargin < 4
    makeNewFigure = 0;
end
if nargin < 5
    whatStat = [];
end
if nargin < 6
    beVocal = 1;
end

%-------------------------------------------------------------------------------
% Load data:
[TS_DataMat,TimeSeries,Operations,whatDataSource] = TS_LoadData(whatData);
% Get groupNames if it exists:
groupNames = TS_GetFromData(whatData,'groupNames');
if isempty(groupNames)
    error('You must assign groups to data to use TS_SingleFeature. Use TS_LabelGroups.');
end

%-------------------------------------------------------------------------------
timeSeriesGroup = [TimeSeries.Group]'; % Use group form
numClasses = max(timeSeriesGroup);
op_ind = find([Operations.ID]==featID);

if isempty(op_ind)
    error('Operation with ID %u not found in %s',featID,whatDataSource);
end

if beVocal
    fprintf(1,'[%u] %s (%s)\n',featID,Operations(op_ind).Name,Operations(op_ind).Keywords);
end

%-------------------------------------------------------------------------------
% Plot this stuff:
if makeNewFigure
    f = figure('color','w');
end
hold on
ax = gca;
colors = GiveMeColors(numClasses);

if makeViolin
    dataCell = cell(numClasses,1);
    for i = 1:numClasses
        dataCell{i} = (TS_DataMat(timeSeriesGroup==i,op_ind));
    end

    % Re-order groups by mean (excluding any NaNs, descending):
    meanGroup = cellfun(@nanmean,dataCell);
    [~,ix] = sort(meanGroup,'descend');

    extraParams = struct();
    extraParams.theColors = colors(ix);
    extraParams.customOffset = -0.5;
    extraParams.offsetRange = 0.7;
    BF_JitteredParallelScatter(dataCell(ix),1,1,0,extraParams);

    % Adjust appearance:
    ax = gca;
    ax.XLim = [0.5+extraParams.customOffset,numClasses+0.5+extraParams.customOffset];
    ax.XTick = extraParams.customOffset+(1:numClasses);
    ax.XTickLabel = groupNames(ix);
    ylabel('Output')
    ax.TickLabelInterpreter = 'none';
    if makeNewFigure
        f.Position(3:4) = [402,159];
    end

    % Annotate rectangles:
    BF_AnnotateRect('diaglinear',TS_DataMat(:,op_ind),timeSeriesGroup,numClasses,colors,ax,'left');

    % Trim y-limits (with 2% overreach)
    ax.YLim(1) = min(TS_DataMat(:,op_ind))-0.02*range(TS_DataMat(:,op_ind));
    ax.YLim(2) = max(TS_DataMat(:,op_ind))+0.02*range(TS_DataMat(:,op_ind));

else
    linePlots = cell(numClasses,1);
    for i = 1:numClasses
        featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
        [~,~,linePlots{i}] = BF_plot_ks(featVector,colors{i},0,2,20,1);
    end
    % Trim x-limits (with 2% overreach)
    ax.XLim(1) = min(TS_DataMat(:,op_ind))-0.02*range(TS_DataMat(:,op_ind));
    ax.XLim(2) = max(TS_DataMat(:,op_ind))+0.02*range(TS_DataMat(:,op_ind));

    % Add a legend:
    legend([linePlots{:}],groupNames,'interpreter','none','Location','best')
    ylabel('Probability density')

    % Annotate rectangles:
    BF_AnnotateRect('diaglinear',TS_DataMat(:,op_ind),timeSeriesGroup,numClasses,colors,ax,'under');

    % Add x-label:
    xlabel('Output')

    % Adjust position
    if makeNewFigure
        f.Position(3:4) = [405,179];
    end
end

%-------------------------------------------------------------------------------
% Get cross-validated accuracy for this single feature using a Naive Bayes linear classifier:
if isempty(whatStat)
    numFolds = 10;
    accuracy = GiveMeCfn('diaglinear',TS_DataMat(:,op_ind),timeSeriesGroup,...
                            [],[],numClasses,1,[],[],numFolds);
    fprintf(1,'%u-fold cross validated balanced accuracy: %.2f +/- %.2f%%\n',...
                            numFolds,mean(accuracy),std(accuracy));
    statText = sprintf('%.1f%%',mean(accuracy));
else
    if isnumeric(whatStat)
        statText = sprintf('%.1f',whatStat);
    else % assume text otherwise
        statText = whatStat;
    end
end
title({sprintf('[%u] %s: %s',Operations(op_ind).ID,Operations(op_ind).Name,statText);...
                        ['(',Operations(op_ind).Keywords,')']},'interpreter','none')

end
