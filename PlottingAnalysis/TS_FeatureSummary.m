function TS_FeatureSummary(opID,whatData,doViolin,doInspect,annotateParams)
% TS_FeatureSummary   How a given feature behaves across a time-series dataset
%
% Plots the distribution of outputs of an operation across the given dataset
% and allows the user to annotate time series onto the plot to visualize
% how the operation is behaving.
%
%---INPUTS:
% opID, the operation ID to plot
% whatData, the data to visualize (HCTSA.mat by default; cf. TS_LoadData)
% doViolin, (logical) show distributions as a violin plot (instead of
%           conventional kernel-smoothed distributions).
% doInspect, (logical) whether to use an interactive version of the plot.
% annotateParams, a structure of custom plotting options.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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
if nargin < 3 || isempty(doViolin) % annotation parameters
    doViolin = true;
end
if nargin < 4
    doInspect = true;
end
if nargin < 5 || isempty(annotateParams) % annotation parameters
    annotateParams = struct();
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);
theOp = (Operations.ID==opID);
if ~any(theOp)
    error('No matches for operation ID %u',opID);
end

dataVector = TS_DataMat(:,theOp); % the outputs of interest
notNaN = find(~isnan(dataVector));
dataVector = dataVector(notNaN); % remove bad values
TimeSeries = TimeSeries(notNaN,:); % remove bad values
numTimeSeries = height(TimeSeries);
theOperation = table2struct(Operations(theOp,:));

if isempty(dataVector)
    error('No data for %s',theOperation.Name);
end

% Retrieve group names also:
[timeSeriesGroup,classLabels,groupLabelsInteger,numGroups] = TS_ExtractGroupLabels(TimeSeries);
annotateParams.groupColors = BF_GetColorMap('set1',numGroups,true);

%-------------------------------------------------------------------------------
% Apply default plotting settings in the annotateParams structure
if ~isfield(annotateParams,'n')
    annotateParams.n = min(15,height(TimeSeries));
end
if annotateParams.n < 1
    error('You need to specify at least one time series to annotate with TS_FeatureSummary');
end
if ~isfield(annotateParams,'maxL')
    annotateParams.maxL = 1000;
end

%-------------------------------------------------------------------------------
% Plot the distribution of feature-values across the dataset
%-------------------------------------------------------------------------------
f = figure('color','w');
box('on');
hold('on');

if doViolin
    %-------------------------------------------------------------------------------
    % Violin plot of a distribution for each labeled time-series class

    if ~doInspect
        % Spaced time-series annotations shown on the right of the plot.
        % Determine a subset, highlightInd, of time series to highlight:
        [~,ix] = sort(TS_DataMat(:,theOp),'ascend');
        highlightInd = ix(round(linspace(1,length(ix),annotateParams.n)));
        rainbowColors = [BF_GetColorMap('set1',5,1); BF_GetColorMap('dark2',5,1)];
    end

    % First the distribution(s):
    if doInspect
        h_violin = subplot(3,1,1:2);
        h_violin.ButtonDownFcn = {@annotateFigure_Callback};
        h_violin.HitTest = 'on';
    else
        h_violin = subplot(1,4,1:2);
    end

    % Set up parameters:
    extraParams = struct();
    if doInspect
        extraParams.makeHorizontal = true;
        extraParams.dontHitMe = true;
    end

    if numGroups > 1
        dataCell = cell(numGroups+1,1);
        dataCell{1} = TS_DataMat(:,theOp); % global distribution
        for i = 1:numGroups
            dataCell{i+1} = TS_DataMat(timeSeriesGroup==classLabels{i},theOp);
        end

        myColors = cell(numGroups+1,1);
        myColors{1} = ones(3,1)*0.5; % gray for combined
        myColors(2:numGroups+1) = GiveMeColors(numGroups);
        extraParams.theColors = myColors;
        extraParams.customOffset = -0.5;
        extraParams.offsetRange = 0.7;
        [ff,xx,xScatter,yScatter] = BF_JitteredParallelScatter(dataCell,true,true,false,extraParams);

        % Add lines denoting each annotated time-series in the distribution:
        if ~doInspect
            for i = 1:annotateParams.n
                ri = find(xx{1} >= TS_DataMat(highlightInd(i),theOp),1);
                groupColor = myColors{1 + groupLabelsInteger(highlightInd(i))};
                plot(0.5 + 0.35*[-ff{1}(ri),ff{1}(ri)],ones(2,1)*xx{1}(ri),'color',rainbowColors{rem(i-1,10) + 1},'LineWidth',2)
                plot(0.5 + 0.35*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',groupColor,'MarkerEdgeColor',groupColor)
                plot(0.5 - 0.35*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',groupColor,'MarkerEdgeColor',groupColor)
            end
        end

        % Label axes:
        axisLabels = cell(numGroups + 1,1);
        axisLabels{1} = 'all';
        axisLabels(2:end) = classLabels;
        if doInspect
            h_violin.YTick = 0.5 + (0:numGroups);
            h_violin.YTickLabel = axisLabels;
            h_violin.YTickLabelRotation = 20;
        else
            h_violin.XTick = 0.5 + (0:numGroups);
            h_violin.XTickLabel = axisLabels;
            h_violin.XTickLabelRotation = 20;
        end
    else
        % No groups: just show the global distribution
        extraParams.theColors = {ones(3,1)*0.5};
        [ff,xx,xScatter,yScatter] = BF_JitteredParallelScatter({TS_DataMat(:,theOp)},false,true,false,extraParams);

        % Annotate lines for each feature in the distribution:
        if ~doInspect
            for i = 1:annotateParams.n
                ri = find(xx{1} >= TS_DataMat(highlightInd(i),theOp),1);
                rainbowColor = rainbowColors{rem(i - 1,10)+1};
                plot(1 + 0.25*[-ff{1}(ri),ff{1}(ri)],ones(2,1)*xx{1}(ri),'color',rainbowColor,'LineWidth',2)
                plot(1 + 0.25*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',rainbowColor,'MarkerEdgeColor',rainbowColor)
                plot(1 - 0.25*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',rainbowColor,'MarkerEdgeColor',rainbowColor)
            end
            h_violin.XTick = [];
        end
    end

    h_violin.TickLabelInterpreter = 'none';
    title(sprintf('[%u]%s (%s)',theOperation.ID,theOperation.Name,theOperation.Keywords),...
                                'interpreter','none')

    if ~doInspect
        ylabel('Feature value');
    else
        xlabel('Feature value');
    end

    %-------------------------------------------------------------------------------
    % Time series annotations using TS_PlotTimeSeries
    % (cycling through groups of 10 rainbow colors):
    if doInspect
        h_TimeSeries = subplot(3,1,3);
        tHandle = text(0.2,0.5,'Click near a point above to inspect a time series!');
    else
        h_TimeSeries = subplot(1,4,3:4);
        plotOptions.newFigure = false;
        plotOptions.colorMap = cell(annotateParams.n,1);
        for i = 1:annotateParams.n
            plotOptions.colorMap{i} = rainbowColors{rem(i-1,10)+1};
        end
        plotOptions.colorMap = flipud(plotOptions.colorMap);

        TS_PlotTimeSeries(TimeSeries,annotateParams.n,flipud(highlightInd),annotateParams.maxL,plotOptions);

        % Put rectangles if data is grouped
        if ~isempty(timeSeriesGroup)
            rectHeight = 1/annotateParams.n;
            rectWidth = 0.1;
            for i = 1:annotateParams.n
                rectangle('Position',[-rectWidth*1,(i - 1)*rectHeight,rectWidth,rectHeight],...
                                'FaceColor',myColors{1 + groupLabelsInteger(highlightInd(i))});
            end
            h_TimeSeries.XLim = [-rectWidth,1];
        end
        f.Position(3:4) = [1151,886];
    end

else
    % Horizontal kernel-smoothed distribution(s)
    if numGroups > 1
        % Repeat for each group
        fx = cell(numGroups,1);
        lineHandles = cell(numGroups+1,1);
        tsInd = cell(numGroups,1); % keeps track of indices from TimeSeries structure

        % Global distribution:
        [~,~,lineHandles{1}] = BF_plot_ks(dataVector,ones(1,3)*0.5,0,1,8);

        % Distribution for each group:
        for k = 1:numGroups
            [fr,xr,lineHandles{k+1}] = BF_plot_ks(dataVector(timeSeriesGroup==classLabels{k}),...
                                annotateParams.groupColors{k},0,2,12);
            fx{k} = [xr',fr'];
            tsInd{k} = find(timeSeriesGroup==classLabels{k});
        end
        xy = vertcat(fx{:});
        % Now make sure that elements of TimeSeries matches ordering of xy
        tsInd = vertcat(tsInd{:});
        % ix = arrayfun(@(x)find(x==tsInd),1:height(TimeSeries));
        TimeSeries = TimeSeries(tsInd,:);

        % Set up legend:
        legendText = cell(length(classLabels)+1,1);
        legendText{1} = 'combined';
        legendText(2:end) = classLabels;
        legend(horzcat(lineHandles{:}),legendText)
    else
        % Just run a single global one (black)
        [fr,xr] = BF_plot_ks(dataVector,'k',0,1.5,10);
        xy = [xr',fr'];
    end

    f.Position = [f.Position(1:2),649,354];

    %-------------------------------------------------------------------------------
    % Annotate time series:
    xlabel(theOperation.Name,'Interpreter','none');
    ylabel('Probability Density')
    BF_AnnotatePoints(xy,TimeSeries,annotateParams);
    title(sprintf('[%u] %s (%u-sample annotations)',opID,theOperation.Name, ...
                    annotateParams.maxL),'Interpreter','none');
    xlabel('Outputs','Interpreter','none');
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
function annotateFigure_Callback(hObject,eventData)
    % Bits and pieces from BF_AnnotatePoints(lowDimComponents,TimeSeries,annotateParams);

    % Get the coordinates of the current cursor position
    coordinates = get(h_violin,'CurrentPoint');
    point = coordinates(1,1:2);

    if ~isempty(hObject.UserData)
        pC = hObject.UserData;
        delete(pC);
    else
        delete(tHandle)
    end
    if numGroups == 1
        xyData = [xScatter{1},yScatter{1}];
    else
        xyData = arrayfun(@(x)[xScatter{x},yScatter{x}],1:numGroups+1,'UniformOutput',false);
    end

    % Match the point
    if numGroups == 1
        iPlot = BF_ClosestPoint_ginput(xyData,point,true);
        plotPoint = xyData(iPlot,:);
    else
        iPlots = zeros(numGroups + 1,1);
        minDists = zeros(numGroups + 1,1);
        for g = 1:numGroups+1
            [iPlots(g),minDists(g)] = BF_ClosestPoint_ginput(xyData{g},point,true);
        end
        [~,theIndex] = min(minDists);
        plotPoint = xyData{theIndex}(iPlots(theIndex),:);
        if theIndex==1 % from the full distribution
            iPlot = iPlots(1);
        else % from a specific group
            % Find its index from TimeSeries
            theIndices = 1:numTimeSeries;
            theIndicesGroup = theIndices(timeSeriesGroup==classLabels{theIndex-1});
            iPlot = theIndicesGroup(iPlots(theIndex));
            % iPlot should be an index of the TimeSeries table
        end
    end
    if numGroups > 1
        theGroupIndex = groupLabelsInteger(iPlot);
        theColor = myColors{1 + theGroupIndex};
        if theGroupIndex==0
            % Can only happen in the rare case of a partially group-labeled dataset
            thisClassLabel = '<Unlabeled>';
        else
            thisClassLabel = classLabels{theGroupIndex};
        end
        titleText = sprintf('%s (%s) [%.3f]',TimeSeries.Name{iPlot},thisClassLabel,plotPoint(1));
    else
        theColor = zeros(1,3);
        titleText = sprintf('%s [%.3f]',TimeSeries.Name{iPlot},plotPoint(1));
    end

    % Plot a circle around the annotated point:
    pC = plot(plotPoint(1),plotPoint(2),'o','MarkerSize',10,'MarkerEdgeColor',theColor,...
                    'MarkerFaceColor',brighten(theColor,0.5),'Parent',h_violin);
    hObject.UserData = pC;

    % Plot the time series in the inspector plot
    timeSeriesData = TimeSeries.Data{iPlot};
    plot(timeSeriesData,'-','color',theColor,'LineWidth',2,'Parent',h_TimeSeries);
    h_TimeSeries.Title.Interpreter = 'none';
    h_TimeSeries.XLabel.String = 'Time';
    h_TimeSeries.XLim = [0,annotateParams.maxL];
    h_TimeSeries.Title.String = titleText;
end

end
