function TS_FeatureSummary(opID, whatData, doViolin, annotateParams)
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

if nargin < 3 || isempty(doViolin) % annotation parameters
    doViolin = 0;
end

if nargin < 4 || isempty(annotateParams) % annotation parameters
    annotateParams = struct();
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
groupNames = TS_GetFromData(whatDataFile,'groupNames');
if isempty(groupNames)
    groupNames = {};
    timeSeriesGroup = [];
else
    timeSeriesGroup = [TimeSeries.Group]'; % Use group form
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
    annotateParams.n = min(15,length(TimeSeries));
end
if ~isfield(annotateParams,'maxL')
    annotateParams.maxL = 500;
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
    numGroups = length(unique(timeSeriesGroup));
    annotateParams.groupColors = BF_getcmap('set1',numGroups,1);
end

if doViolin
    % Violin plots

    rainbowColors = [BF_getcmap('set1',5,1); BF_getcmap('dark2',5,1)];

    if isfield(TimeSeries,'Group')
        dataCell = cell(numGroups+1,1);
        dataCell{1} = TS_DataMat(:,theOp); % global distribution
        for i = 1:numGroups
            dataCell{i+1} = (TS_DataMat(timeSeriesGroup==i,theOp));
        end

        myColors = cell(numGroups+1,1);
        myColors{1} = ones(3,1)*0.5; % gray for combined
        myColors(2:numGroups+1) = GiveMeColors(numGroups);
        extraParams = struct();
        extraParams.theColors = myColors;
        extraParams.customOffset = -0.5;
        extraParams.offsetRange = 0.7;

        ax = subplot(1,4,1:2);
        [ff,xx] = BF_JitteredParallelScatter(dataCell,1,1,0,extraParams);

        % Annotate lines for each feature in the distribution:
        [~,ix] = sort(TS_DataMat(:,theOp),'ascend');
        r = ix(round(linspace(1,length(ix),annotateParams.n)));
        for i = 1:annotateParams.n
            ri = find(xx{1}>=TS_DataMat(r(i),theOp),1);
            plot(0.5+0.35*[-ff{1}(ri),ff{1}(ri)],ones(2,1)*xx{1}(ri),'color',rainbowColors{rem(i-1,10)+1},'LineWidth',2)
            groupColor = myColors{1+timeSeriesGroup(r(i))};
            plot(0.5+0.35*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',groupColor,'MarkerEdgeColor',groupColor)
            plot(0.5-0.35*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',groupColor,'MarkerEdgeColor',groupColor)
        end
        ax.XTick = 0.5+(0:numGroups);
        axisLabels = cell(numGroups+1,1);
        axisLabels{1} = 'all';
        axisLabels(2:end) = groupNames;
        ax.XTickLabel = axisLabels;
    else
        % Just run a single global one
        dataCell = {TS_DataMat(:,theOp)};
        extraParams = struct();
        extraParams.theColors = {ones(3,1)*0.5};

        ax = subplot(1,4,1:2);
        [ff,xx] = BF_JitteredParallelScatter(dataCell,1,1,0,extraParams);

        % Annotate lines for each feature in the distribution:
        [~,ix] = sort(TS_DataMat(:,theOp),'ascend');
        r = ix(round(linspace(1,length(ix),annotateParams.n)));
        for i = 1:annotateParams.n
            ri = find(xx{1}>=TS_DataMat(r(i),theOp),1);
            rainbowColor = rainbowColors{rem(i-1,10)+1};
            plot(1+0.25*[-ff{1}(ri),ff{1}(ri)],ones(2,1)*xx{1}(ri),'color',rainbowColor,'LineWidth',2)
            plot(1+0.25*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',rainbowColor,'MarkerEdgeColor',rainbowColor)
            plot(1-0.25*ff{1}(ri),xx{1}(ri),'o','MarkerFaceColor',rainbowColor,'MarkerEdgeColor',rainbowColor)
        end
        ax.XTick = [];
    end
    title(sprintf('[%u]%s (%s)',theOperation.ID,theOperation.Name,theOperation.Keywords),'interpreter','none')
    ylabel('Feature value');

    % Time series annotations:
    ax = subplot(1,4,3:4);
    plotOptions.newFigure = 0;
    numIts = ceil(annotateParams.n/10);
    plotOptions.colorMap = cell(annotateParams.n,1);
    for i = 1:numIts
        rMap = (i-1)*10+1:i*10;
        rMap = rMap(rMap <= annotateParams.n);
        flipped = flipud(rainbowColors);
        plotOptions.colorMap(rMap) = flipped(1:length(rMap));
    end
    TS_plot_timeseries(whatData,annotateParams.n,fliplr([Operations(r).ID]),annotateParams.maxL,plotOptions);

    % Put rectangles if data is grouped
    if isfield(TimeSeries,'Group')
        rectHeight = 1/annotateParams.n;
        rectWidth = 0.1;
        for i = 1:annotateParams.n;
            rectangle('Position',[-rectWidth*1,(i-1)*rectHeight,rectWidth,rectHeight],...
                                    'FaceColor',myColors{1+timeSeriesGroup(r(i))});
        end
        ax.XLim = [-rectWidth,1];
    end

    fig.Position(3:4) = [1151,886];

else % kernel distributions
    if isfield(TimeSeries,'Group')

        % Repeat for each group
        fx = cell(numGroups,1);
        lineHandles = cell(numGroups+1,1);
        tsInd = cell(numGroups,1); % keeps track of indices from TimeSeries structure

        % Global distribution:
        [~,~,lineHandles{1}] = BF_plot_ks(dataVector,ones(1,3)*0.5,0,1,8);

        % Distribution for each group:
        for k = 1:numGroups
            [fr,xr,lineHandles{k+1}] = BF_plot_ks(dataVector(timeSeriesGroup==k),...
                                annotateParams.groupColors{k},0,2,12);
            fx{k} = [xr',fr'];
            tsInd{k} = find(timeSeriesGroup==k)';
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

end
