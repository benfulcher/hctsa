function TS_LowDimInspect(whatData,whatAlgorithm,maxLength)
% LowDimGUI Plot time series in low-dimensional projection.
%
%---INPUTS:
% whatData: The HCTSA data file to use
% whatAlgorithm: The dimensionality-reduction algorithm to use ('pca','tsne')
% maxLength: maximum length of time series to display in the inspector

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

if nargin < 1 || isempty(whatData)
    whatData = 'norm'; % load in normalized data by default, from HCTSA_N.mat
end
if nargin < 2 || isempty(whatAlgorithm)
    whatAlgorithm = 'pca';
end
if nargin < 3
    maxLength = 1500;
end
%-------------------------------------------------------------------------------

%  Create and then hide the UI as it is being constructed.
f = figure('Color','w','Visible','off','Position',[360,500,420,450]);
h_TimeSeries = axes('Units','pixels','Position',[50,30,350,70]);
tHandle = text(0.2,0.5,'Click near a point above to inspect a time series!');
h_LowDim = axes('Units','pixels','Position',[50,150,350,290],...
                'ButtonDownFcn',{@annotateFigure_Callback},'HitTest','on');
align([h_TimeSeries,h_LowDim],'Left','None');

% Load data
[TS_DataMat,TimeSeries] = TS_LoadData(whatData);
[groupLabels,classLabels,groupLabelsInteger,numGroups] = TS_ExtractGroupLabels(TimeSeries);
groupColors = GiveMeColors(numGroups); % Set colors

%-------------------------------------------------------------------------------
% Do dimensionality reduction
[lowDimComponents,componentLabels] = BF_LowDimProject(TS_DataMat,whatAlgorithm,2);

%-------------------------------------------------------------------------------
% Plot two-dimensional projection
annotateParams = struct('n',0,'makeFigure',false);
[~,handles] = TS_Plot2d(lowDimComponents,TimeSeries,componentLabels,annotateParams,false,struct());
% Avoid hitting data points:
for i = 1:length(handles)
    handles{i}.HitTest = 'off';
end
% Title:
title(sprintf('Low-dimensional feature-based projection (%u time series; %u features)',...
                    size(TS_DataMat,1),size(TS_DataMat,2)))

%-------------------------------------------------------------------------------
% Initialize the UI
f.Units = 'normalized';
h_LowDim.Units = 'normalized';
h_TimeSeries.Units = 'normalized';
f.Name = 'Low-Dimensional Inspector';
movegui(f,'center')
f.Visible = 'on';

%-------------------------------------------------------------------------------
% function annotateButton_Callback(source,eventdata)
function annotateFigure_Callback(hObject,eventdata)
    % Bits and pieces from BF_AnnotatePoints(lowDimComponents,TimeSeries,annotateParams);

    if ~isempty(hObject.UserData)
        pC = hObject.UserData;
        delete(pC);
    else
        delete(tHandle)
    end

    xy_std = std(lowDimComponents);
    xy_mean = mean(lowDimComponents);
    xy_zscore = zscore(lowDimComponents);

    coordinates = get(h_LowDim,'CurrentPoint');
    point = coordinates(1,1:2);
    point_z = (point - xy_mean)./xy_std;
    iPlot = BF_ClosestPoint_ginput(xy_zscore,point_z);
    plotPoint = lowDimComponents(iPlot,:);
    theGroupIndex = groupLabelsInteger(iPlot);

    % Plot a circle around the annotated point:
    pC = plot(plotPoint(1),plotPoint(2),'o','MarkerSize',10,'MarkerEdgeColor',groupColors{theGroupIndex},...
                    'MarkerFaceColor',brighten(groupColors{theGroupIndex},0.5),'Parent',h_LowDim);
    hObject.UserData = pC;

    % Plot the time series in the inspector plot
    timeSeriesData = TimeSeries.Data{iPlot}; % (1:min(maxLength,end))
    plot(timeSeriesData,'-','color',groupColors{theGroupIndex},'LineWidth',2,...
                    'Parent',h_TimeSeries);
    h_TimeSeries.Title.String = sprintf('%s (%s)',...
                    TimeSeries.Name{iPlot},classLabels{theGroupIndex});
    h_TimeSeries.Title.Interpreter = 'none';
    h_TimeSeries.XLabel.String = 'Time';
    h_TimeSeries.XLim = [0,maxLength];
end

end
