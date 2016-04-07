function TS_SingleFeature(whatData,featID)
% TS_SingleFeature  Plot distributions for a single feature given an ID
%
%
%---INPUTS:
% whatData: the data to load in (cf. TS_LoadData)
% featID: the ID of the feature to plot

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Load data:
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
load(whatDataFile,'groupNames')

%-------------------------------------------------------------------------------
timeSeriesGroup = [TimeSeries.Group]'; % Use group form
numClasses = max(timeSeriesGroup);
op_ind = find([Operations.ID]==featID);

%-------------------------------------------------------------------------------
% Plot this stuff:
f = figure('color','w'); hold on
ax = gca;
colors = GiveMeColors(numClasses);

% Plot distributions first for the sake of the legend
linePlots = cell(numClasses,1);
for i = 1:numClasses
    featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
    [~,~,linePlots{i}] = BF_plot_ks(featVector,colors{i},0,2,20,1);
end
% Trim x-limits
ax.XLim(1) = min(TS_DataMat(:,op_ind))-0.02*range(TS_DataMat(:,op_ind));
ax.XLim(2) = max(TS_DataMat(:,op_ind))+0.02*range(TS_DataMat(:,op_ind));

% Add a legend:
legend([linePlots{:}],groupNames,'interpreter','none')
ylabel('Probability density')

% Annotate rectangles:
BF_AnnotateRect('fast_linear',TS_DataMat(:,op_ind),timeSeriesGroup,numClasses,colors,ax);

% Add x-label:
xlabel(Operations(op_ind).Name,'interpreter','none')

% Adjust position
f.Position(3:4) = [405,179];

% Get 10-fold cross-validated accuracy using a Naive Bayes linear classifier:
accuracy = GiveMeCfn('diaglinear',TS_DataMat(:,op_ind),timeSeriesGroup,[],[],numClasses,1,[],[],10);
fprintf(1,'10-fold cross validated balanced accuracy: %.2f +/- %.2f\n',mean(accuracy),std(accuracy));

end
