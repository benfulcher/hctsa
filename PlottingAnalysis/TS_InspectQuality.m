function hadProblem = TS_InspectQuality(inspectWhat,customFile)
% TS_InspectQuality   Statistics of quality labels from an hctsa analysis.
%
% This function loads the calculation quality information from HCTSA.mat,
% and plots a visualization of where different special-valued outputs are occurring.
%
% Useful for checking where errors/special-valued outputs are occurring
%
%---INPUTS:
%
% inspectWhat: (i) 'summary' (default), summarize the proportion of each operation's
%                   outputs that correspond to each type of special-valued output
%              (ii) 'full' or 'all', show the full data matrix
%              (iii) 'reduced', only show operations that produce special-valued outputs
%              (iv) 'master', show master operations that produce special-valued outputs
%
% customFile: run on a custom HCTSA file (HCTSA.mat by default), can also be
%             strings like 'raw' or 'norm' (cf. TS_LoadData)

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

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(inspectWhat)
    inspectWhat = 'summary'; % only show operations that had at least one problem
end
if nargin < 2
    customFile = 'raw'; % use un-normalized data by default, 'HCTSA.mat'
end

% ------------------------------------------------------------------------------
% Load data:
% ------------------------------------------------------------------------------
[~,TimeSeries,Operations,whatDataFile] = TS_LoadData(customFile);
TS_Quality = TS_GetFromData(whatDataFile,'TS_Quality');
if isempty(TS_Quality)
    if ischar(whatDataFile)
        error('Quality labels, TS_Quality, not found in %s',whatDataFile);
    else
        error('Quality labels, TS_Quality, not found in input data structure');
    end
end
if all(isnan(TS_Quality(:)))
    if ischar(whatDataFile)
        error('No good quality labels in %s',whatDataFile);
    else
        error('No good quality labels in the input data structure');
    end
end
MasterOperations = TS_GetFromData(whatDataFile,'MasterOperations');
if isempty(MasterOperations)
    error('Must provide a data input containing MasterOperations');
end

% ------------------------------------------------------------------------------

switch inspectWhat
case {'full','all'}
    % Plot all quality labels (for all time series and all operations)

    % Check that the size is not too large to plot:
    checkSize(TS_Quality);

    % Plot:
    % Get handles for figure (f) and axes (ax):
    f = figure('color','w');
    ax = gca;

    % Give uncalculated entries the label 8
    TS_Quality(isnan(TS_Quality)) = 8;
    imagesc(TS_Quality)

    xlabel('Operations (op_id)','interpreter','none')
    ax.XTick = 1:length(Operations);
    ax.XTickLabel = [Operations.ID];
    ax.XTickLabelRotation = 90;

    formatYAxisColorBar;

case 'reduced'
    % First find where problems exist, and only show these columns
    qualityMean = nanmean(TS_Quality);
    hadProblem = (qualityMean > 0);

    if sum(hadProblem)==0
        fprintf(1,'No operations have problems! Nothing to inspect.\n');
        return
    end

    % Check that the size is not too large to plot:
    checkSize(TS_Quality(:,hadProblem));

    % Plot:
    % Get handles for figure (f) and axes (ax):
    f = figure('color','w');
    ax = gca;

    % Give uncalculated entries the label 8
    TS_Quality(isnan(TS_Quality)) = 8;

    imagesc(TS_Quality(:,hadProblem));

    ax.XTick = 1:sum(hadProblem);
    ax.XTickLabel = [Operations(hadProblem).ID];
    ax.XTickLabelRotation = 90;

    title(sprintf('Displaying %u x %u (displaying %u/%u operations with some special values',...
                    size(TS_Quality,1),sum(hadProblem),sum(hadProblem),size(TS_Quality,2)),...
                    'interpreter','none')

    formatYAxisColorBar;

case 'master'
    % Summarize at the level of master operations
    % Each row is now a different quality label

    numMasters = length(MasterOperations);
    opMIDs = [Operations.MasterID]; % save time by getting this first
    qualities = zeros(7,numMasters);
    for k = 1:numMasters
        masterID = [MasterOperations(k).ID];
        pointInd = (opMIDs==masterID); % indices of operations using this master function
        qualityLabels = TS_Quality(:,pointInd);
        qualityLabels = qualityLabels(qualityLabels > 0); % only want to look at special ones
        if isempty(qualityLabels) % no problems
            continue
        else
            % some problems -- work out which
            qualityLabels = unique(qualityLabels);
            qualities(:,k) = ismember(1:7,qualityLabels);
        end
    end

    % Focus on master operations with problems:
    hadProblem = (mean(qualities) > 0);

    % ------------
    % Text output:
    % ------------
    problemMops = MasterOperations(hadProblem);
    for i = 1:length(problemMops)
        fprintf('[%u] %s -- %s\n', problemMops(i).ID, ...
                problemMops(i).Label, problemMops(i).Code);
    end

    % ------------
    % Plotting:
    % ------------
    % Get handles for figure (f) and axes (ax):
    f = figure('color','w');
    ax = gca;
    imagesc(qualities(:,hadProblem));

    colormap([0.1686,0.5137,0.7294; 0.8431,0.0980,0.1098]);

    ax.YTick = 1:8;
    qLabels = {'error','NaN','Inf','-Inf','complex','empty','link error'};
    ax.YTickLabel = qLabels;

    ax.XTick = 1:sum(hadProblem);
    ax.XTickLabel = {MasterOperations(hadProblem).Code};
    ax.XTickLabelRotation = 90;

    % Get rid of tex interpreter format (for strings with underscores)
    set(ax,'TickLabelInterpreter','none');

    title('Master operations producing special-valued outputs.','interpreter','none')

case 'summary'
    % Summary as a line plot for operations that had some bad values

    % First find where problems exist, and only show these columns
    qualityMean = nanmean(TS_Quality,1);
    hadProblem = (qualityMean > 0);

    % Stop if there are no special values:
    if sum(hadProblem)==0
        fprintf(1,'No operations have problems! Nothing to inspect.\n');
        return
    else
        fprintf(1,'%u operations have special values:\n',sum(hadProblem));
    end

    % Resize TS_Quality
    TS_Quality = TS_Quality(:,hadProblem);

    % Give uncalculated entries the label 8
    TS_Quality(isnan(TS_Quality)) = 8;

    % Compute the proportion of each label per column:
    whatLabel = zeros(9,sum(hadProblem));
    for i = 1:9
        whatLabel(i,:) = mean(TS_Quality==i-1,1);
    end

    % Sort by the proportion of good values
    propGood = whatLabel(1,:);
    [~,ix] = sort(propGood,'ascend');
    whatLabel = whatLabel(:,ix);

    % ------------
    % Text output:
    % ------------
    problemOps = Operations(hadProblem);
    for i = 1:length(problemOps)
        ind = ix(i);
        fprintf('[%u] %s (%s) -- %.1f%% special values.\n', problemOps(ind).ID, problemOps(ind).Name, ...
                problemOps(ind).CodeString, 100*(1 - propGood(ind)));
    end

    % ------------
    % Plot:
    % ------------
    % Get handles for figure (f) and axes (ax):
    f = figure('color','w');
    ax = gca;

    bar(whatLabel','stacked');

    ax.XTick = 1:sum(hadProblem);
    unSortedTicks = [Operations(hadProblem).ID];
    ax.XTickLabel = unSortedTicks(ix);
    ax.XTickLabelRotation = 90;

    title(sprintf('Displaying %u/%u operations with some special values',...
                    sum(hadProblem),length(hadProblem)),...
                    'interpreter','none')

    formatYAxisColorBar(0,1);

    xlabel('Operations (op_id)','interpreter','none')
    ylabel(sprintf('Proportion of outputs (across %u time series)',size(TS_Quality,1)))

otherwise
    error('Unknown input option ''%s''',inspectWhat);
end

if nargout==0
    clear all
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function formatYAxisColorBar(doYaxis,offSet)
    % For visualizing where you have time series on y-axis and color shows what quality
    % label for each computation
    if nargin < 1
        doYaxis = 1;
    end
    if nargin < 2
        offSet = 0;
    end

    if doYaxis
        % Format the y axis
        ax.YTick = 1:length(TimeSeries);
        ylabel('Time series')
        ax.YTickLabel = {TimeSeries.Name};
    end

    % Get rid of tex interpreter format (for strings with underscores)
    set(ax,'TickLabelInterpreter','none');

    % Add a color bar:
    allLabels = {'good','error','NaN','Inf','-Inf','complex','empty','link error','(uncalculated)'};
    labelsExist = unique(TS_Quality(~isnan(TS_Quality)));
    maxShow = max(labelsExist); % maximum quality label to plot

    cb = colorbar(ax);
    caxis([-0.5+offSet,maxShow+0.5+offSet]);
    cb.Ticks = (0+offSet:1:maxShow+offSet);
    cb.TickLabels = allLabels(1:maxShow+1);

    % Set the colormap:
    allColors = [BF_getcmap('spectral',8,0); 0,0,0];
    allColors = allColors([8,1,2,3,4,5,6,7,9],:);
    colormap(allColors(1:maxShow+1,:))
end


% ------------------------------------------------------------------------------
function checkSize(A)
    % Checks size of a matrix, and sends a warning if too large
    if size(A,1)*size(A,2) > 100*9000
        fprintf('Attempting to plot each value in a large matrix (%ux%u)...\n',size(A,1),size(A,2))
        input('[Press any key to continue (or cntrl-C to abort)]','s');
    end
end


end
