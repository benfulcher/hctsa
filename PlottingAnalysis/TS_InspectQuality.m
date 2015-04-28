% --------------------------------------------------------------------------
% TS_InspectQuality
% --------------------------------------------------------------------------
% 
% This function loads the matrix in HCTSA_loc.mat, plots it, showing the 
% quality labels of each entry.
% 
% Most useful for checking where errors/special-valued outputs are occuring
% 
%---INPUTS:
%
% inspectWhat: 'full', show the full data matrix
%              ''
% 
% Ben Fulcher, 2015-03-29
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
function TS_InspectQuality(inspectWhat)

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(inspectWhat)
    inspectWhat = 'reduced'; % only show operations that had at least one problem
end

% ------------------------------------------------------------------------------
% Load data:
% ------------------------------------------------------------------------------
load('HCTSA_loc.mat','TS_Quality','TimeSeries','Operations','MasterOperations');

if ~exist('TS_Quality')
    error('Quality labels not found in HCTSA_loc.mat');
end

if ~any(~isnan(TS_Quality(:)))
    error('No good quality labels in HCTSA_loc.mat');
end

% ------------------------------------------------------------------------------
% Get handles for figure (f) and axes (ax):
f = figure('color','w');
ax = gca;

switch inspectWhat
case 'full'
    % Plot the quality labels:
    imagesc(TS_Quality)
    
    xlabel('Operations (op_id)','interpreter','none')
    ax.XTick = 1:length(Operations);
    ax.XTickLabel = [Operations.ID];

    formatYAxisColorBar;
    
case 'reduced'
    % First find where problems exist, and only show these columns
    qualityMean = mean(TS_Quality);
    hadProblem = (qualityMean > 0);
    imagesc(TS_Quality(:,hadProblem));
    
    ax.XTick = 1:sum(hadProblem);
    ax.XTickLabel = [Operations(hadProblem).ID];
    
    title(sprintf('Displaying %u x %u (keeping %u/%u operations with some special values)',...
                    size(TS_Quality,1),sum(hadProblem),sum(hadProblem),size(TS_Quality,2)),...
                    'interpreter','none')
                    
    formatYAxisColorBar;
    
case 'master'
    % Summarize at the level of master operations
    % Each row is now a different quality label
    
    % First get 
    numMasters = length(MasterOperations);
    opMIDs = [Operations.MasterID]; % save time by getting this first
    qualities = zeros(7,numMasters);
    for k = 1:numMasters
        masterID = [MasterOperations(k).ID];
        pointInd = (opMIDs==masterID); % indices of operations using this master function
        qualityLabels = TS_Quality(:,pointInd);
        qualityLabels = qualityLabels(qualityLabels>0); % only want to look at special ones
        if isempty(qualityLabels) % no problems
            continue
        else
            % some problems -- work out which
            qualityLabels = unique(qualityLabels);
            qualities(:,k) = ismember(1:7,qualityLabels);
        end
    end
    
    % Now keep only master operations with problems:
    hadProblem = (mean(qualities)>0)
    
    imagesc(qualities(:,hadProblem));
    
    colormap([0.1686,0.5137,0.7294; 0.8431,0.0980,0.1098]);
    
    ax.YTick = 1:8;
    ax.YTickLabel = {'error','NaN','Inf','-Inf','complex','empty','link error'};
    
    ax.XTick = 1:sum(hadProblem);
    ax.XTickLabel = {MasterOperations(hadProblem).Label};

    % Get rid of tex interpreter format (for strings with underscores)
    set(ax,'TickLabelInterpreter','none');
    
    title('Master operations producing special-valued outputs.','interpreter','none')
end


function formatYAxisColorBar()
    % For visualizing where you have time series on y-axis and color shows what quality
    % label for each computation
    
    % Format the y axis
    ax.YTick = 1:length(TimeSeries);
    ylabel('Time series')
    ax.YTickLabel = {TimeSeries.FileName};

    % Get rid of tex interpreter format (for strings with underscores)
    set(ax,'TickLabelInterpreter','none');
    
    % Add a color bar:
    allLabels = {'good','error','NaN','Inf','-Inf','complex','empty','link error'};
    labelsExist = unique(TS_Quality(~isnan(TS_Quality)));
    maxShow = max(labelsExist); % maximum quality label to plot

    cb = colorbar(ax);
    caxis([-0.5,maxShow+0.5]);
    cb.Ticks = [0:1:maxShow];
    cb.TickLabels = allLabels(1:maxShow+1);

    % Set the colormap:
    allColors = BF_getcmap('spectral',8,0);
    allColors = allColors([8,1,2,3,4,5,6,7],:);
    colormap(allColors(1:maxShow+1,:))
end

end