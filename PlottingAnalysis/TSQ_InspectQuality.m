% --------------------------------------------------------------------------
% TSQ_InspectQuality
% --------------------------------------------------------------------------
% 
% This function loads the matrix in HCTSA_loc.mat, plots it, showing the 
% quality labels of each entry.
% 
% Most useful for checking where errors/special-valued outputs are occuring
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
function TSQ_InspectQuality()

load('HCTSA_loc.mat','TS_Quality','TimeSeries','Operations');

maxShow = max(TS_Quality(~isnan(TS_Quality))); % maximum quality label to plot
if isempty(maxShow)
    error('No good quality labels in HCTSA_loc.mat');
end

f = figure('color','w'); box('on');
imagesc(TS_Quality)
ax = gca;

ylabel('Time series')
ax.YTick = 1:length(TimeSeries);
ax.YTickLabel = {TimeSeries.FileName};
xlabel('Operations (op_id)','interpreter','none')
ax.XTick = 1:length(Operations);
ax.XTickLabel = [Operations.ID];

% Add a color bar:
allLabels = {'good','error','NaN','Inf','-Inf','complex','empty','link error'};
cb = colorbar(ax);
caxis([-0.5,maxShow+0.5]);
cb.Ticks = [0:1:maxShow];
cb.TickLabels = allLabels(1:maxShow+1);
allColors = BF_getcmap('spectral',8,0);
allColors = allColors([8,1,2,3,4,5,6,7],:);
colormap(allColors(1:maxShow+1,:))

end