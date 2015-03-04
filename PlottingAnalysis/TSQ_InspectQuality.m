
load('HCTSA_loc.mat','TS_Quality','TimeSeries','Operations');
f = figure('color','w'); box('on');
imagesc(TS_Quality)
ax = gca;

ylabel('Time series')
ax.YTick = 1:length(TimeSeries);
ax.YTickLabel = {TimeSeries.FileName};
xlabel('Operations')

% Add a color bar:
allLabels = {'good','error','NaN','Inf','-Inf','complex'};
maxShow = max(TS_Quality(:));
cb = colorbar(ax);
caxis([-0.5,maxShow+0.5]);
cb.Ticks = [0:1:maxShow];
cb.TickLabels = allLabels(1:maxShow+1);
allColors = BF_getcmap('spectral',6,0);
allColors = allColors([6,1,2,3,4,5],:);
colormap(allColors(1:maxShow+1,:))