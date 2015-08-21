function BF_PlotCorrMat(D_corr,rangeHow)
% BF_PlotCorrMat    Visualization of a pairwise similarity matrix
%
% Attempts to determine a set of smaller clusters of objects showing similar
% behavior.
%
%---INPUTS:
% D_corr, a pairwise distance matrix
% rangeHow, the colorbar extent: (i) '' (automatic)
%                                (ii) '-1to1' (range from -1 to 1)
%                                (iii) '0to1' (range from 0 to 1)

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

if nargin < 2 || isempty(rangeHow)
    rangeHow = ''; %'-1to1': range from -1 to 1, or '0to1' for unit interval
end
%-------------------------------------------------------------------------------

% Make a matrix
if any(size(D_corr)==1)
    D_corr = squareform(D_corr);
end

imagesc(D_corr)

axis square

% Set color limits and colormap
switch rangeHow
case '-1to1'
    caxis([-1,1])
    colormap([flipud(BF_getcmap('blues',9,0));BF_getcmap('reds',9,0)])
case '0to1' % assume [0,1] (a normalized distance metric)
    caxis([0,1])
    colormap(BF_getcmap('reds',9,0))
otherwise
    colormap(flipud(BF_getcmap('reds',9,0)))
end

% ------------------------------------------------------------------------------
% Superimpose green/yellow rectangles over NaN values
% ------------------------------------------------------------------------------
if any(isnan(D_corr(:)))
    [theNaNs_i,theNaNs_j] = find(isnan(D_corr));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],'FaceColor','k', ...
                        'EdgeColor','k')
    end
end

% Remove ticks:
set(gca,'TickLength',[0,0])

end
