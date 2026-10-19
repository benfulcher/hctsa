function BF_PlotCorrMat(D_corr,rangeHow,makeFigure)
% BF_PlotCorrMat    Visualization of a pairwise similarity matrix
%
% Plots a square matrix (e.g., of correlation coefficients) as an image, with
% a colormap chosen by the range of values, and NaNs drawn in black.
%
%---INPUTS:
% D_corr, a square matrix (or its vector form, cf. squareform) of pairwise
%               similarities or correlations
% rangeHow, the colorbar extent: (i) 'auto' (automatic: grayscale, or 'balanced' when
%                                    there are both positive and negative values)
%                                (ii) 'n1to1' (range from -1 to 1; '-1to1' is
%                                    accepted as a synonym)
%                                (iii) '0to1' (range from 0 to 1)
%                                (iv) 'balanced' (range from -X to X, for max possible X)
%                                (v) 'positive' (range from 0 to X, for max possible X)
% makeFigure, whether to generate a new figure

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check inputs, set defaults:
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(rangeHow)
    rangeHow = 'auto';
end
if nargin < 3
    makeFigure = false;
end
if strcmp(rangeHow,'-1to1')
    rangeHow = 'n1to1';
end

%-------------------------------------------------------------------------------
% Make a matrix
if any(size(D_corr)==1)
    D_corr = squareform(D_corr);
end

% Open a new figure window:
if makeFigure
    figure('color','w');
end

% Plot the image:
imagesc(D_corr)

axis square

%-------------------------------------------------------------------------------
% Set color limits and colormap

% If auto setting, can set to balanced if both positive and negative
if strcmp(rangeHow,'auto')
    if any(D_corr(:) > 0) && any(D_corr(:) < 0)
        rangeHow = 'balanced';
    end
end

switch rangeHow
case 'positive'
    maxDev = max(D_corr(:));
    clim([0,maxDev])
    colormap(BF_GetColorMap('greenblue',9,0))
case 'balanced'
    maxDev = max(abs(D_corr(:)));
    clim([-maxDev,maxDev])
    colormap([flipud(BF_GetColorMap('blues',9,0));[1,1,1],;BF_GetColorMap('reds',9,0)])
case 'n1to1'
    clim([-1,1])
    colormap([flipud(BF_GetColorMap('blues',9,0));BF_GetColorMap('reds',9,0)])
case '0to1' % assume 0-1 range (a normalized distance metric)
    clim([0,1])
    colormap([1,1,1;BF_GetColorMap('reds',9,0)])
case 'auto'
    colormap(flipud(gray(64)))
otherwise
    error('Unknown colormap range option: ''%s''',rangeHow);
end

% ------------------------------------------------------------------------------
% Superimpose black over NaN values (a single black image that is transparent
% everywhere except at the NaNs, rather than a rectangle for each NaN)

if any(isnan(D_corr(:)))
    % (a low-level image call, which adds to the axes rather than replacing the matrix)
    image('Parent',gca,'XData',[1,size(D_corr,2)],'YData',[1,size(D_corr,1)], ...
            'CData',zeros([size(D_corr),3]),'AlphaData',double(isnan(D_corr)))
end

% Remove ticks:
set(gca,'TickLength',[0,0])

end
