function [fr,xr,h_line,h_points] = BF_plot_ks(dataVector,whatColor,doSwap,lineWidth,markerSize,~)
% Plot a kernel smoothed distribution with individual datapoints annotated
%
%---INPUTS:
% dataVector, a vector of data to plot (NaNs are not plotted)
% whatColor, the color to plot line and points
% doSwap, (logical) whether to swap axes to be horizontal
% lineWidth, width of line to plot
% markerSize, size of markers annotating points
% (sixth input), trimRange, accepted for backwards compatibility but has no effect: the density
%           is only ever evaluated between the minimum and maximum of the data
%
%---OUTPUTS:
% fr, the height of the density at each datapoint, in the same order as dataVector
%       (a row vector, NaN where dataVector is NaN)
% xr, the position of each datapoint on the density's grid, in the same order
%       as dataVector (so that [xr(i),fr(i)] is the plotted point for dataVector(i))
% h_line, handle to the plotted density curve
% h_points, handle to the plotted datapoints
%
% For constant data, the density is not defined, so a single stem (of unit
% height) at the value is plotted.
%
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
% Check inputs and set defaults:
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(whatColor)
    whatColor = 'k';
end
if nargin < 3 || isempty(doSwap)
    doSwap = false;
end
if nargin < 4 || isempty(lineWidth)
    lineWidth = 1;
end
if nargin < 5 || isempty(markerSize)
    markerSize = 14;
end

numPoints = 1000; % points for the ks density

%-------------------------------------------------------------------------------
% Estimate a probability density:
%-------------------------------------------------------------------------------
isGood = ~isnan(dataVector(:)');
goodData = dataVector(isGood);
if isempty(goodData)
    error('No (non-NaN) data to plot');
end

fr = nan(1,length(dataVector));
xr = nan(1,length(dataVector));

if max(goodData) - min(goodData) < 10*eps(max(abs(goodData))) % ~ constant data
    if length(goodData) > 1
        warning('Constant data: plotting a single stem')
    end
    x = goodData(1)*ones(1,2);
    f = [0,1];
    xr(isGood) = goodData(1);
    fr(isGood) = 1;
else
    [f,x] = ksdensity(goodData,linspace(min(goodData),max(goodData),numPoints),...
                        'function','pdf');

    %---------------------------------------------------------------------------
    % Match each datapoint to a point on the distribution (the first grid point
    % at or above it), keeping the order of dataVector:
    %---------------------------------------------------------------------------
    ind = discretize(goodData,[-Inf,x],'IncludedEdge','right');
    xr(isGood) = x(ind);
    fr(isGood) = f(ind);
end

%-------------------------------------------------------------------------------
% Plot the line and matched points:
%-------------------------------------------------------------------------------
if ~doSwap
    h_line = plot(x,f,'color',whatColor,'LineWidth',lineWidth); % the curve
    h_points = plot(xr,fr,'.','color',whatColor,'MarkerSize',markerSize); % individual TimeSeries as points
else % (plot horizontally rather than vertically)
    h_line = plot(f,x,'color',whatColor,'LineWidth',lineWidth); % the curve
    h_points = plot(fr,xr,'.','color',whatColor,'MarkerSize',markerSize); % individual TimeSeries as points
end

end
