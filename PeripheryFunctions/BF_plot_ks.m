function [fr,xr,h_line,h_points] = BF_plot_ks(dataVector,whatColor,doSwap,lineWidth,markerSize,trimRange)
% Plot a kernel smoothed distribution with individual datapoints annotated
%
%---INPUTS:
% dataVector, a vector of data to plot
% whatColor, the color to plot line and points
% doSwap, (logical) whether to swap axes to be horizontal
% lineWidth, width of line to plot
% markerSize, size of markers annotating points
% trimRange, (logical) whether to trim the range of plot to that of the data points
%
% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check inputs and set defaults:
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(whatColor)
    whatColor = 'k';
end
if nargin < 3 || isempty(doSwap)
    doSwap = 0;
end
if nargin < 4 || isempty(lineWidth)
    lineWidth = 1;
end
if nargin < 5 || isempty(markerSize)
    markerSize = 14;
end
if nargin < 6 || isempty(trimRange)
    trimRange = 0;
end

numPoints = 1000; % points for the ks density

%-------------------------------------------------------------------------------
% Estimate a probability density:
%-------------------------------------------------------------------------------
if all(dataVector-dataVector(1)<10*eps); % ~ constant data vector
    warning('Constant data')
    x = dataVector(1);
    f = 1;
end

[f,x] = ksdensity(dataVector,linspace(min(dataVector),max(dataVector),numPoints),...
                    'function','pdf');

%-------------------------------------------------------------------------------
% Match each datapoint to a point on the distribution:
%-------------------------------------------------------------------------------
getIndex = @(m) find(x>=dataVector(m),1,'first');
ind = arrayfun(getIndex,1:length(dataVector));
ind = sort(ind,'ascend');

% The matched points, [xr,fr]
xr = x(ind);
fr = f(ind);

%-------------------------------------------------------------------------------
% Trim
%-------------------------------------------------------------------------------
if trimRange
    rPlot = (x>=min(dataVector) & x<=max(dataVector));
    if sum(rPlot) > 0
        x = x(rPlot);
        f = f(rPlot);
    else
        x = ones(2,1)*xr(1);
        f = [0,fr(1)];
    end
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
