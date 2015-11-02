function [fr,xr,h_line,h_points] = BF_plot_ks(dataVector,whatColor,doSwap,lineWidth,markerSize)
% Plot a kernel smoothed distribution with individual datapoints annotated
%
%---INPUTS:
% dataVector, a vector of data to plot

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

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2, whatColor = 'k'; end
if nargin < 3, doSwap = 0; end
if nargin < 4, lineWidth = 1; end
if nargin < 5, markerSize = 8; end

numPoints = 1000; % points for the ks density

%-------------------------------------------------------------------------------
% Estimate a probability density:
%-------------------------------------------------------------------------------
[f,x] = ksdensity(dataVector,linspace(min(dataVector),max(dataVector),numPoints),...
                    'function','pdf');

if all(dataVector==dataVector(1))
    warning('Attempting to estimate a kernel-smoothed probability distribution using constant data')
end

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
