% ------------------------------------------------------------------------------
% ClosestPoint_ginput
% ------------------------------------------------------------------------------
% Returns the closest InputPoint to the one given
% Ben Fulcher 20/10/2010
%
%---INPUTS:
% xy can be a cell, each component of which is a different group plotted,
%                       contains a Nx2 vector of co-ordinates
% InputPoint is the output of a ginput
% 
%---OUTPUT:
% PlotMe, the closest point.
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

function PlotMe = ClosestPoint_ginput(xy,InputPoint)

if nargin < 2 || isempty(InputPoint)
    InputPoint = ginput(1);
end

if iscell(xy)
    NumGroups = length(xy);
else
    NumGroups = 1;
end

% Calculate distances from each point to input InputPoint
dpxy = cell(NumGroups,1);
for i = 1:NumGroups
    dpxy{i} = sum((xy{i} - repmat(InputPoint,size(xy{i},1),1)).^2,2); % Euclidean distances to the input point
end
if iscell(xy)
    mins = cellfun(@(x)min(x),dpxy);
    thisgroup = find(mins==min(mins),1); % This group contains closest element to InputPoint
else
    thisgroup = 1;
end

PlotMe = find(dpxy{thisgroup} == min(mins),1); % this is the point to plot
if iscell(xy)
    PlotMe = [thisgroup,PlotMe]; % Also output the group
end

end