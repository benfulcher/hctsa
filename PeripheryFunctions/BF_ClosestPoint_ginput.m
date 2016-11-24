function iPlot = BF_ClosestPoint_ginput(xy,inputPoint)
% ClosestPoint_ginput   The closest point in a dataset to the input co-ordinates given.
%
%---INPUTS:
% xy is a Nx2 vector of x-y co-ordinates
% inputPoint is the output of a ginput
%
%---OUTPUT:
% plotMe, the closest point.

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
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(inputPoint)
    inputPoint = ginput(1);
end

%-------------------------------------------------------------------------------
% Calculate distances from each point to input inputPoint
%-------------------------------------------------------------------------------
dpxy = sum((xy - repmat(inputPoint,size(xy,1),1)).^2,2); % Euclidean distances to the input point
[~,iPlot] = min(dpxy);

end
