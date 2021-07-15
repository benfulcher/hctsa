function [iPlot,minDist] = BF_ClosestPoint_ginput(xy,inputPoint,doNormalize)
% ClosestPoint_ginput   The closest point in a dataset to the input co-ordinates given.
%
%---INPUTS:
% xy, A Nx2 vector of x-y co-ordinates.
% inputPoint, The output of a ginput.
% doNormalize, whether to normalize the input space so that relative differences
%               in x and y are treated (rather than raw Euclidean distances).
%
%---OUTPUT:
% plotMe, the closest point.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(inputPoint)
    inputPoint = ginput(1);
end
if nargin < 3
    doNormalize = true;
end

%-------------------------------------------------------------------------------
% Normalize
%-------------------------------------------------------------------------------
% Normalizing the space to a unit square treats distances in both axes similarly
% relative to their ranges.
if doNormalize
    minMaxx = [min(xy(:,1)),max(xy(:,1))];
    minMaxy = [min(xy(:,2)),max(xy(:,2))];
    xy(:,1) = (xy(:,1) - minMaxx(1))/(minMaxx(2)-minMaxx(1));
    inputPoint(1) = (inputPoint(1) - minMaxx(1))/(minMaxx(2)-minMaxx(1));
    xy(:,2) = (xy(:,2) - minMaxy(1))/(minMaxy(2)-minMaxy(1));
    inputPoint(2) = (inputPoint(2) - minMaxy(1))/(minMaxy(2)-minMaxy(1));
end

%-------------------------------------------------------------------------------
% Calculate distances from each point to input inputPoint
%-------------------------------------------------------------------------------
dpxy = sum((xy - repmat(inputPoint,size(xy,1),1)).^2,2); % Euclidean distances to the input point
[minDist,iPlot] = min(dpxy);

end
