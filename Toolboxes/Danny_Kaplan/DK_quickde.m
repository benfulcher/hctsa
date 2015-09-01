function de = DK_quickde(ts,dim,lag,nmin)
% DK_quickde   Quick-and-dirty characterization of determinism using de
%
%---INPUTS:
% ts -- the time series
% dim -- the embedding dimension
% lag -- the embedding lag
% nmin -- optional: number of points to use for delta-eps fitting
%           default value: 500
%
%
% Tweaked by B. D. Fulcher
% ------------------------------------------------------------------------------
% Copyright (C) 1996, D. Kaplan <kaplan@macalester.edu>
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(dim)
    dim = 2;
end
if nargin < 3 || isempty(lag)
    lag = 1;
end
if strcmp(lag,'ac')
    lag = CO_FirstZero(y,'ac');
end

if nargin < 4
    nmin = 500;
end

%-------------------------------------------------------------------------------

% Embed the data
xx = DK_lagembed(ts,dim,lag);

% Get the 'forecasting' pre-image and image
[pre, post] = DK_getimage(xx,1);

% Find a 'typical' distance to use for delta-epsilon
% We want to have at least nmin points, so we pick

perc = (2*nmin)/(length(ts).^2);
deltamax = DK_disttyp(pre,perc);
[delta, epsilon] = DK_deltaeps(pre,post);
[de, b] = DK_defit(delta,epsilon,deltamax);

end
