% DK_quickde
% 
% Does a quick-and-dirty characterization of determinism using de
% ts -- the time series
% dim -- the embedding dimension
% lag -- the embedding lag
% nmin -- optional: number of points to use for delta-eps fitting
%           default value: 500
% 
% 
% Tweaked ever so slightly by B. D. Fulcher
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

function de = DK_quickde(ts,dim,lag,nmin)

if nargin < 4
    nmin = 500;
end

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