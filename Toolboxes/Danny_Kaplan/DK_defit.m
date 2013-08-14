% DK_defit
% 
% [a,b] = DK_defit(delta,epsilon,maxdelta)
% 
% linear fitting routine for delta-epsilon
% delta   -- distances between pre-images: output by delta-epsilon
% epsilon -- distances between images: output by delta-epsilon
% maxdelta-- optional - largest delta to consider.
% 
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
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function [a, b] = DK_defit(delta,epsilon,maxdelta)

if nargin < 3
  maxdelta = 10e100;
end

% get rid of points whose pre-images are too far apart.
foo = (delta < maxdelta);
delta = delta(foo);
epsilon = epsilon(foo);

% Sort the pairs into ascending order on delta
[x, i] = sort(delta);
y = epsilon(i);

s = length(x);
sx = sum(x);
sy = sum(y);
sxx = sum( x.*x);
sxy = sum( x.*y);
del = s*sxx - sx*sx;
a = (sxx*sy - sx*sxy)/del;
b = (s*sxy - sx*sy)/del;

end