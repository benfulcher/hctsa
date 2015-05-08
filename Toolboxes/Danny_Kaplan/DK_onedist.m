% DK_onedist(z,pt) calculates the distance between point pt and each
% row in matrix z
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
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function res = DK_onedist(z,pt)

[r,c] = size(z);
if c ~= length(pt)
  error('pt and z must have same number of columns');
end

sum = zeros(r,1);
for n=1:c
  foo = z(:,n) - pt(n);
  sum = sum + foo.*foo;
end
res = sqrt(sum);
