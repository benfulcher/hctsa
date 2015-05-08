% DK_lagembed(x,dim,lag) constructs an embedding of a time series on a vector
% DK_lagembed(x,dim) makes an m-dimensional embedding with lag 1
% DK_lagembed(x,dim,lag) uses the specified lag
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

function y = DK_lagembed(x,M,lag)

if nargin < 3
	lag = 1;
end
if nargin < 4
	advance=0;
end

%convert x to a column
[xr,xc] = size(x);
if xr == 1	
    x = x';
end


lx = length(x);
	
newsize = lx - lag*(M-1);
y = zeros(newsize,M);
i=1;

for j = 0:-lag:lag*(-(M-1))

	first=1+lag*(M-1)+j;

	last=first+newsize-1;


	if last > lx

		last = lx;

	end

	y(:,i) = x(first:last, 1);

	i = i+1;

end

