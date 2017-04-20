function out = DK_crinkle(x)
% DK_crinkle  Computes James Theiler's crinkle statistic
%
% Calculates the "crinkle statistic" on a vector x
% 	<(x_{t-1}-2*x_t+x_{t+1})^4> / < ( x_t^2 ) >^2
%	as proposed by James Theiler
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
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% Subtract out the mean
x = x - mean(x);
x2 = mean(x.*x)^2;

if x2 == 0
	out = 0;
	return
end

d2 = 2*x(2:end-1) - x(1:end-2) - x(3:end);
out = mean(d2.^ 4)/x2;

end
