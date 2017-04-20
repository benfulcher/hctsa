function Q = DK_theilerQ(x)
% DK_theilerQ Computes Theiler's Q statistic
%
% theilerQ calculates Q=<(x_t + x_{t+1})^3> normalized by <x^2>^{3/2}
% on a vector x
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

x2 = mean(x.^2)^(3/2);

if x2 == 0
	Q = 0;
else
	d2 = x(1:end-1) + x(2:end);
	Q = mean(d2.^ 3)/x2;
end

end
