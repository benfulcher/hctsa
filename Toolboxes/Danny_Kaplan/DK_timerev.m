function res = DK_timerev(x,timeLag)
% DK_timerev Time reversal asymmetry statistic
%
% Calculates a time reversal asymmetry statistic
% x -- the time series
% timeLag -- a time scale (in samples) default 1
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

if nargin < 2
    timeLag = 1;
end

foo = DK_lagembed(x,3,timeLag);
a = foo(:,1);
b = foo(:,2);
c = foo(:,3);

res = mean(a.*a.*b - b.*c.*c);

end
