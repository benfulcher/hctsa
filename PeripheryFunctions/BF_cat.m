% BF_cat
% 
% Converts a cell of strings (or vector of numbers) into a single, concatonated
% string listing each string in the cell.
% 
% INPUTS:
% s, the cell of strings to be concatinated
% d, the delimiter (a comma by default)
% surr, an optional string to surround each element of s (empty by default, but
%       can be handy for surrounding each string by inverted commas, for
%       example)
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
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

function out = BF_cat(s,d,surr)

if nargin < 2
    d = ', '; % default delimiter is a comma
else
    d = [d, ' '];
end

if nargin < 3 || isempty(surr)
    sumstr = [];
    if iscellstr(s)
        for i = 1:length(s), sumstr = [sumstr, s{i}, d]; end
    elseif isnumeric(s)
        for i = 1:length(s), sumstr = [sumstr, num2str(s(i)), d]; end
    end
else
    sumstr = [];
    if iscellstr(s)
        for i = 1:length(s), sumstr = [sumstr, surr, s{i}, surr, d]; end
    elseif isnumeric(s)
        for i = 1:length(s), sumstr = [sumstr, surr, num2str(s(i)), surr, d]; end
    end    
end

out = sumstr(1:end-length(d));

end