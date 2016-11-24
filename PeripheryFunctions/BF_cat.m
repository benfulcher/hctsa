function out = BF_cat(s,d,surr)
% BF_cat    Concatenate strings.
%
% Converts a cell of strings (or vector of numbers) into a single, concatonated
% string listing each string in the cell.
%
%---INPUTS:
% s, the cell of strings to be concatinated
% d, the delimiter (a comma by default)
% surr, an optional string to surround each element of s (empty by default, but
%       can be handy for surrounding each string by inverted commas, for
%       example).

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    d = ', '; % default delimiter is a comma
else
    d = [d, ' '];
end

if nargin < 3 || isempty(surr)
    sumString = [];
    if iscellstr(s)
        for i = 1:length(s)
            sumString = [sumString, s{i}, d];
        end
    elseif isnumeric(s)
        for i = 1:length(s)
            sumString = [sumString, num2str(s(i)), d];
        end
    end
else
    sumString = [];
    if iscellstr(s)
        for i = 1:length(s)
            sumString = [sumString, surr, s{i}, surr, d];
        end
    elseif isnumeric(s)
        for i = 1:length(s)
            sumString = [sumString, surr, num2str(s(i)), surr, d];
        end
    end
end

out = sumString(1:end-length(d));

end
