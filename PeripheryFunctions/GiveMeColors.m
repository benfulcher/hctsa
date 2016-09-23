function colorCell = GiveMeColors(numColors)
% GiveMeColors outputs a cell of colors for the specified number of colors
%
% Relies on the BF_getcmap function to pull out the appropriate colorbrewer
% colormaps

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if numColors == 1
    colorCell = {[0,0,0]}; % Just use black...
elseif numColors==4
    colorCell = BF_getcmap('set2',numColors,1); % set1, accent, set2
elseif numColors <= 5
    % colorCell = BF_getcmap('set1',5,1);
    colorCell = BF_getcmap('set1',numColors,1); % set1, accent, set2
    % colorCell = BF_getcmap('set2',5,1);
    % if numColors==2, colorCell = colorCell([2,4]); end
elseif numColors < 10
    colorCell = BF_getcmap('dark2',numColors,1);
elseif numColors <= 12
    colorCell = BF_getcmap('set3',numColors,1);
elseif numColors <= 22
    colorCell = [BF_getcmap('set1',numColors,1); ...
                BF_getcmap('set3',numColors,1)];
elseif numColors <= 50
    colorCell = mat2cell(jet(numColors));
else
    error('There aren''t enough colors in the rainbow to plot this many groups!')
end

end
