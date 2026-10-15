function colorCell = GiveMeColors(numColors)
% GiveMeColors outputs a cell of colors for the specified number of colors
%
% Relies on the BF_GetColorMap function to pull out the appropriate colorbrewer
% colormaps (qualitative palettes for up to 21 colors, then the jet colormap).
%
%---INPUT:
% numColors, the number of colors (groups) required
%
%---OUTPUT:
% colorCell, a numColors x 1 cell of RGB triplets

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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
    colorCell = BF_GetColorMap('set2',numColors,true);
elseif numColors <= 5
    colorCell = BF_GetColorMap('set1',numColors,true);
elseif numColors <= 8
    colorCell = BF_GetColorMap('dark2',numColors,true); % dark2 has 8 colors
elseif numColors == 9
    colorCell = BF_GetColorMap('set1',numColors,true); % set1 has 9 colors
elseif numColors <= 12
    colorCell = BF_GetColorMap('set3',numColors,true); % set3 has 12 colors
elseif numColors <= 21
    % Combine set1 (up to 9 colors) and set3 (up to 12 colors):
    numSet1 = min(ceil(numColors/2),9);
    colorCell = [BF_GetColorMap('set1',numSet1,true);...
                 BF_GetColorMap('set3',numColors-numSet1,true)];
else
    % Too many for distinct qualitative colors: sample the jet colormap
    colorCell = num2cell(jet(numColors),2);
end

end
