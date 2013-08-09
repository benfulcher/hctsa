% BF_cell2cellcell
% 
% Inputs a cell with some delimiter, and outputs a cell of cells using this
% delimiter.
% 
% INPUTS:
% cellin, the cell
% delimiter, the delimiter
% 
% OUTPUT:
% cellcell, the cell of cells, using this delimiter
% 
% Used for some tasks involving mySQL.
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

function cellcell = BF_cell2cellcell(cellin,delimiter)
% Ben Fulcher, 24/3/2010

if nargin < 2 || isempty(delimiter)
    delimiter = ','; % comma as the default delimiter
end

Nelements = length(cellin); % number of elements in the input cell
cellcell = cell(Nelements,1);

for i = 1:Nelements
    cellcell{i} = regexp(cellin{i}, delimiter, 'split');
end

end