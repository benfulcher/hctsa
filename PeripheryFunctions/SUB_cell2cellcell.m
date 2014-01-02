% SUB_cell2cellcell
% 
% Inputs a cell with some delimiter, and outputs a cell of cells using this
% delimiter.
% 
% Ben Fulcher 24/3/2010
% 

function cellcell = SUB_cell2cellcell(cellin,delimiter)

if nargin < 2 || isempty(delimiter)
    delimiter = ',';
end

Nelements = length(cellin); % number of elements in the input cell
cellcell = cell(Nelements,1);

for i = 1:Nelements
    cellcell{i} = regexp(cellin{i}, delimiter, 'split');
end


end