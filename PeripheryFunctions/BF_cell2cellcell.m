function cellcell = BF_cell2cellcell(cellin,delimiter)
% Inputs a cell with some delimiter, and outputs a cell of cells using this
% delimiter.
% Ben Fulcher 24/3/2010

if nargin < 2 || isempty(delimiter)
    delimiter = ','; % comma as the default delimiter
end

Nelements = length(cellin); % number of elements in the input cell
cellcell = cell(Nelements,1);

for i = 1:Nelements
    cellcell{i} = regexp(cellin{i}, delimiter, 'split');
end

end