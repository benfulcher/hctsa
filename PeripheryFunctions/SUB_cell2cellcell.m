function outputCellCell = SUB_cell2cellcell(inputCell,delimiter)
% SUB_cell2cellcell     Turn a cell of delimited strings into a cell of cells.
%
% Inputs a cell with some delimiter, and outputs a cell of cells using this
% delimiter.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% Comma delimited by default:
if nargin < 2 || isempty(delimiter)
    delimiter = ',';
end

% ------------------------------------------------------------------------------

numElements = length(inputCell); % number of elements in the input cell
outputCellCell = cell(numElements,1);

for i = 1:numElements
    outputCellCell{i} = regexp(inputCell{i}, delimiter, 'split');
end

end
