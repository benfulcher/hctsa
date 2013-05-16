function [result, index] = cellsort(x)

% tstool/cellsort
% sort an cell-array x of strings
% result is an cell array containing the strings in ascending order
% index contains the permutation order

if nargin < 1, help(mfilename), return, end 

if iscellstr(x)
	temp = char(x); 	% convert cell-array to fixed-length array, each string in a row
	[temp, index] = sortrows(temp);	
	result = cellstr(temp);	% convert array back to cell-array
else
	error('Argument is not a cell array of strings');
end
