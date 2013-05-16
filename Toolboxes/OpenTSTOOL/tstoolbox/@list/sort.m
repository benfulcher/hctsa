function [l, index] = sort(l)

%tstoolbox/@list/sort
%   sort list l in increasing order.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

temp = char(l.data);         % convert cell-array to fixed-length array, each string in a row
[temp, index] = sortrows(temp);
l.data = (cellstr(temp))';
