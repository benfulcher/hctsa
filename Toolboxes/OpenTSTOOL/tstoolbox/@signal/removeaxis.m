function s = removeaxis(s, dim)

%tstoolbox/@signal/removeaxis
%   Syntax:
%     * s = removeaxis(s, dim)
%
%   Remove axis one of the current xaxes. No bound checking for dim.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);


% no bound checking for dim

s.xaxes(dim) = [];	
