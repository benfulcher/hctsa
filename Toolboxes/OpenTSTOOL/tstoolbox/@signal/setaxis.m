function s = setaxis(s, dim, a)

%tstoolbox/@signal/setaxis
%   Syntax:
%     * s = setaxis(s, dim, achse)
%
%   Change one of the current xaxes.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

% change one of the current xaxes
%
% s = setaxis(s, dim, achse)	% set axis number dim to new value achse
%
% C.Merkwirth,U.Parlitz,W.Lauterborn  DPI Goettingen 1998

narginchk(3,3);

s.xaxes{dim} = a;	
