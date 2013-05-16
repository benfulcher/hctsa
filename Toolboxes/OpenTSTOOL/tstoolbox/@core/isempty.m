function r = isempty(s)

%tstoolbox/@core/isempty
%   Syntax:
%     * r = isempty(s)
%
%   Input Arguments:
%     * s - core object
%
%   test if core contains no (valid) data
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

r = isempty(s.data);
