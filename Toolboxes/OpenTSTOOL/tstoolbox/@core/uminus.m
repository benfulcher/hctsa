function r = uminus(c)

%tstoolbox/@core/uminus
%   Syntax:
%     * r = uminus(c)
%
%   Input Arguments:
%     * c - core object
%
%   negate time series
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

r = core(-data(c));
