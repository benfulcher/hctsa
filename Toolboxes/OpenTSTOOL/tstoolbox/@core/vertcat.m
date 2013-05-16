function r = vertcat(c1,c2)

%tstoolbox/@core/vertcat
%   Syntax:
%     * r = vertcat(c1,c2)
%
%   Input Arguments:
%     * c1,c2 - core objects
%
%   catenate two timeseries verticaly
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


r = core([data(c1) ; data(c2)]);
