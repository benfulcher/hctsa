%tstoolbox/mex/henon
%   Generate time series by iterating the henon map.
%
%   Syntax:
%
%     * x = henon(length, [a b xo yo])
%
%   Input arguments:
%
%     * length - number of samples to generate
%     * [a b xo yo] - vector of parameters and initial conditions
%
%   Output arguments:
%
%     * x - vector of size D
%
%   Example:
%
%x =  henon(500, [-1.4 0.3 0.2 0.12]);
%plot(x(:,1), x(:,2), '.');
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


