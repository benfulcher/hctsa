%tstoolbox/mex/baker
%   Generate time-series from the iterated Baker map .
%
%   Syntax:
%
%     * x = baker(length, [eta l1 l2 x0 y0])
%
%   Input arguments:
%
%     * length - number of samples to generate
%     * [eta l1 l2 x0 y0] - vector of parameters and initial conditions
%
%   Output arguments:
%
%     * x - time series
%
%   Example:
%
%x = baker(2000, [0.6 0.25 0.4 rand(1,1) rand(1,1)]);
%plot(x(1:end-1,2), x(2:end,2), '.')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

