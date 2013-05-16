%tstoolbox/mex/tentmap
%   Generate samples of the generalized iterated tentmap.
%
%   Syntax:
%
%     * x = tentmap(length, [h e s x0])
%
%   Input arguments:
%
%     * length - number of samples to generate
%     * [h e s x0] - vector of parameters and initial conditions
%
%   Output arguments:
%
%     * x - time series
%
%   Example:
%
%x = tentmap(500, [0 1 0.97 rand(1,1)]);
%plot(x)
%plot(x(1:end-1), x(2:end), '.')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


