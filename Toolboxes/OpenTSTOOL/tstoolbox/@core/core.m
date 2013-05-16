function c = core(arg)

%tstoolbox/@core/core
%   core
%   class constructor Syntax:
%     * c = core(arg)
%
%   Input Arguments:
%     * arg double array
%
%   A core object contains the pure data part of a signal object.
%   Methods: ndim dlens data
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if nargin < 1	
	help(mfilename)
	return
end

if ~isa(arg, 'double')
	help(mfilename)
	return
end

c.data = arg;	 
c.dlens = dims(arg);

c = class(c, 'core');

